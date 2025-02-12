import re
import os
import time
import requests
import pandas as pd
import datetime
import gzip
import orjson
import multiprocessing
import argparse

import xml.etree.ElementTree as ET

from random import shuffle
from pathlib import Path


tag_pattern = '{https://uniprot.org/uniprot}'
# rest_uniprot_pattern_compressed = 'https://rest.uniprot.org/uniprotkb/stream?compressed=true&format=fasta&query=%28%28proteome%3A'
# rest_uniprot_pattern_uncompressed = 'https://rest.uniprot.org/uniprotkb/stream?compressed=false&format=fasta&query=%28%28proteome%3A'

########################################################
parser = argparse.ArgumentParser()
parser.add_argument(
    '--proteomes_table', type=str,
    help='reference proteomes table obtained from UniProt')
parser.add_argument(
    '--out_dir', type=str,
    help='output directory')
parser.add_argument(
    '--N_process', type=str,
    help='number of processes')
parser.add_argument(
    '--cache_size', type=str,
    help='number of taxids in each process'
)

FLAGS = parser.parse_args()
########################################################

N_processes = int(FLAGS.N_process)
suffix = str(datetime.datetime.now()).replace(':', '')
log_path = Path(FLAGS.out_dir) / f'fetch_annotation_log_{suffix}.txt'
database_log = Path(FLAGS.out_dir) / 'database_log.xlsx'
referece_proteomes_table_path = Path(FLAGS.proteomes_table)
proteome_table = pd.read_csv(referece_proteomes_table_path, sep = "\t")
proteome_table = proteome_table.set_index('Organism Id')

### schema: https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot.xsd

################## Functions to parse xml ##################

def parse_accessions(entry, index_dict, entry_dict):
    assert 'accession' in index_dict.keys()
    primary_accession = entry[index_dict['accession'][0]].text
    entry_dict['accession'] = dict()
    entry_dict['accession']['primary_accession'] = primary_accession
    entry_dict['accession']['all_accessions'] = []
    for index in index_dict['accession']:
        entry_dict['accession']['all_accessions'] += [entry[index].text]
    return entry_dict


def parse_names(entry, index_dict, entry_dict):
    assert 'name' in index_dict.keys()
    primary_name = entry[index_dict['name'][0]].text
    entry_dict['name'] = dict()
    entry_dict['name']['primary_name'] = primary_name
    entry_dict['name']['all_names'] = []
    for index in index_dict['name']:
        entry_dict['name']['all_names'] += [entry[index].text]
    return entry_dict

def parse_protein_gene(entry, index_dict, entry_dict):
    ## We only grab information from the 'recommendedName' item.
    tag_pattern = '{https://uniprot.org/uniprot}'
    assert 'protein' in index_dict.keys()
    protein_name = entry[index_dict['protein'][0]]
    assert protein_name[0].tag.replace(tag_pattern, '') in ['recommendedName', 'submittedName']
    assert protein_name[0][0].tag.replace(tag_pattern, '') == 'fullName'
    full_name = protein_name[0][0].text

    entry_dict['recommended_name'] = dict()
    entry_dict['recommended_name']['full_name'] = full_name
    entry_dict['recommended_name']['ec_numbers'] = []
    entry_dict['recommended_name']['short_names'] = []
    for subitem in protein_name[0]:
        if subitem.tag.replace(tag_pattern, '') == 'shortName':
            entry_dict['recommended_name']['short_names'] += [subitem.text]
        if subitem.tag.replace(tag_pattern, '') == 'ecNumber':
            entry_dict['recommended_name']['ec_numbers'] += [subitem.text]
    
    if 'gene' in index_dict.keys():
        ## We only grab information from the primary gene item.
        gene_dict = dict()
        gene_name = entry[index_dict['gene'][0]]
        for sub_item in gene_name:
            assert len(sub_item) == 0
            assert 'type' in sub_item.attrib.keys()
            if sub_item.attrib['type'] not in gene_name.keys():
                gene_dict[sub_item.attrib['type']] = []
            gene_dict[sub_item.attrib['type']] += [sub_item.text]
        
        entry_dict['gene_names'] = gene_dict
    return entry_dict

def extract_ref(ref_entry, ref_dict):
    db = ref_entry.attrib['type']
    if db not in ref_dict.keys():
        ref_dict[db] = dict()
    id_ = ref_entry.attrib['id']
    ref_dict[db][id_] = dict({'id': id_})
    for subitem in ref_entry:
        assert len(subitem) == 0
        if subitem.tag.replace(tag_pattern, '') == 'property':
            property_type = subitem.attrib['type']
            value = subitem.attrib['value']
        else:
            assert subitem.tag.replace(tag_pattern, '') == 'molecule'
            assert set(subitem.attrib.keys()) == {'id'}
            property_type = 'mocule_id'
            value = subitem.attrib['id']
        ref_dict[db][id_][property_type] = value
    return ref_dict

def parse_dbref(entry, index_dict, entry_dict):
    ref_dict = dict()
    ref_indeces = index_dict['dbReference']
    for index in ref_indeces:
        ref_dict = extract_ref(entry[index], ref_dict)
    entry_dict['db_references'] = ref_dict
    return entry_dict


def parse_keyword(entry, index_dict, entry_dict):
    keyword_dict = dict()
    kw_indeces = index_dict['keyword']
    for index in kw_indeces:
        assert 'id' in entry[index].attrib.keys()
        assert len(entry[index]) == 0
        keyword_dict[entry[index].attrib['id']] = entry[index].text
    entry_dict['keyword'] = keyword_dict
    return entry_dict


def parse_entry(entry):
    entry_tags = [item.tag.replace(tag_pattern, '') for item in entry]
    index_dict = dict()
    entry_dict = dict()

    for tag in set(entry_tags):
        index_dict[tag] = [i for i in range(len(entry_tags)) if entry_tags[i] == tag]
    
    entry_dict = parse_accessions(entry, index_dict, entry_dict)
    entry_dict = parse_names(entry, index_dict, entry_dict)
    entry_dict = parse_protein_gene(entry, index_dict, entry_dict)
    if 'dbReference' in index_dict.keys():
        entry_dict = parse_dbref(entry, index_dict, entry_dict)
    if 'keyword' in index_dict.keys():
        entry_dict = parse_keyword(entry, index_dict, entry_dict)

    entry_dict['entry_attributes'] = entry.attrib
    assert 'sequence' in index_dict.keys()
    entry_dict['sequence_attributes'] = entry[index_dict['sequence'][0]].attrib

    return entry_dict       

def parse_request_xml(request):
    proteome_dict = dict()
    with gzip.open(request.raw) as file:
        tree = ET.parse(file)
        root = tree.getroot()
        for entry in root:
            if entry.tag.replace(tag_pattern, '') == 'entry':
                new_dict = parse_entry(entry)
                accession = new_dict['accession']['primary_accession']
                proteome_dict[accession] = new_dict
            else:
                assert entry.tag.replace('{http://uniprot.org/uniprot}', '') == 'copyright'
                proteome_dict['copyright'] = entry.text

    return proteome_dict

## Add ko annotations using the KEGG species annotation
def fetch_ko_linktable(all_annotations, log_buffer):
    ko_api = 'https://rest.kegg.jp/link/ko/'
    kegg_entries = dict()
    def extract_kegg_(entry):
        if isinstance(entry, dict):
            if 'db_references' in entry.keys():
                if 'KEGG' in entry['db_references'].keys():
                    return list(entry['db_references']['KEGG'].keys())
                else:
                    return None
            else:
                return None
        else:
            None
        
    kegg_entries = [extract_kegg_(entry) for entry in all_annotations.values()]
    kegg_entries = [k for k in kegg_entries if k is not None]
    kegg_entries = [x for k in kegg_entries for x in k if k is not None]
    kegg_dbs = list(set([x.split(':')[0] for x in kegg_entries]))
    link_tables_dict, link_tables_dict_ = dict(), dict()
    if kegg_dbs == []:
        log_buffer[len(log_buffer)] = f'No KEGG annotations in the xml file.'
    for db in kegg_dbs:
        url = f'{ko_api}{db}'
        linkdf_path = Path(FLAGS.out_dir) / f'ko_link_tables/ko_{db}_link.txt'
        if not linkdf_path.exists():
            linkdf_path.parent.mkdir(parents=True, exist_ok=True)
            request = requests.get(url, stream = True)
            if request.status_code == 404:
                log_buffer[len(log_buffer)] = f'Link table not found (404) at url {url}. Please check manually.'
            elif request.status_code != 200:
                log_buffer[len(log_buffer)] = f'Encountered {request.status_code} error when fetching link table from {url}.'
            else:
                log_buffer[len(log_buffer)] = f'Found the ko and {db} link table. Writing to {linkdf_path.name}.'
                link_tables_dict[db] = (request.content, linkdf_path)
                link_tables_dict_[db] = request.content
        else:
            log_buffer[len(log_buffer)] = f'Already downloaded the ko to {db} link table.'
            with linkdf_path.open('rb') as file:
                link_tables_dict_[db] = file.read()
    ko_dict = dict()
    for db in link_tables_dict_.keys():
        links = [x.split(b'\t') for x in link_tables_dict_[db].split(b'\n')][:-1]
        ko_dict = ko_dict | {x[0].decode('utf-8') : x[1].decode('utf-8') for x in links}
    for accession, entry in all_annotations.items():
        if isinstance(entry, dict):
            if 'db_references' in entry.keys():
                if 'KEGG' in entry['db_references'].keys():
                    kegg_anns = [x for x in list(entry['db_references']['KEGG'].keys()) if x in ko_dict.keys()]
                    addition = {ko_dict[kegg_ann] : {"id" : ko_dict[kegg_ann], "mapped_from" : kegg_ann} for kegg_ann in kegg_anns}
                    if len(addition) > 0:
                        entry['db_references']['ko'] = addition
                        all_annotations[accession] = entry
    return log_buffer, link_tables_dict, all_annotations

################## Functions to parse fasta ##################

def download_proteome(urls, log_buffer):
    start_time = time.time() 
    fasta_content = []
    accessions = []
    request = requests.get(urls[0], stream = True)
    if request.status_code == 200:
        log_buffer[len(log_buffer)] = f'Downloading proteome from {urls[0]}.'
        with gzip.open(request.raw, 'rb') as file:
            read_content = file.read()
            fasta_content += [read_content]
            # xx = read_content.split(b'|')
            # accessions += [xx[i].decode('utf-8') for i in range(len(xx)) if i % 2 == 1]
            split_fasta = read_content.split(b'\n>')
            accession = split_fasta[0].split(b'|')[1].decode('utf-8')
            accessions += [accession] + [fasta_data.split(b'|')[1].decode('utf-8') for fasta_data in split_fasta[1:]]
            log_buffer[len(log_buffer)] = f'Took {time.time()-start_time} seconds to download proteome from {urls[0]}.'
    else:
        log_buffer[len(log_buffer)] = f'Proteome not found at {urls[0]}!!'
    start_time = time.time()
    request = requests.get(urls[1], stream = True)
    if request.status_code == 200:
        log_buffer[len(log_buffer)] = f'Downloading proteome from {urls[1]}.'
        with gzip.open(request.raw, 'rb') as file:
            read_content = file.read()
            fasta_content += [read_content]
            # xx = read_content.split(b'|')
            # accessions += [xx[i].decode('utf-8') for i in range(len(xx)) if i % 2 == 1]
            split_fasta = read_content.split(b'\n>')
            accession = split_fasta[0].split(b'|')[1].decode('utf-8')
            accessions += [accession] + [fasta_data.split(b'|')[1].decode('utf-8') for fasta_data in split_fasta[1:]]
            log_buffer[len(log_buffer)] = f'Took {time.time()-start_time} seconds to download proteome from {urls[1]}.'
    else:
        log_buffer[len(log_buffer)] = f'No additional additional fasta found.'
    return fasta_content, accessions, log_buffer

################## Functions to download fasta and xml in parallel ##################

def get_proteome_and_annotations(taxid, process_name, number_tried):
    log_buffer = dict()
    log_buffer[len(log_buffer)] = f'Process {process_name}. Working on taxa {taxid}. This is try number {number_tried+1}.'
    url_pattern = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/reference_proteomes'
    lineage = proteome_table['Taxonomic lineage'][taxid].split(', ')
    proteome_id = proteome_table['Proteome Id'][taxid]

    ## Viruses show first in lineage, all others show second.
    if lineage[0] == 'Viruses':
        lineage = 'Viruses'
    else:
        lineage = lineage[1]

    url1 = f'{url_pattern}/{lineage}/{proteome_id}/{proteome_id}_{taxid}.fasta.gz'
    url2 = f'{url_pattern}/{lineage}/{proteome_id}/{proteome_id}_{taxid}_additional.fasta.gz'
    urls = [url1, url2]
    proteome_path = Path(FLAGS.out_dir) / f'{proteome_id}_taxaid_{taxid}_proteome.fasta'
    log_buffer[len(log_buffer)] = f'Fetching proteome of {taxid}'
    if proteome_path.exists():
        log_buffer[len(log_buffer)] = f'The file {str(proteome_path)} already exists.'
        fasta_content = []
        accessions = proteome_path
    else:
        fasta_content, accessions, log_buffer = download_proteome(urls, log_buffer)

    annotations_path = Path(FLAGS.out_dir) / f'{proteome_id}_taxaid_{taxid}_annotations.json'
    all_annotations = dict()
    if annotations_path.exists():
        log_buffer[len(log_buffer)] = f'Process {process_name}. Annotation path {annotations_path} already exists'
        all_annotations = "Path already exists"
    else:
        url = f'{url_pattern}/{lineage}/{proteome_id}/{proteome_id}_{taxid}.xml.gz'
        log_buffer[len(log_buffer)] = f'Process {process_name}. Fetching annotations for taxid {taxid} from url {url}'

        request = requests.get(url, stream = True)
        if request.status_code == 404:
            log_buffer[len(log_buffer)] = f'Process {process_name}. Annotations not found at url {url}. Probably {taxid} is not a reference proteome any longer.'
        elif request.status_code != 200:
            log_buffer[len(log_buffer)] = f'Process {process_name}. Encountered {request.status_code} error when fetching annotations for {taxid} from url {url}'
        else:
            start_time = time.time()
            all_annotations = parse_request_xml(request)
            log_buffer[len(log_buffer)] = f'Process {process_name}. Took {time.time()-start_time} seconds to parse the main xml.'

        additional_annotations = dict()
        url = f'{url_pattern}/{lineage}/{proteome_id}/{proteome_id}_{taxid}_additional.xml.gz'
        log_buffer[len(log_buffer)] = f'Process {process_name}. Fetching additional annotations for taxid {taxid} from url {url}'
        request = requests.get(url, stream = True)
        if request.status_code == 404:
            log_buffer[len(log_buffer)] = f'Process {process_name}. No additional annotation for {taxid} found.'
        elif request.status_code != 200:
            log_buffer[len(log_buffer)] = f'Process {process_name}. Encountered error {request.status_code} when fetching additional annotations for {taxid} from url {url}\n {request.reason}'

        if request.status_code == 200:
            start_time = time.time()
            additional_annotations = parse_request_xml(request)
            log_buffer[len(log_buffer)] = f'Process {process_name}. Took {time.time()-start_time} seconds to parse the additional xml.'
        if additional_annotations:
            assert set(additional_annotations.keys()).intersection(set(all_annotations.keys())) == {'copyright'} 
            assert additional_annotations['copyright'] == all_annotations['copyright']
            for new_accession in additional_annotations.keys():
                all_annotations[new_accession] = additional_annotations[new_accession]

    if annotations_path.exists():
        with annotations_path.open('rb') as f:
            data = f.read()
            all_annotations = orjson.loads(data)

    log_buffer, link_tables_dict, all_annotations = fetch_ko_linktable(all_annotations, log_buffer)
    log_buffer = check_integrity(accessions, all_annotations, log_buffer)
    out = (taxid, all_annotations, annotations_path, fasta_content, proteome_path, link_tables_dict, log_buffer)
    return out

def check_integrity(accessions, all_annotations, log_buffer):
    if not isinstance(accessions, list):
        proteome_path = accessions
        with proteome_path.open('rb') as proteome_fasta:
            loaded_fasta = proteome_fasta.read()
            split_fasta = loaded_fasta.split(b'\n>')
            accession = split_fasta[0].split(b'|')[1].decode('utf-8')
            accessions = [accession] + [fasta_data.split(b'|')[1].decode('utf-8') for fasta_data in split_fasta[1:]]
    log_buffer[len(log_buffer)] = f'Found {len(accessions)} accessions in fasta.'
    annotation_accessions = set(all_annotations.keys()) - {'copyright'}
    log_buffer[len(log_buffer)] = f'Found {len(annotation_accessions)} accessions in json.'
    if annotation_accessions != set(accessions):
        log_buffer[len(log_buffer)] = f'Failed integrity check. Annotation accessions and fasta accessions are not the same.'
        log_buffer[len(log_buffer)] = f'These are the accessions which are not in the intersection: {annotation_accessions-set(accessions)} {set(accessions)-annotation_accessions}.'
    else:
        log_buffer[len(log_buffer)] = f'The same accessions are present in both.'
    return log_buffer

def download_(taxids, out_dict, process_name):
    this_out = dict()
    (taxids, failed_taxa_) = taxids
    for taxid in taxids:
        if taxid in failed_taxa_.keys():
            out = get_proteome_and_annotations(taxid, process_name, failed_taxa_[taxid])
        else:
            out = get_proteome_and_annotations(taxid, process_name, 0)
        this_out[taxid] = out
    out_dict[process_name] = this_out    


################## Functions to handle the splitting of tasks for multiprocess ##################


def next_name(N_processes, current_names):
    all_names = [str(i) for i in range(N_processes) if str(i) not in current_names]
    return current_names + [all_names[0]]

class taxa_cache():
    def __init__(self, taxids = [], cache_size = 10, N_retries = 10):
        self.taxids = taxids
        self.cache_size = cache_size
        self.total_size = len(taxids)
        self.all_caches = [taxids[i*cache_size:i*cache_size+cache_size] for i in range(int(len(taxids)/cache_size) + 1)]
        self.all_caches = [x for x in self.all_caches if x != []]
        self.current_cache = 0
        self.finished = False
        self.maybe_finished = False
        self.failed_taxa_ = {}
        self.super_failed = set({})
        self.N_retries = N_retries
        self.past_caches = dict()

    def next_cache(self, processname):
        if self.current_cache != len(self.all_caches):
            out = (self.all_caches[self.current_cache], self.failed_taxa_)
            self.current_cache += 1
            self.past_caches[processname] = out[0]
        else:
            out = None
            if all([x > self.N_retries for x in self.failed_taxa_.values()]) and self.maybe_finished:
                self.finished = True
            self.maybe_finished = True
        return out
    
    def failed_taxa(self, taxid):
        if taxid in self.failed_taxa_.keys():
            self.failed_taxa_[taxid] += 1
            if self.failed_taxa_[taxid] > self.N_retries:
                self.super_failed = self.super_failed | {taxid}
        else:
            self.failed_taxa_[taxid] = 1
        if taxid not in self.super_failed:
            last_cache = self.all_caches[-1]
            if len(last_cache) < self.cache_size and self.current_cache < len(self.all_caches):
                last_cache += [taxid]
                self.all_caches[-1] = last_cache
            else:
                self.all_caches += [[taxid]]


################## Main loop ##################

all_process = []
process_names = []
taxa_list = list(proteome_table.index)
shuffle(taxa_list)
# taxa_list = [3702]
# taxa_list = [2695836]
# taxa_list = [1892558, 2053570]
if __name__ == '__main__':
    still_working = True
    manager = multiprocessing.Manager()
    out_dict = manager.dict()
    taxa_cache_ = taxa_cache(taxa_list, cache_size=int(FLAGS.cache_size))
    N_parsed = 0

    while still_working:
        if len(all_process) < N_processes:
            process_name = next_name(N_processes, process_names)[-1]
            database_taxids = taxa_cache_.next_cache(process_name)
            process_args = (database_taxids, out_dict, process_name)
            if database_taxids is not None:
                new_process = multiprocessing.Process(target = download_, args = process_args)
                all_process = all_process + [new_process]
                process_names = process_names + [process_name]
                new_process.start()
            else:
                # print('Transitioning to failed taxids once other processes finish. N_process set to 1.')
                N_processes = 1
        for process in all_process:
            process.join(timeout = 0.1)
        for process, process_name in zip(all_process, process_names):
            if process.exitcode is not None:
                all_process = [x for x in all_process if x is not process]
                process_names = [x for x in process_names if x is not process_name]
                this_out = out_dict[process_name]
                if process.exitcode != 0:
                    taxids = taxa_cache_.past_caches[process_name]
                    for taxid in taxids:
                        taxa_cache_.failed_taxa(taxid)
                        with log_path.open('a') as output:
                            output.write(f'A subprocess failed. The taxids it was working on ({taxids}) will be retried.\n')
                            output.write('======================================================\n')
                else:
                    for out in this_out.values():
                        taxa_failed = False
                        (taxid, all_annotations, annotations_path, fasta_content, proteome_path, link_tables_dict, log_buffer) = out
                        N_parsed = N_parsed + 1
                        ## Write annotations
                        if annotations_path.exists():
                            print(f'Annotations from {taxid} already exist. This is proteome number {N_parsed}.')
                        elif len(all_annotations) != 0:
                            print(f'Saving annotations from {taxid} into {str(annotations_path)}. This is proteome number {N_parsed}.')
                            with open(annotations_path, "wb") as outfile: 
                                outfile.write(orjson.dumps(all_annotations)) 
                        else:
                            print(f'Annotations from {taxid} not found. This is proteome number {N_parsed}.')
                            taxa_failed = True
                        ## Write proteome
                        if proteome_path.exists():
                            print(f'Proteome from {taxid} already exists. This is proteome number {N_parsed}.')
                        elif fasta_content != []:
                            print(f'Saving proteome from {taxid} into {str(proteome_path)}. This is proteome number {N_parsed}.')
                            with open(proteome_path, 'wb') as outfile:
                                for read_content in fasta_content:
                                    outfile.write(read_content)
                        else:
                            print(f'Proteome from {taxid} not found. This is proteome number {N_parsed}.')
                            taxa_failed = True
                        ## Write KO link tables
                        if taxa_failed:
                            taxa_cache_.failed_taxa(taxid)
                        for x in link_tables_dict.values():
                            (content, linkdf_path) = x
                            with linkdf_path.open('wb') as file:
                                file.write(content)       
                        ## Write log
                        with log_path.open('a') as output:
                            for new_line in log_buffer.values():
                                output.write(f'{new_line}\n')
                            output.write('======================================================\n')
        if len(all_process) == 0 and taxa_cache_.finished:
                still_working = False
                failed_taxa = list(taxa_cache_.super_failed)
                if failed_taxa == []:
                    print(f'All taxa downloaded!')
                    with log_path.open('a') as output:
                            output.write(f'\n')
                            output.write('======================================================\n')
                            output.write(f'All taxa downloaded!\n')
                            output.write(f'The following taxids failed to downloaded initially, but all were successfully downloaded in under {taxa_cache_.N_retries} retries.\n')
                            output.write(f'{taxa_cache_.failed_taxa_}\n')
                            output.write('======================================================\n')
                else:
                    print(f'These taxa failed to download after {taxa_cache_.N_retries} retries! {failed_taxa}')
                    with log_path.open('a') as output:
                            output.write(f'\n')
                            output.write('======================================================\n')
                            output.write(f'The following taxids failed to download after {taxa_cache_.N_retries+1} attempts! {failed_taxa}\n')
                            output.write('======================================================\n')

