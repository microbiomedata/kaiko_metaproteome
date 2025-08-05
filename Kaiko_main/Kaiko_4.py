import gzip
import time
import os
import orjson

# import gffpandas.gffpandas
import pandas as pd
import numpy as np
import indexed_gzip as igzip

from pathlib import PureWindowsPath, Path

byte2gigabyte = 1073741824
EXCLUDED_RANKS = ['family','order','class','phylum','kingdom','superkingdom']
# @profile
def aggregate_fasta(ref_fasta, kaiko_tally, member_csv, output_fasta_path, output_gff_path, sheet_name, target_coverage,
                    kingdom_list, DB, mode):
    df = pd.read_excel(kaiko_tally, sheet_name=sheet_name)
    rank_conditions = ~df['rank'].isin(EXCLUDED_RANKS)

    if len(kingdom_list) > 0:
        conditions = df.superkingdom.isin(kingdom_list)
        conditions |= df.kingdom.isin(kingdom_list)
        tdf = df[rank_conditions & conditions]
    else:
        tdf = df[rank_conditions]

    if type(target_coverage) is str:
        target_size = round(float(target_coverage.replace('Gb', '')), 2)
    else:
        target_coverage = round(float(target_coverage), 1)
        tdf = tdf[tdf['running_coverage'] < target_coverage/100]
        target_size = None

    print(tdf.head(20))
    coverage_steps = [*set(tdf['running_coverage'].values)]
    coverage_steps.sort()
    taxids = []

    for fasta_addition in coverage_steps:
        fasta_addition = tdf[tdf['running_coverage'] == fasta_addition]
        primary_species = ['Primary' in str for str in fasta_addition['notes']]
        secondary_strains = ['Secondary' in str for str in fasta_addition['notes']]
        if any(secondary_strains):
            fasta_addition = fasta_addition[secondary_strains].nlargest(5, 'hits')
        else:
            fasta_addition = fasta_addition[primary_species].nlargest(1, 'hits')
        taxids = taxids + [int(selected_taxid) for selected_taxid in fasta_addition['taxid'].values]

    print(f'collecting fasta from {len(taxids)} taxa\n')

    if DB == "Uniref100":
        uniref_fasta(taxids, output_fasta_path, ref_fasta, member_csv, mode)
    else:
        write_log = output_fasta_path.parent / 'test_log.txt'
        with output_fasta_path.open('wb') as output_fasta:
            with output_gff_path.open('w') as output_gff:
                with write_log.open('w') as log:
                    for taxid in taxids:
                        log_dict = write_proteome_and_annotations(taxid, ref_fasta, output_fasta, output_gff, dict())
                        for line in log_dict.values():
                            log.write(f'{line}\n')


def write_proteome_and_annotations(taxid, ref_fasta, output_fasta, output_gff, log_dict):
    fasta_paths = list(ref_fasta.glob(f'*_taxaid_{taxid}_proteome.fasta'))
    if len(fasta_paths) == 0:
        log_dict[len(log_dict)] = f'No fasta found for taxid {taxid}.'
        fasta_accessions = set({})
    if len(fasta_paths) >= 1:
        if len(fasta_paths) > 1:
            log_dict[len(log_dict)] = f'Found more than one fasta for taxid {taxid}.'
        loaded_fasta = fasta_paths[0].open('rb').read()
        split_fasta = loaded_fasta.split(b'\n>')
        fasta_dict = dict()

        accession = split_fasta[0].split(b'|')[1].decode('utf-8')
        fasta_dict[accession] = split_fasta[0] + b'\n'
        for fasta_data in split_fasta[1:]:
            accession = fasta_data.split(b'|')[1].decode('utf-8')
            fasta_dict[accession] = b'>' + fasta_data + b'\n'
        fasta_accessions = set(fasta_dict.keys())

    json_paths = list(ref_fasta.glob(f'*_taxaid_{taxid}_annotations.json'))
    if len(json_paths) == 0:
        log_dict[len(log_dict)] = f'No json found for taxid {taxid}.'
        annotation_accessions = set({})
    if len(json_paths) >= 1:
        if len(fasta_paths) > 1:
            log_dict[len(log_dict)] = f'Found more than one json for taxid {taxid}.'
        with json_paths[0].open('rb') as f:
            data = f.read()
            all_annotations = orjson.loads(data)
            annotation_accessions = set(all_annotations.keys()) - {'copyright'}

    if annotation_accessions != fasta_accessions:
        log_dict[len(log_dict)] = f'Warning, mismatch between annotations and fasta'
    if annotation_accessions - fasta_accessions:
        log_dict[len(log_dict)] = f'The set {annotation_accessions - fasta_accessions} is in the annotations and not the fasta.'
    if fasta_accessions - annotation_accessions:
        log_dict[len(log_dict)] = f'The set {fasta_accessions - annotation_accessions} is in the fasta and not the annotations.'
    
    common_accessions = annotation_accessions.intersection(fasta_accessions)
    for accession in common_accessions:
        fasta_data = fasta_dict[accession]
        output_fasta.write(fasta_data)
        gff_line = prepare_gff_line(all_annotations[accession])
        output_gff.write(gff_line)
    log_dict[len(log_dict)] = f'Finished writing data for {taxid}.'
    return log_dict


def prepare_gff_line(accession_annotations):
    accession = accession_annotations['accession']['primary_accession']
    protein_name = accession_annotations['recommended_name']['full_name'].replace(';', ' ')
    gff_line = f'{accession}\tKaiko2_Reference_Proteomes\tCDS\t.\t.\t.\t.\t.'
    gff_line = f'{gff_line}\tID={accession};product={protein_name};product_source=UniProt'
    references = accession_annotations['db_references']

    def add_to_line(gff_line, category, name):
        if category == 'cog' and 'eggNOG' in references.keys():
            annotations = list(references['eggNOG'].keys())
            annotations = [x for x in annotations if 'COG' in x]
            if len(annotations) > 0:
                gff_line = f'{gff_line};{name}={annotations[0]}'
                for ann in annotations[1:]:
                    gff_line = f'{gff_line},{ann}'
            else:
                gff_line = f'{gff_line};{name}=NA'
        elif category in references.keys():
            annotations = list(references[category].keys())
            if category == 'ko':
                    annotations = [ann.replace('ko', 'KO') for ann in annotations]
            gff_line = f'{gff_line};{name}={annotations[0]}'
            for ann in annotations[1:]:
                gff_line = f'{gff_line},{ann}'
        else:
            gff_line = f'{gff_line};{name}=NA'
        return gff_line

    gff_line = add_to_line(gff_line, 'Pfam', 'pfam')
    gff_line = add_to_line(gff_line, 'ko', 'ko')
    gff_line = add_to_line(gff_line, 'EC', 'ec_number')
    gff_line = add_to_line(gff_line, 'cog', 'cog')
    gff_line = add_to_line(gff_line, 'KEGG', 'kegg')
    gff_line = add_to_line(gff_line, 'eggNOG', 'eggnog')
    gff_line = gff_line.replace('\n', '')
    gff_line = f'{gff_line}\n'
    return gff_line


def uniref_fasta(taxids, output_fasta_path, ref_fasta, member_csv, mode):
    ofile = open(output_fasta_path, 'w')
    key_for_taxid = 'TaxID'
    num_seqs = 0
    start_time = time.time()
    if mode == "member":
        uniref_ids = get_uniref_ids(member_csv, [str(x) for x in taxids])
    
    with gzip.open(ref_fasta, 'rb') as file:
        try:
            is_selected = False
            for bline in file:
            
                # if line.startswith('>UniRef'):
                if bline[0] == 62:  # to find `>`
                    num_seqs += 1
                    line = bline.decode("utf-8")
                    if mode == "member":
                        uniref_id = line.split(' ')[0].split('>')[1]
                        if uniref_id in ["", "N/A"]:
                            is_selected = False
                        elif uniref_id in uniref_ids.keys():
                            is_selected = True
                            ofile.write(line)
                        else:
                            is_selected = False
                    elif mode == "common":
                        taxid = line.split(key_for_taxid)[1].split(' ')[0]
                        if taxid in ["", "N/A"]:
                            is_selected = False
                        elif int(taxid) in taxids:
                            is_selected = True
                            ofile.write(line)
                        else:
                            is_selected = False
                else:
                    if is_selected:
                        line = bline.decode("utf-8")
                        ofile.write(line)
                if (num_seqs % 10000000) == 0:
                    print("{}M sequences has been parsed. {:.1f}min".format(num_seqs//1e6, (time.time()-start_time)/60))
        except Exception as e:
            print(line, e)

    ofile.close()


def get_uniref_ids(member_tbl_file, taxids):
    chunksize = 1000000

    stime = time.time()
    uniref_ids = dict()
    num_iters = 0
    for chunk in pd.read_csv(member_tbl_file, chunksize=chunksize):
        chunk_ids = chunk.uid.tolist()
        for ii, members in enumerate(chunk.members.tolist()):
            members = members.split(':')
            taxids_ = [x in taxids for x in members]
            if any(taxids_):
                uniref_ids[chunk_ids[ii]] = taxids_
        num_iters += 1
    with open(Path("protein_to_taxaid"), "wb") as outfile: 
        outfile.write(orjson.dumps(uniref_ids)) 
    return uniref_ids

