import gzip
import time
import os
import json

import pandas as pd
import numpy as np
import indexed_gzip as igzip

from pathlib import PureWindowsPath, Path

def extract_single_feature_annotations(annotation_dict, keywords):
    allowed_keywords = ['GO', 'KEGG', 'Reactome']
    # allowed = [x in allowed_keywords for x in keywords]
    # assert(all(allowed))
    ## Required, always in output
    simplified_dict = dict()
    simplified_dict['feature_id'] = annotation_dict['Identification (ID)']['entry_ID']
    simplified_dict['primary_accession'] = annotation_dict['Accession (AC)']['primary_accession']
    simplified_dict['sequence_length'] = annotation_dict['Identification (ID)']['seq_len']
    simplified_dict['version'] = annotation_dict['Date (DT)']['version_number']
    simplified_dict['version_date'] = annotation_dict['Date (DT)']['version_date']
    simplified_dict['review_status'] = annotation_dict['Identification (ID)']['review_status']

    if 'Gene name (GN)' in annotation_dict.keys():
        simplified_dict['gene_name'] = annotation_dict['Gene name (GN)']['official_uniprot_gene']
    else:
        simplified_dict['gene_name'] = 'no_genename_in_annotations'

    ## RecName exists for reviewed entries (recommended name, SwissProt)
    if 'RecName_full' in annotation_dict['Description (DE)'].keys():
        simplified_dict['submission_name'] = annotation_dict['Description (DE)']['RecName_full']
    ## SubName exists for unreviewed entries (name provided by submitter, TrEMBL)
    elif 'SubName_full' in annotation_dict['Description (DE)'].keys():
        simplified_dict['submission_name'] = annotation_dict['Description (DE)']['SubName_full']

    for keyword in keywords:
        if keyword == 'GO' and 'Gene Ontology annotations (DR)' in annotation_dict.keys():
            simplified_dict['GO_annotations'] = []
            go_annotations = annotation_dict['Gene Ontology annotations (DR)']
            for go_term in go_annotations.keys():
                go_code = go_annotations[go_term]['GO_code']
                go_annotation = f'{go_term} -- go_code={go_code}'
                simplified_dict['GO_annotations'] = simplified_dict['GO_annotations'] + [go_annotation]
        if keyword == 'Reactome' and 'Reactome annotations (DR)' in annotation_dict.keys():
            simplified_dict['Reactome_annotations'] = []
            reactome_annotations = annotation_dict['Reactome annotations (DR)']
            for rt_term in reactome_annotations.keys():
                reactome_code = reactome_annotations[rt_term]['Reactome_code']
                reactome_annotation = f'{rt_term} -- reactome_code={reactome_code}'
                simplified_dict['Reactome_annotations'] = simplified_dict['Reactome_annotations'] + [reactome_annotation]
        if keyword == 'KEGG' and 'KEGG identifier (DR)' in annotation_dict.keys():
            simplified_dict['KEGG_annotations'] = [annotation_dict['KEGG identifier (DR)']]
    return simplified_dict


def create_annotation_file(ref_proteome_log, ref_fasta, taxids, output_annotation_path, requested_annotations):
    proteome_df = pd.read_excel(ref_proteome_log)
    proteome_df.index = proteome_df['Organism Id']
    total_proteins = 0
    selected_annotations = dict()

    for taxid in taxids:
        proteome_id = proteome_df.loc[taxid]['Proteome Id']
        annotation_path = ref_fasta / f'{proteome_id}_taxaid_{taxid}_annotations.JSON'
        if annotation_path.exists():
            proteome_annotation = json.load(annotation_path.open())
            total_proteins = total_proteins + len(proteome_annotation.keys())
            for feature in proteome_annotation.keys():
                annotation_dict = proteome_annotation[feature]
                selected_annotations[feature] = extract_single_feature_annotations(annotation_dict, requested_annotations)
                selected_annotations[feature]['taxa_id'] = taxid
    with output_annotation_path.open('w') as output_annotations:
        json.dump(selected_annotations, output_annotations)
    output_annotation_gff_path = output_annotation_path
    output_annotation_gff_path = output_annotation_gff_path.parent / (output_annotation_gff_path.stem + '.gff')

    with output_annotation_gff_path.open('w') as output_annotations:
        feature_source = 'kaiko_1.2'
        feature_type = '.'
        feature_start = '.'
        feature_end = '.'
        feature_score = '.'
        feature_strand = '.'
        feature_phase = '.'
        for feature in selected_annotations.keys():
            if selected_annotations[feature]['review_status'] == 'Reviewed':
                product_source = 'UniProt Swiss-Prot'
            else:
                product_source = 'UniProt TrEMBL'
            product = selected_annotations[feature]['submission_name']
            feature_id = selected_annotations[feature]['feature_id']
            feature_gene = selected_annotations[feature]['gene_name']
            feature_ver = selected_annotations[feature]['version']
            feature_date = selected_annotations[feature]['version_date']
            taxid = selected_annotations[feature]['taxa_id']
            feature_attrributes = f'ID={feature};product={product};product_source={product_source}'
            feature_attrributes = f'{feature_attrributes};uniprot_id={feature_id};uniprot_gene={feature_gene}'
            feature_attrributes = f'{feature_attrributes};version={feature_ver};version_date={feature_date}'
            feature_attrributes = f'{feature_attrributes};taxa_id={taxid}'
            if 'GO_annotations' in selected_annotations[feature].keys():
                go_annotaions = selected_annotations[feature]['GO_annotations']
                feature_attrributes = f'{feature_attrributes};go_annotations={go_annotaions}'
            if 'KEGG_annotations' in selected_annotations[feature].keys():
                kegg_annotations = selected_annotations[feature]['KEGG_annotations']
                feature_attrributes = f'{feature_attrributes};kegg_annotations={kegg_annotations}'
            if 'Reactome_annotations' in selected_annotations[feature].keys():
                reactome_annotations = selected_annotations[feature]['Reactome_annotations']
                feature_attrributes = f'{feature_attrributes};reactome_annotations={reactome_annotations}'
            
            line = f'{feature}\t{feature_source}\t{feature_type}\t{feature_start}\t{feature_end}\t{feature_score}\t{feature_strand}\t{feature_phase}\t{feature_attrributes}\n'
            output_annotations.write(line)

# @profile
def aggregate_fasta(ref_fasta, ref_proteome_log, kaiko_tally, output_fasta_path, output_annotation_path, 
                    coverage_target, top_strains, ref_fasta_igzip_index, index_path, index_s_path, kingdom_list, mode):
    df = pd.read_excel(kaiko_tally, sheet_name='pident at least 100 percent')
    # df, _ = rank_to_lineage(df)

    rank_conditions = ~df['rank'].isin(EXCLUDED_RANKS)

    if len(kingdom_list) > 0:
        conditions = df.superkingdom.isin(kingdom_list)
        conditions |= df.kingdom.isin(kingdom_list)

        # tdf = df[rank_conditions & conditions].nlargest(coverage_target, 'hits', keep="all")
        tdf = df[rank_conditions & conditions]
    else:
        # tdf = df[rank_conditions].nlargest(coverage_target, 'hits', keep="all")
        tdf = df[rank_conditions]
    
    coverage_target = min(tdf[tdf['running_coverage'] > coverage_target]['running_coverage'])
    tdf = tdf[tdf['running_coverage'] <= coverage_target]

    print(tdf.head(20))
    coverage_steps = [*set(tdf['running_coverage'].values)]
    coverage_steps.sort()
    taxids = []

    for fasta_addition in coverage_steps:
        fasta_addition = tdf[tdf['running_coverage'] == fasta_addition]
        primary_species = ['Primary' in str for str in fasta_addition['notes']]
        secondary_strains = ['Secondary' in str for str in fasta_addition['notes']]
        if any(secondary_strains):
            fasta_addition = fasta_addition[secondary_strains].nlargest(top_strains, 'hits')
        else:
            fasta_addition = fasta_addition[primary_species].nlargest(top_strains, 'hits')
        taxids = taxids + [int(selected_taxid) for selected_taxid in fasta_addition['taxid'].values]

    print(f'collecting fasta from {len(taxids)} taxa\n')

    if mode == 'uniref100':
        with igzip.IndexedGzipFile(str(ref_fasta), index_file = str(ref_fasta_igzip_index)) as database_file:
            write_taxa(database_file, taxids, output_fasta_path, index_s_path, index_path)
    elif mode == 'ref_prot':
        write_taxa_ref_proteome(ref_proteome_log, ref_fasta, taxids, output_fasta_path)
        requested_annotations = ['GO', 'Reactome', 'KEGG']
        create_annotation_file(ref_proteome_log, ref_fasta, taxids, output_annotation_path, requested_annotations)


def get_taxa_proteome(member_tbl_file, unique_taxids, fout_proteome):
    chunksize = 1000000
    unique_taxids = set([str(x) for x in unique_taxids])

    print('Collecting protein names in ' + fout_proteome.name + '\n')
    nchunks = 0
    output = fout_proteome.open('w')
    for chunk in pd.read_csv(member_tbl_file, chunksize=chunksize):
        print(nchunks)
        proteome_index = [any([x in set(member_taxa.split(':')) for x in unique_taxids]) for member_taxa in chunk['members']]
        proteome = chunk[proteome_index]['uid'].values.tolist()
        for protein in proteome:
            output.write(f'{protein}\n')
        
        nchunks = nchunks + 1
    return proteome


def rank_to_lineage(df):
    coverage = {}
    for c in ['species','genus','family','order','class','phylum','kingdom','superkingdom']:
        idx = df['rank']==c
        df.loc[idx, c] = df.loc[idx, 'tax_name']
    for c in ['species','genus','family','order','class','phylum','superkingdom']:
        coverage[c] = 100*df[c].dropna().shape[0]/df.shape[0]
        print("{}: {:.2f}%".format(c, coverage[c]))
    return df, coverage

EXCLUDED_RANKS = ['family','order','class','phylum','kingdom','superkingdom']

def write_taxa(database_file, taxa_list, output_fasta_path, index_s_path, index_path, max_ram_lines = 8000):
    index_s = pd.read_csv(index_s_path, sep = '\t', header = None)
    index_s.index = index_s[0]
    index = index_path.open('r')
    all_pos = []

    print("Gathering positions")
    start_time = time.perf_counter()
    for taxa in taxa_list:
        index_pos = int(index_s.loc[f'taxid_{taxa}_positions'][1])
        index.seek(index_pos)
        protein_pos = index.readline()
        protein_pos = protein_pos.split(';;')[1].split(';')[:-1]
        all_pos = all_pos + protein_pos
    # Unique positions
    all_pos = list(set(all_pos))
    all_pos = [int(x) for x in all_pos]
    all_pos.sort()
    print(all_pos[1:100])
    print("Finished collecting locations")
    print(time.perf_counter() - start_time)

    with output_fasta_path.open('a') as output_fasta:
        lines = []
        count = 0
        start_time = time.perf_counter()
        for pos in all_pos:
            lines = lines + get_single_protein(database_file, pos)

            if (len(lines) > max_ram_lines):
                print("Writing to FASTA" + str(count))
                count = count + 1
                for line in lines:
                    output_fasta.write(line)
                lines = []
        print("Writing to FASTA (last time)")
        for line in lines:
            output_fasta.write(line)
        print(time.perf_counter() - start_time)

def write_taxa_ref_proteome(ref_proteome_log, ref_fasta, taxids, output_fasta_path, max_ram_lines = 8000):
    proteome_df = pd.read_excel(ref_proteome_log)
    proteome_df.index = proteome_df['Organism Id']
    total_proteins = 0
    with output_fasta_path.open('w') as output_fasta:
        for taxid in taxids:
            proteome_id = proteome_df.loc[taxid]['Proteome Id']
            fasta_path = ref_fasta / f'{proteome_id}_taxaid_{taxid}.fasta'
            with fasta_path.open() as taxa_proteome:
                for line in taxa_proteome:
                    if '>' in line:
                        total_proteins = total_proteins + 1
                    output_fasta.write(line)
    print(f'Finished writing total of {total_proteins} proteins')
        
        

def get_single_protein(database_file, pos):
    database_file.seek(pos)
    
    line = database_file.readline()
    line = line.decode("utf-8")
    lines = []

    reading_protein = True
    while reading_protein:
        lines = lines + [line]
        line = database_file.readline()
        if line[0] == 62:  # to find `>`
            reading_protein = False
        else:
            line = line.decode("utf-8")
    return(lines)

# from pathlib import Path, PureWindowsPath
# import indexed_gzip as igzip
# import time


# ref_fasta = Path(PureWindowsPath('Kaiko_volume/Kaiko_stationary_files/uniref100.fasta.gz'))
# ref_fasta_igzip_index = Path(PureWindowsPath('Kaiko_volume/Kaiko_stationary_files/uniref100_fasta_gzindex.gzidx'))
# # output_fasta_path = Path(PureWindowsPath('Kaiko_volume/ecoli_test_sept.FASTA'))
# index_s_path = Path('Kaiko_volume/Kaiko_stationary_files/uniref100_index_s.txt')
# index_path = Path('Kaiko_volume/Kaiko_stationary_files/uniref100_index.txt')

# # with igzip.IndexedGzipFile(str(ref_fasta), index_file = str(ref_fasta_igzip_index)) as file:
# #     write_taxa(file, ['562'], output_fasta_path, index_s_path, index_path)


# # prefix = "S1_NM0001_NMDC_MixedCommunities_R3_mgf"
# prefix = "Kansas_soil_no_gly"
# diamond_search_out = Path("Kaiko_volume/Kaiko_intermediate/" + prefix + "_diamond_search_output.dmd")
# kaiko_tally = Path("Kaiko_volume/Kaiko_intermediate/" + prefix + "_kaiko_prediction_top_taxa.csv")
# ncbi_taxa_folder = Path(PureWindowsPath("Kaiko_volume/Kaiko_stationary_files/ncbi_taxa").as_posix())
# nprot = '{:.5e}'.format(int(151000))
# kaiko_tally = Path("Kaiko_volume/Kaiko_intermediate/" + prefix + "_kaiko_prediction" + f'_top_taxa_nprot_{nprot}_top_{1}_strains.csv')
# ref_fasta = Path(PureWindowsPath('Kaiko_volume/Kaiko_stationary_files/uniref100.fasta.gz').as_posix())
# kaiko_final_output = Path("Kaiko_volume/Kaiko_output/" + prefix + "_kaiko_output.fasta")
# taxa_member_path = Path(PureWindowsPath('Kaiko_volume/Kaiko_stationary_files/ncbi_taxa/uniref100_member_taxa_tbl.csv'))

# aggregate_fasta(ref_fasta,
#                 kaiko_tally,
#                 kaiko_final_output,
#                 0.66,
#                 1,
#                 ncbi_taxa_folder,
#                 'TaxID',
#                 [])

# get_taxa_proteome(taxa_member_path, [562], Path(PureWindowsPath('Kaiko_volume/Kaiko_stationary_files/ecoli_sanity_check.txt')))

# aggregate_fasta(ref_fasta, diamon_tally, fout, coverage_target, top_strains, ncbi_taxa_folder, taxa_key, kingdom_list = [])