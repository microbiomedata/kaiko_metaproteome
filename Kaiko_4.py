import gzip
import time
import pandas as pd
import numpy as np
import os
import indexed_gzip as igzip

# @profile
def aggregate_fasta(ref_fasta, diamon_tally, output_fasta_path, coverage_target, top_strains, ref_fasta_igzip_index, index_path, index_s_path, kingdom_list = []):
    """
    python -u ExtractProteinSequenceByTaxID.py --ref_fasta uniref100.fasta.gz --taxafile Weintraub_kaiko_25p_UniRef100_toptaxa.csv --fout Weintraub_kaiko_25p_UniRef100_Bacteria_top100.fasta --ntops 100 -l Bacteria --key TaxID
    """
    df = pd.read_csv(diamon_tally)
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
    
    with igzip.IndexedGzipFile(str(ref_fasta), index_file = str(ref_fasta_igzip_index)) as database_file:
        write_taxa(database_file, taxids, output_fasta_path, index_s_path, index_path)    


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