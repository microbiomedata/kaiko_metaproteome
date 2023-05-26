import gzip
import time
import pandas as pd
import numpy as np
import os


# @profile
def aggregate_fasta(ref_fasta, diamon_tally, fout, coverage_target, top_strains, ncbi_taxa_folder, taxa_key, kingdom_list = []):
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
    taxids = []

    for fasta_addition in coverage_steps:
        fasta_addition = tdf[tdf['running_coverage'] == fasta_addition]
        primary_species = ['Primary' in str for str in fasta_addition['notes']]
        secondary_strains = ['Strain' in str for str in fasta_addition['notes']]
        if any(secondary_strains):
            fasta_addition = fasta_addition[secondary_strains].nlargest(top_strains, 'hits')
        else:
            fasta_addition = fasta_addition[primary_species].nlargest(top_strains, 'hits')
        taxids = taxids + [int(fasta_addition['taxid'].iloc[0])]

    fout_proteome = fout.parent / (fout.stem + '_proteome.txt')
    print("NCBI-TaxIDs:", taxids)
    proteome = get_taxa_proteome(ncbi_taxa_folder / 'uniref100_member_taxa_tbl.csv', taxids, fout_proteome)

    key_for_taxid = '{}='.format(taxa_key)
    print("Key for parsing Tax IDs:", key_for_taxid)
    # ofile = open(fout, 'w')
    ofile = fout.open('w')

    num_seqs = 0
    start_time = time.time()
    proteome = pd.read_csv(fout_proteome, sep = "\t", header = None)
    proteome = set(proteome[0].values)
    
    with gzip.open(ref_fasta, 'rb') as file:
        try:
            is_selected = False
            for bline in file:
            
                # if line.startswith('>UniRef'):
                if bline[0] == 62:  # to find `>`
                    num_seqs += 1
                    line = bline.decode("utf-8")
                    uid = line.split('>')[1].split(' ')[0]
                    taxid = line.split(key_for_taxid)[1].split(' ')[0]
                    if uid in ["", "N/A"]:
                        is_selected = False
                    elif uid in proteome:
                        is_selected = True
                        ofile.write(line)
                    else:
                        is_selected = False
                else:
                    if is_selected:
                        line = bline.decode("utf-8")
                        ofile.write(line)
                if (num_seqs % 1000000) == 0:
                    print("{}M sequences has been parsed. {:.1f}min".format(num_seqs//1e6, (time.time()-start_time)/60))
        except Exception as e:
            print(line, e)

    ofile.close()


def get_taxa_proteome(member_tbl_file, unique_taxids, fout_proteome):
    chunksize = 1000000
    unique_taxids = set([str(x) for x in unique_taxids])

    print('Collecting protein names in ' + fout_proteome.name + '\n')
    nchunks = 1
    output = fout_proteome.open('w')
    for chunk in pd.read_csv(member_tbl_file, chunksize=chunksize):
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

# def write_proteins(database_file, taxa, fout):
#     index_s_path = Path('Kaiko_volume/Kaiko_stationary_files/uniref100_index_s.txt')
#     index_path = Path('Kaiko_volume/Kaiko_stationary_files/uniref100_index.txt')
#     index_s = pd.read_csv(index_s_path, sep = '\t', header = None)
#     index_s.index = index_s[0]
#     index_pos = int(index_s.loc[f'taxid_{taxa}_positions'][1])
#     index = index_path.open('r')
#     index.seek(index_pos)
#     protein_pos = index.readline()
#     protein_pos = protein_pos.split(';;')[1].split(';')[:-1]

#     with fout.open('a') as output_fasta:
#         for pos in protein_pos:
#             database_file.seek(int(pos))
    
#             line = database_file.readline()
#             line = line.decode("utf-8")

#             reading_protein = True
#             while reading_protein:
#                 output_fasta.write(line)
#                 line = database_file.readline()
#                 if line[0] == 62:  # to find `>`
#                     reading_protein = False
#                 else:
#                     line = line.decode("utf-8")

# def write_single_protein(database_file, pos, output_fasta):
#     database_file.seek(pos)
    
#     line = database_file.readline()
#     line = line.decode("utf-8")

#     reading_protein = True
#     while reading_protein:
#         output_fasta.write(line)
#         line = database_file.readline()
#         if line[0] == 62:  # to find `>`
#             reading_protein = False
#         else:
#             line = line.decode("utf-8")

# from pathlib import Path, PureWindowsPath

# prefix = "S1_NM0001_NMDC_MixedCommunities_R3_mgf"
# diamond_search_out = Path("Kaiko_volume/Kaiko_intermediate/" + prefix + "_diamond_search_output.dmd")
# kaiko_tally = Path("Kaiko_volume/Kaiko_intermediate/" + prefix + "_kaiko_prediction_top_taxa.csv")
# ncbi_taxa_folder = Path(PureWindowsPath("Kaiko_volume/Kaiko_stationary_files/ncbi_taxa").as_posix())
# nprot = '{:.5e}'.format(int(150000))
# kaiko_tally = Path("Kaiko_volume/Kaiko_intermediate/" + prefix + "_kaiko_prediction" + f'_top_taxa_nprot_{nprot}_top_{5}_strains.csv')
# ref_fasta = Path(PureWindowsPath('Kaiko_volume/Kaiko_stationary_files/uniref100.fasta.gz').as_posix())
# kaiko_final_output = Path("Kaiko_volume/Kaiko_output/" + prefix + "_kaiko_output.fasta")

# aggregate_fasta(ref_fasta,
#                 kaiko_tally,
#                 kaiko_final_output,
#                 0.66,
#                 'TaxID',
#                 [])