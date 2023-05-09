import gzip
import time
import pandas as pd
import numpy as np
import os


# @profile
def aggregate_fasta(ref_fasta, diamon_tally, fout, ntops, taxa_key, kingdom_list = []):
    """
    python -u ExtractProteinSequenceByTaxID.py --ref_fasta uniref100.fasta.gz --taxafile Weintraub_kaiko_25p_UniRef100_toptaxa.csv --fout Weintraub_kaiko_25p_UniRef100_Bacteria_top100.fasta --ntops 100 -l Bacteria --key TaxID
    """
    df = pd.read_csv(diamon_tally)
    df, _ = rank_to_lineage(df)

    rank_conditions = ~df['rank'].isin(EXCLUDED_RANKS)

    if len(kingdom_list) > 0:
        conditions = df.superkingdom.isin(kingdom_list)
        conditions |= df.kingdom.isin(kingdom_list)

        tdf = df[rank_conditions & conditions].nlargest(ntops, 'hits', keep="all")
    else:
        tdf = df[rank_conditions].nlargest(ntops, 'hits', keep="all")

    print(tdf.head(20))
    
    taxids = set(tdf.taxid.tolist())
    print("NCBI-TaxIDs:", taxids)

    key_for_taxid = '{}='.format(taxa_key)
    print("Key for parsing Tax IDs:", key_for_taxid)
    # ofile = open(fout, 'w')
    ofile = fout.open('w')

    num_seqs = 0
    start_time = time.time()
    
    with gzip.open(ref_fasta, 'rb') as file:
        try:
            is_selected = False
            for bline in file:
            
                # if line.startswith('>UniRef'):
                if bline[0] == 62:  # to find `>`
                    num_seqs += 1
                    line = bline.decode("utf-8")
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


# from pathlib import Path, PureWindowsPath

# prefix = "S1_NM0001_NMDC_MixedCommunities_R3_mgf"
# diamond_search_out = Path("Kaiko_volume/Kaiko_intermediate/" + prefix + "_diamond_search_output.dmd")
# kaiko_tally = Path("Kaiko_volume/Kaiko_intermediate/" + prefix + "_kaiko_prediction_top_taxa.csv")
# ncbi_taxa_folder = Path(PureWindowsPath("Kaiko_volume/Kaiko_stationary_files/ncbi_taxa").as_posix())
# nprot = '{:.5e}'.format(int(300000))
# kaiko_tally = Path("Kaiko_volume/Kaiko_intermediate/" + prefix + "_kaiko_prediction" + f'_top_taxa_nprot_{nprot}.csv')
# ref_fasta = Path(PureWindowsPath('Kaiko_volume/Kaiko_stationary_files/uniref100.fasta.gz').as_posix())
# kaiko_final_output = Path("Kaiko_volume/Kaiko_output/" + prefix + "_kaiko_output.fasta")



# aggregate_fasta(ref_fasta,
#                 kaiko_tally,
#                 kaiko_final_output,
#                 5,
#                 'TaxID',
#                 [])