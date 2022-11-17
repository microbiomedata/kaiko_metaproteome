import gzip
import time
import pandas as pd
import numpy as np
import argparse
import os

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

###############################################################
# parser = argparse.ArgumentParser()

# parser.add_argument(
#     '--taxafile', type=str, required=True,
#     help='top-rank taxa file.')

# parser.add_argument(
#     '--ref_fasta', type=str, required=True,
#     help='reference fasta file')

# parser.add_argument(
#     '--fout', type=str, required=True,
#     help='output fasta file.')

# parser.add_argument(
#     '--key', type=str, default="TaxID",
#     help='key value for the taxa ID. `TaxID` for UniRef fasta and `OX` for UniProtKB fasta')

# parser.add_argument('-l','--list', nargs='+', required=False)

# # parser.add_argument('-r','--rank', nargs='+', required=True)

# parser.add_argument(
#     '--ntops', type=int, default=5, required=True,
#     help='top n taxa. -1 for all')

# FLAGS = parser.parse_args()
# print(FLAGS)

###############################################################



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
    
    taxids = set(tdf.tax_id.tolist())
    print("NCBI-TaxIDs:", taxids)

    key_for_taxid = '{}='.format(taxa_key)
    print("Key for parsing Tax IDs:", key_for_taxid)
    ofile = open(fout, 'w')

    num_seqs = 0
    start_time = time.time()
    
    print(ref_fasta)
    print("\n")
    print(os.getcwd())
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
