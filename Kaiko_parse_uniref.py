import gzip
import time
import re
import gc

import pandas as pd
import numpy as np

from pathlib import Path, PureWindowsPath

## Moving this function here, as the lineage info can be included in the taxa stats file.
## taxon info
def read_ncbi_taxa_lineage(rankedlineage_file, nodes_file):
    taxon_info = pd.read_csv(rankedlineage_file, sep='|', header=None,
                             names=['tax_id','tax_name','species','genus','family','order',
                                    'class','phylum','kingdom','superkingdom','null'])
    taxon_info = taxon_info.replace('\t\t', np.nan)
    for col in ['tax_name','species','genus','family','order','class','phylum','kingdom','superkingdom']:
        taxon_info[col] = taxon_info[col].str.extract(r'^\t(.*)\t$')[0].tolist()
        # taxon_info[col] = taxon_info[col].str.extract(r'^\t(.*)\t$').tolist()
    del taxon_info['null']
    print('  taxon_info:', taxon_info.shape, taxon_info.tax_id.drop_duplicates().shape)
    
    rank = pd.read_csv(nodes_file, sep='|', usecols=[0,2], header=None,
                     names=['tax_id','rank'])
    rank = rank.replace('\t\t', np.nan)
    rank['rank'] = rank['rank'].str.extract(r'^\t(.*)\t$')[0].tolist()
    # rank['rank'] = rank['rank'].str.extract(r'^\t(.*)\t$').tolist()
    print('  rank:', rank.shape, rank.tax_id.drop_duplicates().shape)
    
    final_df = taxon_info.merge(rank, left_on='tax_id', right_on='tax_id')
    print('  final_df:', final_df.shape, final_df.tax_id.drop_duplicates().shape)
    return final_df.set_index('tax_id')

## Part 1 of the index process. 
## This makes a partial index. Takes a few hours. A complete index in one pass would take too much time. 
## We are keeping track of the location of each protein in thge FASTA file, and the taxa it belongs to.
def gather_taxa_stats_1(ref_fasta, fout, taxa_key = 'TaxID', protein_key = '>UniRef100_'):
    """
    python -u ExtractProteinSequenceByTaxID.py --ref_fasta uniref100.fasta.gz --taxafile Weintraub_kaiko_25p_UniRef100_toptaxa.csv --fout Weintraub_kaiko_25p_UniRef100_Bacteria_top100.fasta --ntops 100 -l Bacteria --key TaxID
    """

    key_for_taxid = '{}='.format(taxa_key)
    print("Key for parsing Tax IDs:", key_for_taxid)
    # ofile = open(fout, 'w')
    ofile = fout.open('w')
    ofile.write('protein_id' + '\t' + 'taxid' + '\t' + 'seq_length' + '\n')

    num_seqs = 0
    start_time = time.time()
    
    with gzip.open(ref_fasta, 'rb') as file:
        try:
            reading_protein = False
            current_position = 0
            protein_position = 0
            output = open("UniRef100_index.txt", "w")
            universe = {}
            for bline in file:
            
                # if line.startswith('>UniRef'):
                if bline[0] == 62:  # to find `>`
                    if reading_protein:
                        if taxid in universe.keys():
                            universe[taxid] = universe[taxid] + f':{protein_id};{protein_position};{seq_length}'
                        else:
                            universe[taxid] = f'{protein_id};{protein_position};{seq_length}'
                    if (num_seqs % 100000) == 0:
                        for k, v in universe.items():
                            output.write(f'taxid_{k}\t{v}\n')
                        print("{}M sequences has been parsed. {:.1f}min".format(num_seqs//1e6, (time.time()-start_time)/60))
                        universe = {}
                    protein_position = current_position
                    seq_length = 0
                    num_seqs += 1
                    line = bline.decode("utf-8")
                    taxid = line.split(key_for_taxid)[1].split(' ')[0]
                    protein_id = line.split(protein_key)[1].split(' ')[0]
                    if taxid in ["", "N/A"]:
                        reading_protein = False
                    else:
                        reading_protein = True
                else:
                    if reading_protein:
                        seq_length += (len(bline) - 1)
                current_position += len(bline)
            
        except Exception as e:
            print(line, e)

    output.close()

def rank_to_lineage(df):
    coverage = {}
    for c in ['species','genus','family','order','class','phylum','kingdom','superkingdom']:
        idx = df['rank']==c
        df.loc[idx, c] = df.loc[idx, 'tax_name']
    for c in ['species','genus','family','order','class','phylum','superkingdom']:
        coverage[c] = 100*df[c].dropna().shape[0]/df.shape[0]
        print("{}: {:.2f}%".format(c, coverage[c]))
    return df, coverage


## Part two of the index process. This makes the complete index from the partial index above. Takes about ~10 mins.
## The result is:
# taxid    protein1;protein2;...
# taxid_position    pos1;pos2;...
# taxid_size    (n_proteins, n_AA)
def gather_taxa_stats_2(parse1_output, parse2_output):
    start_time = time.time()
    universe = dict()
    index = 0
    with parse1_output.open('r') as file:
        for line in file:
            index = index + 1
            taxid = line.split('_')[1].split('\t')[0]
            if not taxid in universe.keys():
                universe[taxid] = 'uniref100;;'
                ## first coordinate is protein number, second is total protein lenth (fasta size)
                universe[taxid + '_positions'] = 'uniref100.fasta;;'
                universe[taxid + '_size'] = (0, 0)
            info_list = re.split('[: ;]', line.split('_')[1].split('\t')[1])

            # protein is index of the form 3n (start at zero)
            # protein position in fasta is 3n + 1  
            # protein length is 3n + 2 
            num_proteins = len(info_list)//3
            new_proteins = [info_list[3*i] for i in range(num_proteins)]
            prot_pos = [info_list[3*i+1] for i in range(num_proteins)]
            new_len = sum([int(info_list[3*i+2]) for i in range(num_proteins)])
            new_proteins = ';'.join(new_proteins)
            universe[taxid] = universe[taxid] + new_proteins + ';'
            universe[taxid + '_size'] = (universe[taxid + '_size'][0] + num_proteins, universe[taxid + '_size'][1] + new_len)
            universe[taxid + '_positions'] = universe[taxid + '_positions'] + ';'.join(prot_pos) + ';'
            if (index % 100000) == 0:
                print("{} total partial taxa found. {:.1f}min".format(len(universe.keys())//3, (time.time()-start_time)/60))
    with parse2_output.open('w') as file:
        for k, v in universe.items():
            file.write(f'taxid_{k}\t{v}\n')

## The index is quite large on its own.  
## To speed up parsing and allow for fast and dynamic making of FASTA files,
## we make an index for the index. Yes.
def gather_taxa_stats_3(parse2_output, parse3_output):
    index = 0
    pos = 0
    with parse2_output.open('r') as file:
        output = parse3_output.open("w")
        for line in file:
            if index % 3 == 0:
                taxid = line.split('\t')[0]
                output.write(f'{taxid}\t{pos}\n')
            ## have to add +1 to account for the \n character.
            pos = pos + len(line) + 1

## Get taxa stats in an easy to open table
def gather_taxa_stats_4(parse2_output, parse3_output, stats_output):
    index_file = parse2_output.open('r')
    output = stats_output.open('w')
    output.write('taxid\tn_protein\tn_AA\n')
    with parse3_output.open('r') as file:
        for line in file:
            line_start = line.split('\t')[0].split('_')
            if len(line_start) > 2:
                taxid = line_start[1]
                if line_start[2] == "size":
                    pos = int(line.split('\t')[1])
                    index_file.seek(pos)
                    xx = index_file.readline()
                    xx = xx.split('\t')[1]
                    n_protein = re.sub(r'\(([0-9]+), ([0-9]+)\)\n$', "\\1", xx)
                    n_AA = re.sub(r'\(([0-9]+), ([0-9]+)\)\n$', "\\2", xx)
                    output.write(f'{taxid}\t{n_protein}\t{n_AA}\n')

def gather_taxa_stats_5(stats_output, ncbi_taxa_folder):
    ############################################################
    # retrieve taxa lineage info such as taxa rank and lineage
    ############################################################
    print("Retrieving taxa lineage info...")
    taxa_table = read_ncbi_taxa_lineage(ncbi_taxa_folder / 'rankedlineage.dmp', 
                                        ncbi_taxa_folder / 'nodes.dmp')
    taxa_stats = pd.read_csv(stats_output, sep = '\t')
    print(len(taxa_stats))
    taxa_stats = taxa_stats.merge(taxa_table, left_on="taxid", right_on="tax_id", how="left")
    taxa_stats = taxa_stats[['taxid', 'tax_name', 'rank', 'n_protein', 'n_AA', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom', 'superkingdom']]
    taxa_stats.to_csv(stats_output, sep = '\t', index = False)
    


parse1_output = Path(PureWindowsPath('Kaiko_volume/Kaiko_stationary_files/uniref100_index_intermediate.txt'))
parse2_output = Path(PureWindowsPath('Kaiko_volume/Kaiko_stationary_files/uniref100_index.txt'))
parse3_output = Path(PureWindowsPath('Kaiko_volume/Kaiko_stationary_files/uniref100_index_s.txt'))
stats_output = Path(PureWindowsPath('Kaiko_volume/Kaiko_stationary_files/uniref100_stats.txt'))
ref_fasta = Path(PureWindowsPath('Kaiko_volume/Kaiko_stationary_files/uniref100.fasta.gz'))
ncbi_taxa_folder = Path(PureWindowsPath('Kaiko_volume/Kaiko_stationary_files/ncbi_taxa').as_posix())

# gather_taxa_stats_1(ref_fasta, parse1_output)

# gc.collect()
# gather_taxa_stats_2(parse1_output, parse2_output)

# gc.collect()
# gather_taxa_stats_3(parse2_output, parse3_output)

gather_taxa_stats_4(parse2_output, parse3_output, stats_output)

gather_taxa_stats_5(stats_output, ncbi_taxa_folder)