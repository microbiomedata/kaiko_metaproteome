import gzip
import time
import re
import gc
import indexed_gzip as igzip

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

# def sanity_check(ref_fasta_path, fout, taxa_key = 'TaxID', protein_key = '>UniRef100_'):

#     key_for_taxid = '{}='.format(taxa_key)
#     print("Key for parsing Tax IDs:", key_for_taxid)
#     start_time = time.time()
    
#     with fout.open('a') as building_index, igzip.IndexedGzipFile(str(ref_fasta_path), index_file = str(ref_fasta_igzip_index)) as file:
#         building_index.write('protein_id' + '\t' + 'common_or_title_taxid' + '\t' + 'position' + '\t' + 'seq_length' + '\n')
#         file.seek(160412426850)
#         universe = {}
#         num_seqs = 0
#         for bline in file:
#             # if line.startswith('>UniRef'):
#             if bline[0] == 62:  # to find `>`
#                 line = bline.decode("utf-8")
#                 taxid = line.split(key_for_taxid)[1].split(' ')[0]
#                 if (taxid == "UPI001FFFB27D"):
#                     print(bline)



## Keeping track of protein locations in fasta.gz
def make_protein_index(ref_fasta_path, fout, taxa_key = 'TaxID', protein_key = '>UniRef100_'):

    key_for_taxid = '{}='.format(taxa_key)
    print("Key for parsing Tax IDs:", key_for_taxid)
    start_time = time.time()
    
    with fout.open('a') as building_index, igzip.IndexedGzipFile(str(ref_fasta_path), index_file = str(ref_fasta_igzip_index)) as file:
        building_index.write('protein_id' + '\t' + 'common_or_title_taxid' + '\t' + 'position' + '\t' + 'seq_length' + '\n')
        was_reading_protein = False
        current_position = 0
        protein_position = 0
        universe = {}
        num_seqs = 0
        for bline in file:
            # if line.startswith('>UniRef'):
            if bline[0] == 62:  # to find `>`
                if (num_seqs % 100000) == 0:
                    for k, v in universe.items():
                        building_index.write(f'{k}\t{v}\n')
                    print("{}M sequences has been parsed. {:.1f}min".format((num_seqs//1e5)/10, (time.time()-start_time)/60))
                    universe = {}
                if was_reading_protein:
                    universe[protein_id] = taxid + '\t' + str(protein_position) + '\t' + str(seq_length)
                protein_position = current_position
                seq_length = 0
                line = bline.decode("utf-8")
                taxid = line.split(key_for_taxid)[1].split(' ')[0]
                protein_id = line.split(protein_key)[1].split(' ')[0]
                if taxid in ["", "N/A"]:
                    was_reading_protein = False
                else:
                    was_reading_protein = True
                    num_seqs += 1
            else:
                if was_reading_protein:
                    seq_length += (len(bline) - 1)
            current_position = file.tell()
        ## Write last proteins
        for k, v in universe.items():
            building_index.write(f'taxid_{k}\t{v}\n')
        print("Writing last {} sequences. {}M sequences has been parsed. {:.1f}min".format(len(universe.keys()), (num_seqs//1e5)/10, (time.time()-start_time)/60))
        building_index.flush()


def load_protein_index(protein_index_path):
    protein_index = protein_index_path.open('r')
    universe = {}
    header = protein_index.readline()
    start_time = time.time()
    nlines = 0
    for line in protein_index:
        nlines = nlines + 1
        line = line.split('\t')
        protein = line[0]
        protein_pos = line[2]
        protein_len = line[3].split('\n')[0]
        universe[protein] = protein_pos + '_' + protein_len
        if nlines % 1e6 == 0:
            print("{}M sequences have been loaded. {:.1f}min".format(nlines//1e6, (time.time()-start_time)/60))

    print("Done!")
    print(time.perf_counter() - start_time)
    return(universe)

def rank_to_lineage(df):
    coverage = {}
    for c in ['species','genus','family','order','class','phylum','kingdom','superkingdom']:
        idx = df['rank']==c
        df.loc[idx, c] = df.loc[idx, 'tax_name']
    for c in ['species','genus','family','order','class','phylum','superkingdom']:
        coverage[c] = 100*df[c].dropna().shape[0]/df.shape[0]
        print("{}: {:.2f}%".format(c, coverage[c]))
    return df, coverage

## Part 1 of the index process. 
## This makes a partial index. Takes a few hours. A complete index in one pass would take too much time. 
## We are keeping track of the location of each protein in the FASTA file, and the taxa it belongs to.
def gather_taxa_stats_1(taxa_member_path, protein_index_path, fout, chunksize = 1000000):
    ofile = fout.open('w')
    start_time = time.time()
    universe = {}
    exceptions = []
    protein_index = load_protein_index(protein_index_path)

    for chunk in pd.read_csv(taxa_member_path, chunksize=chunksize):
        chunk['uid'] = [x.split('_')[1] for x in chunk['uid']]
        for protein_id, taxa_list in zip(chunk['uid'], chunk['members']):
            if protein_id in protein_index.keys():
                protein_position, protein_len = protein_index[protein_id].split('_')[0], protein_index[protein_id].split('_')[1]
                for taxa in set(taxa_list.split(':')):
                    if taxa in universe.keys():
                        universe[taxa] = universe[taxa] + f':{protein_id};{protein_position};{protein_len}'
                    else:
                        universe[taxa] = f'{protein_id};{protein_position};{protein_len}'
            else:
                print("Found exception. Protein " + protein_id + " not found in protein index")
                exceptions = exceptions + [protein_id]
        for k, v in universe.items():
            ofile.write(f'taxid_{k}\t{v}\n')
        universe = {}
    ofile.close()

    exceptions_path = Path(PureWindowsPath('Kaiko_volume/Kaiko_stationary_files/protein_index_exceptions.txt'))
    exceptions_file = exceptions_path.open('w')
    for exception in exceptions:
        exceptions_file.write(f'{exception}\n')
    exceptions_file.close()


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
                universe[taxid + '_positions'] = 'uniref100.fasta;;'
                ## first coordinate is protein number, second is total protein lenth (fasta size)
                universe[taxid + '_size'] = (0, 0)
            ## info_list looks like this -> Q6GZX4;0;256:Q6GZX2;1641;438:Q6GZX1;2179;60:Q6GZW8;3469;128:Q6GZW6;4734;948. 
            ## : separates proteins, Q6GZX4;0;256 refers to protein_name;protein_position;protein_length
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

## The index is very large on its own.  
## To speed up parsing and allow for fast and dynamic making of FASTA files,
## we make an index for the index. Yes.
def gather_taxa_stats_3(parse2_output, parse3_output):
    index = 0
    pos = 0
    with parse2_output.open('r') as file:
        output = parse3_output.open("w")
        line = file.readline()
        while line:
            if index % 3 == 0 or index % 3 == 1:
                taxid = line.split('\t')[0]
                output.write(f'{taxid}\t{pos}\n')
            ## have to add +1 to account for the \n character.
            index = index + 1
            pos = file.tell()
            line = file.readline()

## Get taxa stats in an easy to open table
def gather_taxa_stats_4(parse2_output, parse3_output, stats_output):
    index_file = parse2_output.open('r')
    output = stats_output.open('w')
    output.write('taxid\tn_protein\tn_AA\n')
    with parse3_output.open('r') as file:
        for line in file:
            line_start = line.split('\t')[0].split('_')
            if len(line_start) > 2:
                position = int(line.split('\t')[1].split('\n')[0])
                index_file.seek(position)
                index_file.readline()
                xx = index_file.readline()
                taxid = line_start[1]
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



def write_single_protein(database_file, pos, output_fasta):
    database_file.seek(pos)
    
    line = database_file.readline()
    line = line.decode("utf-8")

    reading_protein = True
    while reading_protein:
        output_fasta.write(line)
        line = database_file.readline()
        if line[0] == 62:  # to find `>`
            reading_protein = False
        else:
            line = line.decode("utf-8")
    


parse1_output = Path(PureWindowsPath('Kaiko_volume/Kaiko_stationary_files/uniref100_protein_index.txt'))
parse2_output = Path(PureWindowsPath('Kaiko_volume/Kaiko_stationary_files/uniref100_index_partial.txt'))
parse3_output = Path(PureWindowsPath('Kaiko_volume/Kaiko_stationary_files/uniref100_index_s.txt'))
parse4_output = Path(PureWindowsPath('Kaiko_volume/Kaiko_stationary_files/uniref100_index.txt'))
protein_index_path = Path(PureWindowsPath('Kaiko_volume/Kaiko_stationary_files/uniref100_protein_index.txt'))
stats_output = Path(PureWindowsPath('Kaiko_volume/Kaiko_stationary_files/uniref100_stats.txt'))
ref_fasta_path = Path(PureWindowsPath('Kaiko_volume/Kaiko_stationary_files/uniref100.fasta.gz'))
ref_fasta_igzip_index = Path(PureWindowsPath('Kaiko_volume/Kaiko_stationary_files/uniref100_fasta_gzindex.gzidx'))
taxa_member_path = Path(PureWindowsPath('Kaiko_volume/Kaiko_stationary_files/ncbi_taxa/uniref100_member_taxa_tbl.csv'))
# ref_fasta = Path(PureWindowsPath('D:/uniref100.fasta.gz'))
# ref_fasta_igzip_index = Path(PureWindowsPath('D:/uniref100_fasta_gzindex.gzidx'))
ncbi_taxa_folder = Path(PureWindowsPath('Kaiko_volume/Kaiko_stationary_files/ncbi_taxa').as_posix())

# output_fasta_path = Path(PureWindowsPath('Kaiko_volume/testing.FASTA'))
# output_fasta = open(output_fasta_path,  "a")

# with igzip.IndexedGzipFile(str(ref_fasta_path), index_file = str(ref_fasta_igzip_index)) as file:
#     write_single_protein(file, 10, output_fasta)

# load_index_partial(parse2_output)

# make_protein_index(ref_fasta_path, parse1_output)
# sanity_check(ref_fasta_path, parse1_output)
# load_protein_index(protein_index_path)
# gather_taxa_stats_1(taxa_member_path, protein_index_path, parse2_output)

# gc.collect()
# gather_taxa_stats_2(parse2_output, parse4_output)

# gc.collect()
# gather_taxa_stats_3(parse4_output, parse3_output)

# gather_taxa_stats_4(parse4_output, parse3_output, stats_output)

gather_taxa_stats_5(stats_output, ncbi_taxa_folder)