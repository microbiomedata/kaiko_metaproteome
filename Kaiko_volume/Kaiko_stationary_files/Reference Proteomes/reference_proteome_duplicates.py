import re
import os
import time
import requests
import pandas as pd
import datetime

from collections import Counter
from pathlib import PureWindowsPath, Path

name_pattern = re.compile(r'^>.+ SV=[0-9]+\n$')

referece_proteomes_table_path = Path(PureWindowsPath('proteomes_AND_proteome_type_1_2023_11_28.tsv'))
proteome_table = pd.read_csv(referece_proteomes_table_path, sep = "\t")
proteome_table = proteome_table.set_index('Organism Id')

all_taxa_ids = proteome_table.index.to_list()
all_proteome_ids = proteome_table['Proteome Id'].to_list()
all_n_proteins = proteome_table['Protein count'].to_list()

depth_of_folders = 2
depth_of_files = 3
def compress_protein(protein_seq):
    # A C D E F G H I K L M N P Q R S T V W Y X B Z J U 
    protein_letters = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y','X','B','Z','J','U']
    #strip any newline characters
    protein_seq = re.sub("\n", "", protein_seq)
    
    all_counts = Counter(protein_seq)
    folder_name = ''
    file_name = ''
    index = 1
    most_common = all_counts.most_common(depth_of_folders + depth_of_files)
    for pair in most_common[0:depth_of_folders]:
        folder_name = f'{folder_name}max{index}_{pair[0]}/'
        index = index + 1
    index = 1
    for pair in most_common:
        file_name = f'{file_name}max{index}_{pair[0]}_'
        index = index + 1

    compression = 'counts'
    for amino_acid in protein_letters:
        # order is A C D E F G H I K L M N P Q R S T V W Y X B Z J U 
        if amino_acid in all_counts.keys():
            ac_count = all_counts[amino_acid]
        else:
            ac_count = 0
        compression = f'{compression}_{ac_count}'
    
    return(f'{folder_name}/{file_name}.txt', compression)


def read_compression(compressed_path):
    start_time = time.time()
    with compressed_path.open() as signatures:
        signatures.readline()
        # for lion
        all_sigs = dict()
        for line in signatures:
            signature = line.split('\t')[3].split('\n')[0]
            if signature in all_sigs.keys():
                all_sigs[signature] = all_sigs[signature] + [line]
            else:
                all_sigs[signature] = [line]
    print(time.time() - start_time)
    return all_sigs

compressed_path = Path(PureWindowsPath('compressed/max1_A/max2_L/max1_A_max2_L_max3_G_max4_V_max5_R_.txt'))
xx = read_compression(compressed_path)

# compressed = Path(PureWindowsPath('compressed/'))
# M = 0
# chosen_file = ''
# for folder in compressed.glob('**/'):
#     for file in folder.glob('*.txt'):
#         file_size = os.path.getsize(file)
#         if file_size > M:
#             M = file_size
#             chosen_file = file
#             print(M)
#             print(chosen_file)

# print(chosen_file)
# print(M)

# for taxa_id, proteome_id, n_proteins in zip(all_taxa_ids, all_proteome_ids, all_n_proteins):
#     suffix = f'_taxaid_{taxa_id}'
#     proteome_path = Path(PureWindowsPath(f'{proteome_id}{suffix}.fasta'))

#     actual_n_proteins = 0
#     with proteome_path.open('r') as proteome_fasta:
#         reading_protein = False
#         line = proteome_fasta.readline()
#         while line:
#             if name_pattern.search(line):
#                 if reading_protein:
#                     file_path, compression = compress_protein(protein_seq)
#                     compressed_path = Path(PureWindowsPath('compressed')) / file_path
#                     compressed_path.parent.mkdir(parents=True, exist_ok=True)
#                     if not compressed_path.exists():
#                         with compressed_path.open('w') as fp:
#                             fp.write('taxa_id\tprotein_name\tposition\tsignature(A_C_D_E_F_G_H_I_K_L_M_N_P_Q_R_S_T_V_W_Y_X_B_Z_J_U)\n')
#                     with compressed_path.open('a') as fp:
#                         fp.write(f'{taxa_id}\t{protein_name}\t{sequence_pos}\t{compression}\n')
#                 sequence_pos = proteome_fasta.tell()
#                 protein_name = line.split('|')[1]
#                 protein_seq = ''
#                 reading_protein = True
#             else:
#                 protein_seq = protein_seq + line
#             line = proteome_fasta.readline()
#     print(taxa_id)





