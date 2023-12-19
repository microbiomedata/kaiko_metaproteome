import re
import os
import time
import requests
import pandas as pd
import datetime


from openpyxl import load_workbook
from pathlib import PureWindowsPath, Path
from download_proteomes import log_proteomes



referece_proteomes_table_path = Path(PureWindowsPath('proteomes_AND_proteome_type_1_2023_11_28.tsv'))
proteome_table = pd.read_csv(referece_proteomes_table_path, sep = "\t")
proteome_table = proteome_table.set_index('Organism Id')

log_file = Path(PureWindowsPath('database_log.xlsx'))



all_taxa_ids = proteome_table.index.to_list()
all_proteome_ids = proteome_table['Proteome Id'].to_list()
all_n_proteins = proteome_table['Protein count'].to_list()

sequence_pattern = re.compile(r'^[ACDEFGHIKLMNPQRSTVWYXBZJU]')
name_pattern = re.compile(r'^>.+ SV=[0-9]+\n$')

def check_integrity(proteome_path, n_proteins):
    actual_n_proteins = 0
    with proteome_path.open('r') as proteome_fasta:
        reading_protein = False
        empty_protein = False
        for line in proteome_fasta:
            if not (name_pattern.search(line) or sequence_pattern.search(line)):
                print('fudge')
            assert name_pattern.search(line) or sequence_pattern.search(line)
            if reading_protein:
                assert sequence_pattern.search(line) or not empty_protein
                empty_protein = False
                if name_pattern.search(line):
                    reading_protein = False
            if not reading_protein:
                reading_protein = True
                actual_n_proteins = actual_n_proteins + 1
                empty_protein = True
        assert actual_n_proteins == n_proteins and not empty_protein

taxa_ids = []
dates = []
proteome_paths = []
all_good = True
for taxa_id, proteome_id, n_proteins in zip(all_taxa_ids, all_proteome_ids, all_n_proteins):
    suffix = f'_taxaid_{taxa_id}'
    proteome_path = Path(PureWindowsPath(f'{proteome_id}{suffix}.fasta'))
    try:
        check_integrity(proteome_path, n_proteins)
    except:
        all_good = False
        taxa_ids = taxa_ids + [taxa_id]
        dates = dates + [datetime.datetime.now()]
        proteome_paths = proteome_paths + [proteome_path]
        print(f'taxa {taxa_id} with proteome id {proteome_id} did not pass integrity check')
        print(f'deleting...')
        # os.remove(proteome_path)
if all_good:
    print(f'All fasta passed integrity check!')
else:
    print("Logging fasta which failed integrity test")
    log_proteomes(taxa_ids, dates, proteome_paths, log_file)







