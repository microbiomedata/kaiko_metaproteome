"""
this script is focused on minimizing the parsing and putting together
a comparison between the target sequence and what the two models,
kaiko and casanovo predicts
"""
from pyteomics import mgf, mztab
import pandas as pd
import numpy as np
import sys


def aggregate_kaiko_casanovo(mgf_path, kaiko_path, casanovo_path, output_path):
    # path below: path\to\mgf_file
    # mgf_file_path = sys.argv[1] # "C:\\Users\\leej179\\git\\kaiko_metaproteome\\Kaiko_volume\\Kaiko_input_files\\mgf_large_unit_test\\Biodiversity_A_cryptum_FeTSB_anaerobic_1_01Jun16_Pippin_16-03-39_unit_test_prop_0.005_seed_6969.mgf"
    # read a mgf file
    mgf_data = mgf.read(mgf_path, use_index=True)

    # extract the precursor info from mgf_data as a dataframe
    precursors = [d['params'] for d in mgf_data]
    df = pd.DataFrame.from_records(precursors)

    # create an index column
    df.reset_index(inplace=True)
    # select the first element in a pepmass tuple
    df['pepmass'] = df['pepmass'].apply(lambda x: x[0])
    # get a integer for a charge value from the charge list
    df['charge'] = df['charge'].apply(lambda x: int(x[0]))

    # path below: path\to\text_file
    # kaiko_output_path = sys.argv[2] # 'C:\\Users\\leej179\\git\\kaiko_metaproteome\\Kaiko_volume\\Kaiko_intermediate\\denovo_output\\mgf_large_unit_test\\Biodiversity_A_cryptum_FeTSB_anaerobic_1_01Jun16_Pippin_16-03-39_unit_test_prop_0.005_seed_6969_out.txt'
    # read a kaiko output file
    kaiko_df = pd.read_csv(kaiko_path, delimiter="\t")[['scan', 'target_seq', 'output_seq', 'exact_match']]
    for i in range(len(kaiko_df)):
        kaiko_seq = kaiko_df['output_seq'][i]
        target_seq = kaiko_df['target_seq'][i]
        kaiko_df.loc[i, 'output_seq'] = kaiko_seq.translate({ord(','): None})
        kaiko_df.loc[i, 'target_seq'] = target_seq.translate({ord(','): None})

    # path below: path\to\mztab_file
    # casanovo_output_path = sys.argv[3] # "C:\\Users\\leej179\\git\\casanovo\\output\\casanovo_test_1_output.mztab"
    # print(sys.argv)
    # read a casanovo output mztab file
    tables = mztab.MzTab(casanovo_path)
    # get a PSM table
    psms = tables.spectrum_match_table
    # extract the scan index from the `spectra_ref` column, e.g., ms_run[1]:index=25
    psms['scan_index'] = psms["spectra_ref"].str.extract("index=(\\d+)")[0]
    # convert the scan index as an object to integer
    psms['scan_index'] = psms['scan_index'].astype(int)
    # create the mztab dataframe with scan index and sequence
    casanovo_df = psms[['scan_index', 'sequence']]

    # modified the casanovo output to change the modified charge amino acid into simipler terms
    list_of_amino_acids = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M',
                           'F', 'P', 'O', 'S', 'U', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', 'J']
    for i in range(len(casanovo_df)):
        modified = []
        casanovo_sequence = []
        sequence = casanovo_df['sequence'][i+1]
        idx = 0
        for j, amino_acid in enumerate(sequence):
            # if sequence[i] in list_of_amino_acids:
            casanovo_sequence.append(amino_acid)
            if amino_acid not in list_of_amino_acids:
                idx = j
                modified.append(amino_acid)
                casanovo_sequence = casanovo_sequence[:idx-len(modified)+1] + ['mod'] + casanovo_sequence[idx+1:]
        casanovo_sequence = ','.join(casanovo_sequence)    
        casanovo_sequence = casanovo_sequence.translate({ord(','): None})
        casanovo_df.loc[i+1, 'sequence'] = casanovo_sequence
    # casanovo_sequence = [casanovo_df['sequence'][1][i] if casanovo_df['sequence'][1][i] in list_of_amino_acids else [casanovo_df['sequence'][1][i]] for i in range(len(casanovo_df['sequence'][1]))]
    

    # rename the column: sequnce to casanovo 
    casanovo_df = casanovo_df.rename(columns={'sequence': 'casanovo_peptide'})

    # merge the mgf df and kaiko df with scan number
    result = df.merge(kaiko_df, left_on="scans", right_on="scan")
    # merge the merged df with the casanovo df with scan index
    comparison_df = result.merge(casanovo_df, left_on="index", right_on="scan_index")
    # get the specific columns that are necessary
    comparison_df = comparison_df[[
        'scans',
        # 'pepmass',
        # 'charge',
        # 'rtinseconds',
        # 'seq',
        'target_seq',
        'output_seq',
        'casanovo_peptide',
        'exact_match'
    ]]
    # rename the columns to specific names
    comparison_df = comparison_df.rename(
        columns={
            'target_seq': 'original_peptide',
            # 'scans': 'scan_num',
            'output_seq': 'kaiko_peptide',
            'exact_match': 'kaiko_match'
            # 'casanovo_sequence': 'casanovo_seq'
        }
    )

    original_peptide = comparison_df.loc[:, 'original_peptide']
    casanovo_peptide = comparison_df.loc[:, 'casanovo_peptide']
    casanovo_exact = {'casanovo_match': []}
    for i, sequence in enumerate(original_peptide):
        if casanovo_peptide[i] == sequence:
            casanovo_exact['casanovo_match'].append('T')
        else:
            casanovo_exact['casanovo_match'].append('F')
    
    casanovo_exact_df = pd.DataFrame.from_records(casanovo_exact)
    # print(casanovo_exact_df)

    comparison_df = comparison_df.join(casanovo_exact_df, how='outer')
    kaiko_exact = comparison_df.loc[:, 'kaiko_match']
    for i, exact_match in enumerate(kaiko_exact):
        if exact_match == 1:
            comparison_df.loc[i, 'kaiko_match'] = 'T'
        else:
            comparison_df.loc[i, 'kaiko_match'] = 'F'

    # save the df into a text file called comparison_file.txt
    # output_path = 'C:\\Users\\leej179\\git\\kaiko_metaproteome\\kaiko_compare\\comparison_output\\comparison_file.txt'
    comparison_df.to_csv(output_path, index=None, sep="\t")
    return comparison_df


if __name__ == '__main__':
    mgf_path = sys.argv[1]
    kaiko_path = sys.argv[2]
    casanovo_path = sys.argv[3]
    output_path = sys.argv[4]
    aggregate_kaiko_casanovo(
        mgf_path=mgf_path,
        kaiko_path=kaiko_path,
        casanovo_path=casanovo_path,
        output_path=output_path
    )