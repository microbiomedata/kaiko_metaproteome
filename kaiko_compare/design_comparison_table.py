"""
this script is focused on minimizing the parsing and putting together
a comparison between the target sequence and what the two models,
kaiko and casanovo predicts
"""
from pyteomics import mgf
from pyteomics import mztab
import pandas as pd
import numpy as np
import sys

# path below: path\to\mgf_file
mgf_file_path = sys.argv[1] # "C:\\Users\\leej179\\git\\kaiko_metaproteome\\Kaiko_volume\\Kaiko_input_files\\mgf_large_unit_test\\Biodiversity_A_cryptum_FeTSB_anaerobic_1_01Jun16_Pippin_16-03-39_unit_test_prop_0.005_seed_6969.mgf"
# read a mgf file
mgf_data = mgf.read(mgf_file_path, use_index=True)

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
kaiko_output_path = sys.argv[2] # 'C:\\Users\\leej179\\git\\kaiko_metaproteome\\Kaiko_volume\\Kaiko_intermediate\\denovo_output\\mgf_large_unit_test\\Biodiversity_A_cryptum_FeTSB_anaerobic_1_01Jun16_Pippin_16-03-39_unit_test_prop_0.005_seed_6969_out.txt'
# read a kaiko output file
kaiko_df = pd.read_csv(kaiko_output_path, delimiter="\t")

# path below: path\to\mztab_file
casanovo_output_path = sys.argv[3] # "C:\\Users\\leej179\\git\\casanovo\\output\\casanovo_test_1_output.mztab"
print(sys.argv)
# read a casanovo output mztab file
tables = mztab.MzTab(casanovo_output_path)
# get a PSM table
psms = tables.spectrum_match_table
# extract the scan index from the `spectra_ref` column, e.g., ms_run[1]:index=25
psms['scan_index'] = psms["spectra_ref"].str.extract("index=(\\d+)")[0]
# convert the scan index as an object to integer
psms['scan_index'] = psms['scan_index'].astype(int)
# create the mztab dataframe with scan index and sequence
casanovo_df = psms[['scan_index', 'sequence']]
# rename the column: sequnce to casanovo 
casanovo_df = casanovo_df.rename(columns={'sequence': 'casanovo_sequence'})

# merge the mgf df and kaiko df with scan number
result = df.merge(kaiko_df, left_on="scans", right_on="scan")
# merge the merged df with the casanovo df with scan index
comparison_df = result.merge(casanovo_df, left_on="index", right_on="scan_index")
# get the specific columns that are necessary
comparison_df = comparison_df[[
    'scans',
    'pepmass',
    'charge',
    'rtinseconds',
    'seq',
    'target_seq',
    'output_seq',
    'casanovo_sequence'
]]
# rename the columns to specific names
comparison_df = comparison_df.rename(
    columns={
        'scans': 'scan_num',
        'output_seq': 'kaiko_seq',
        'casanovo_sequence': 'casanovo_seq'
    }
)


# save the df into a text file called comparison_file.txt
output_path = 'C:\\Users\\leej179\\git\\kaiko_metaproteome\\kaiko_compare\\comparison_output\\comparison_file.txt'
comparison_df.to_csv(output_path, sep="\t")
