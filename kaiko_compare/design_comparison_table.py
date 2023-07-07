from pyteomics import mgf
from pyteomics import mztab
import pandas as pd
import numpy as np


f = mgf.read("C:\\Users\\leej179\\git\\kaiko_metaproteome\\Kaiko_volume\\Kaiko_input_files\\mgf_large_unit_test\\Biodiversity_A_cryptum_FeTSB_anaerobic_1_01Jun16_Pippin_16-03-39_unit_test_prop_0.005_seed_6969.mgf", use_index=True)
sequence = {}
for i in range(0, 40):
    sequence[i] = i
    sequence[i] = f[i]['params']

df = pd.DataFrame.from_dict({i: sequence[i]
                            for i in sequence.keys()},
                            orient='index')
df.reset_index(inplace=True)
df['pepmass'] = df['pepmass'].apply(lambda x: x[0])
df['charge'] = df['charge'].apply(lambda x: int(x[0]))


kaiko_df = pd.read_csv('C:\\Users\\leej179\\git\\kaiko_metaproteome\\Kaiko_volume\\Kaiko_intermediate\\denovo_output\\mgf_large_unit_test\\Biodiversity_A_cryptum_FeTSB_anaerobic_1_01Jun16_Pippin_16-03-39_unit_test_prop_0.005_seed_6969_out.txt',
                   delimiter="\t")
kaiko_df.reset_index(inplace=True)
# print(data)
tables = mztab.MzTab("C:\\Users\\leej179\\git\\casanovo\\output\\casanovo_test_1_output.mztab")
psms = tables.spectrum_match_table
# index = index.apply(lambda x: x[-2:-1])
psms['scan_index'] = psms["spectra_ref"].str.extract("index=(\\d+)")[0]
psms['scan_index'] = psms['scan_index'].astype(int)
mztab_df = psms[['scan_index', 'sequence']]
casanovo_df = mztab_df.rename(columns={'sequence': 'casanovo_sequence'})

result = pd.merge(df, kaiko_df, on="index")
result = result.rename(columns={'index': 'scan_index'})
comparison_df = pd.merge(result, casanovo_df, on="scan_index")
comparison_df = comparison_df[['scans', 'pepmass', 'charge', 'rtinseconds', 'target_seq', 'output_seq', 'casanovo_sequence']]
comparison_df = comparison_df.rename(columns={'scans': 'scan_num', 'output_seq': 'kaiko_seq', 'casanovo_sequence': 'casanovo_seq'})


comparison_np = comparison_df.to_numpy()
np.savetxt('C:\\Users\\leej179\\git\\kaiko_metaproteome\\kaiko_compare\\comparison_output\\comparison_file.txt', comparison_np, fmt="%s")
