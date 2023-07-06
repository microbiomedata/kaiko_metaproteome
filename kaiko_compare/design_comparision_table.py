from pyteomics import mgf
from pyteomics import mztab
import pandas as pd
# import numpy as np


f = mgf.read("C:\\Users\\leej179\\git\\casanovo\\sample_data\\mgf_test_file.mgf", use_index=True)
sequence = {}
for i in range(0, 40):
    #converted the integer index to string || can change
    sequence[i] = i
    sequence[i] = f[i]['params']

df = pd.DataFrame.from_dict({i: sequence[i]
                            for i in sequence.keys()},
                            orient='index')
df.reset_index(inplace=True)
df['pepmass'] = df['pepmass'].apply(lambda x: x[0])
df['charge'] = df['charge'].apply(lambda x: int(x[0]))


data = pd.read_csv('C:\\Users\\leej179\\git\\casanovo\\output\\Biodiversity_A_cryptum_FeTSB_anaerobic_1_01Jun16_Pippin_16-03-39_unit_test_prop_0.005_seed_6969_out.txt',
                   delimiter="\t")
data.reset_index(inplace=True)
# print(data)

result = pd.merge(df, data, on="index")
print(result)

