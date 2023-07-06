from pyteomics import mztab
import numpy as np
import os

tables = mztab.MzTab("C:\\Users\\leej179\\git\\casanovo_output_trial_1_comp.mztab")
psms = tables.spectrum_match_table
sequence = psms["sequence"]
index = psms["spectra_ref"]
# for i in range(40):
#     psms.loc(, ["spectra_ref"]) = i
# print(psms[["sequence", "spectra_ref"]])

# print(sequence)
# mztab._MzTabParserBase("sequence", ["sequence"])
path = r"C:\\Users\\leej179\\git\\mztab_parsed.txt"
# assert os.path.isfile(path)
# with open(path, "r") as f:
#     pass

with open(path, 'a') as f:
    df_string = sequence.to_string(header=False, index=False)
    f.write(df_string)