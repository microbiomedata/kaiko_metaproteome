import glob
from pyteomics import mztab
import os

mgf_files = glob.glob("C:\\Users\\leej179\\git\\kaiko_metaproteome\\Kaiko_volume\\Kaiko_input_files\\mgf_large_unit_test\\*.mgf")
kaiko_files = glob.glob("C:\\Users\\leej179\\git\\kaiko_metaproteome_main\\Kaiko_volume\\Kaiko_intermediate\\denovo_output\\mgf_large_unit_test\\*.txt")
casanovo_files = glob.glob("C:\\Users\\leej179\\git\\casanovo\\output\\*.mztab")
mgf_dict = {}

def get_mgf_names(mgf_files):
    for mgf_path in mgf_files:
        mgf_name = os.path.basename(mgf_path)
        # print(os.path.splitext(mgf_name)[0])
        mgf_dict['mgf_names'] = mgf_name
    print(mgf_dict)

def get_kaiko_names(kaiko_files):
    for kaiko_path in kaiko_files:
        kaiko_name = os.path.basename(kaiko_path)
        # print(os.path.splitext(kaiko_name)[0])

def get_casanovo_names(casanovo_files):
    for casanovo_path in casanovo_files:
        tables = mztab.MzTab(casanovo_path)
        casanovo_name_path = tables.metadata['ms_run[1]-location'][8:]
        casanovo_name = os.path.basename(casanovo_name_path)
        # print(os.path.splitext(casanovo_name)[0])

if __name__ == '__main__':
    # if get_mgf_names(mgf_files) == get_kaiko_names(kaiko_files) == get_casanovo_names(casanovo_files):
    #     print('yes')
    get_mgf_names(mgf_files)
    # get_kaiko_names(kaiko_files)
    # get_casanovo_names(casanovo_files)