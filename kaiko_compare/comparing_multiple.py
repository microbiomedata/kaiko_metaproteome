"""
this is to combine multiple outputs from kaiko and casanovo in directories.

Inputs:
- mgf filenames: input1.mgf, input2.mgf, ...
- kaiko output filenames: input1_out.txt, input2_out.txt, ...
- casanovo output filenames: casanovo_20230710143300.mztab, casanovo_20230710144300.mztab, ...

Result:
[(input1.mgf, input1_out.txt, casanovo_20230710143300.mztab), ...]
"""
import glob
import pandas as pd
from pyteomics import mztab
import os
from design_comparison_table import aggregate_kaiko_casanovo

mgf_files = glob.glob("C:\\Users\\leej179\\git\\kaiko_metaproteome\\Kaiko_volume\\Kaiko_input_files\\mgf_large_unit_test\\*.mgf")
kaiko_files = glob.glob("C:\\Users\\leej179\\git\\kaiko_metaproteome_main\\Kaiko_volume\\Kaiko_intermediate\\denovo_output\\mgf_large_unit_test\\*.txt")
casanovo_files = glob.glob("C:\\Users\\leej179\\git\\casanovo\\output\\*.mztab")
output_path = 'C:\\Users\\leej179\\git\\kaiko_metaproteome\\kaiko_compare\\comparison_output\\comparison.tsv'


def get_mgf_path_by_filename(mgf_paths):
    """get a dict for map a mgf filename to mgf path 
    Args:
        mgf_paths: list
            a list of mgf file paths
    Return:
        mgf_dict: dict
            {
                'input1': 'path/to/input1.mgf',
                'input2': 'path/to/input2.mgf',
                ...
            }
    """
    mgf_dict = {}
    for mgf_path in mgf_paths:
        mgf_dict[os.path.splitext(os.path.basename(mgf_path))[0]] = mgf_path
    # mgf_dict['mgf_names'] = [os.path.splitext(os.path.basename(mgf_path))[0] for mgf_path in mgf_files]
    return mgf_dict


def get_kaiko_path_by_mgf_filename(kaiko_paths):
    """map a mgf filename into a corresponding kaiko output path
    Args:
        mgf_paths: list
            a list of mgf file paths
    Return:
        mgf_dict: dict
            {
                'input1': 'path/to/input1_out.txt',
                'input2': 'path/to/input2_out.txt',
                ...
            }
    """
    kaiko_dict = {}
    for kaiko_path in kaiko_paths:
        kaiko_dict[os.path.splitext(os.path.basename(kaiko_path))[0][:-4]] = kaiko_path
    return kaiko_dict
    # kaiko_name = None
    # for kaiko_path in kaiko_files:
    #     kaiko_name = os.path.basename(kaiko_path)
    # return os.path.splitext(kaiko_name)[0][:-4]
    

def get_casanovo_path_by_mgf_filename(casanovo_paths):
    """map a mgf filename into a corresponding casanovo output path
    Args:
        casanovo_paths: list
            a list of casanovo file paths
    Return:
        mgf_dict: dict
            {
                'input1': 'path/to/casanovo_20230710143300.mztab',
                'input2': 'path/to/casanovo_20230710143400.mztab',
                ...
            }
    """
    csnv_dict = {}
    for casanovo_path in casanovo_paths:
        tables = mztab.MzTab(casanovo_path)
        mgf_name_path = tables.metadata['ms_run[1]-location'][8:]
        mgf_name = os.path.splitext(os.path.basename(mgf_name_path))[0]
        csnv_dict[mgf_name] = casanovo_path
    return csnv_dict
    # casanovo_name = None
    # for casanovo_path in casanovo_files:
    #     os.path.splitext(os.path.basename(mztab.MzTab(casanovo_path).metadata['ms_run[1]-location'][8:]))[0]
    #     tables = mztab.MzTab(casanovo_path)
    #     casanovo_name_path = tables.metadata['ms_run[1]-location'][8:]
    #     casanovo_name = os.path.basename(casanovo_name_path)
    # return os.path.splitext(casanovo_name)[0]
    # mgf_dict['casanovo_names'] = [os.path.splitext(os.path.basename(mztab.MzTab(casanovo_path).metadata['ms_run[1]-location'][8:]))[0] for casanovo_path in casanovo_files]


if __name__ == '__main__':
    # if get_mgf_names(mgf_files) == get_kaiko_names(kaiko_files) == get_casanovo_names(casanovo_files):
    #     print('yes')
    mgf_path_by_mgf_name = get_mgf_path_by_filename(mgf_files)
    # print(mgf_path_by_mgf_name)
    kaiko_path_by_mgf_name = get_kaiko_path_by_mgf_filename(kaiko_files)
    # print(kaiko_path_by_mgf_name)
    casanovo_path_by_mgf_name = get_casanovo_path_by_mgf_filename(casanovo_files)
    # print(casanovo_path_by_mgf_name)

    # mgf 1,2,3
    # kaiko 2,3,4
    # casa 3,4,5
    # (3, 3, 3)

    common_paths = []
    for mgf_name, mgf_path in mgf_path_by_mgf_name.items():
        if mgf_name in kaiko_path_by_mgf_name and mgf_name in casanovo_path_by_mgf_name:
            common_paths.append((
                mgf_path,
                kaiko_path_by_mgf_name[mgf_name],
                casanovo_path_by_mgf_name[mgf_name]
            ))
    # run
    dfs = []
    for (mgf_path, kaiko_path, casanovo_path) in common_paths:
        temp_df = aggregate_kaiko_casanovo(mgf_path, kaiko_path, casanovo_path, output_path)
        temp_df['mgf'] = mgf_path
        temp_df['kaiko'] = kaiko_path
        temp_df['casanovo'] = casanovo_path
        dfs.append(temp_df)
    compare_df = pd.concat(dfs, ignore_index=True)
    print(compare_df)

    compare_df.to_csv(output_path, index=None, sep='\t')

