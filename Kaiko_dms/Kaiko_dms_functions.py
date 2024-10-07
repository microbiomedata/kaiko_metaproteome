# import pyodbc 
import pandas as pd
import time
import gzip
import glob
import os
import sys
import subprocess

from pathlib import PureWindowsPath, Path
from pyteomics import mzml, auxiliary

def get_request_dataset_paths(job_req_id):
    # cnxn = pyodbc.connect("DRIVER={SQL Server};SERVER=gigasax;DATABASE=dms5;")
    cnxn = None
    # cnxn = pyodbc.connect(f"DRIVER={pyodbc.drivers()[0]};SERVER=gigasax;DATABASE=dms5;")
    sql_str_job_req = f"SELECT * FROM v_analysis_job_request_detail_report WHERE [request]={job_req_id}"
    datasets_df = pd.read_sql(sql_str_job_req, cnxn)

    datasets = datasets_df.datasets.to_list()[0].split(', ')
    # sql_str_dataset = f"SELECT * FROM v_dataset_detail_report_ex WHERE [dataset] = '{datasets[0]}'"
    # results_path = pd.read_sql(sql_str_dataset, cnxn).dataset_folder_path.to_list()[0].split(': ')[-1]
    # results_folder = Path(PureWindowsPath(results_path)).parent

    mzml_paths = []
    for dataset in datasets:
        sql_str_dataset = f"SELECT * FROM v_dataset_detail_report_ex WHERE [dataset] = '{dataset}'"
        results_path = pd.read_sql(sql_str_dataset, cnxn).dataset_folder_path.to_list()[0].split(': ')[-1]
        # results_folder = Path(PureWindowsPath(results_path)).parent
        # results_path = results_folder / dataset 
        results_path = Path(PureWindowsPath(results_path))
        assert results_path.exists()

        mzml_results_pointer = [x for x in results_path.glob('MSXML*')]
        assert len(mzml_results_pointer) == 1

        mzml_results_pointer = [x for x in mzml_results_pointer[0].glob(f'{dataset}*')]
        assert len(mzml_results_pointer) == 1

        with mzml_results_pointer[0].open('r') as mzml_pointer:
            mzml_path = mzml_pointer.readline()
            mzml_path = mzml_path.split('\n')[0]

        mzml_path = Path(PureWindowsPath(mzml_path))
        assert mzml_path.exists()

        mzml_paths = mzml_paths + [mzml_path]
        ## To not spam DMS if there are too many datasets
        time.sleep(1)

    return(mzml_paths)


def get_dataset_mzml_mzid(analysis_job_ID, species_log, local_mzml_folder):
    cnxn = pyodbc.connect("DRIVER={SQL Server};SERVER=gigasax;DATABASE=dms5;")
    # sql_str_dataset = f"SELECT * FROM v_dataset_detail_report_ex WHERE [ID] = '{dataset_ID}'"
    sql_str = f"SELECT * FROM v_analysis_job_detail_report_2 WHERE [job] = {analysis_job_ID}"
    sql_res = pd.read_sql(sql_str, cnxn)
    dataset_path = sql_res.dataset_folder_path[0].split(': ')[-1]
    results_path = sql_res.results_folder_path[0]
    dataset = sql_res.dataset[0]
    # results_folder = Path(PureWindowsPath(results_path)).parent
    # results_path = results_folder / dataset 
    results_path = Path(PureWindowsPath(results_path))
    dataset_path = Path(PureWindowsPath(dataset_path))
    assert results_path.exists()

    # mzid_path = results_path / f'{dataset}_msgfplus.mzid'
    mzid_path = [x for x in results_path.glob(f'{dataset}_msgfplus.mzid*')]
    mzid_parameter_path = [x for x in results_path.glob('JobParameters_*')][0]
    assert len(mzid_path) == 1
    mzid_path = mzid_path[0]

    convert_from_raw = False
    mzml_results_pointer = [x for x in dataset_path.glob('MSXML*')]
    if len(mzml_results_pointer) == 1:
        mzml_parameter_path = [x for x in mzml_results_pointer[0].glob('JobParameters_*')][0]
        mzml_parameter_path = Path(PureWindowsPath(mzml_parameter_path))
        mzml_results_pointer = [x for x in mzml_results_pointer[0].glob(f'{dataset}*')]
        assert len(mzml_results_pointer) == 1

        with mzml_results_pointer[0].open('r') as mzml_pointer:
            mzml_path = mzml_pointer.readline()
            mzml_path = mzml_path.split('\n')[0]
        mzml_path = Path(PureWindowsPath(mzml_path))
        if not mzml_path.exists():
            convert_from_raw = True
    else:
        convert_from_raw = True
            
    if convert_from_raw:
        print("No mzML cache found. Creating local mzML from raw file")
        mzml_parameter_path = None
        dataset_raw_path = [x for x in dataset_path.glob(f'{dataset}*.raw')][0]
        start_time = time.time()
        species_log.write(f'\nConverting raw file for dataset {dataset} to mzML, since the cached mzML does not exist.')
        species_log.write(f'\nThe command to MSConvert is: ')
        msconvert_commands = ['MSConvert.exe', str(dataset_raw_path), '--filter', '"peakPicking vendor msLevel=1-"', "--mzML", "--32", "-g", "-o", str(local_mzml_folder), "--outfile", f'{dataset}']
        species_log.write(" ".join(msconvert_commands))
        subprocess.run(msconvert_commands)
        species_log.write(f'\nTook a total of {time.time() - start_time} to convert to mzML.\n')
        species_log.flush()
        mzml_path = local_mzml_folder / f'{dataset}.mzML.gz'

    
    assert mzml_path.exists()

    ## To not spam DMS if there are too many datasets
    time.sleep(1)
    return mzid_path, mzml_path, mzid_parameter_path, mzml_parameter_path


# def convert_raw_to_mzml()


def inspect_mzML_file(fpath, gzipped=True):
    spectra = []
    if gzipped:
        f = gzip.open(fpath, 'rb')
    else:
        f = fpath
    for obj in mzml.read(f):
        spectra.append(obj)
    if gzipped: f.close()
    return spectra

def generate_mgf_without_annotation(mzml_spectra, file_index=0, ntops=500, out_file='out.mgf'):
    num_spectra = 0
    with open(out_file, 'w') as f:
        for spectrum in mzml_spectra:
            if spectrum['ms level'] != 2:
                continue
            scan = int(spectrum['id'].split('scan=')[1])
            try:
            
                mz_arr = spectrum['m/z array']
                int_arr = spectrum['intensity array']
                rtsec = 60.0*(spectrum['scanList']['scan'][0]['scan start time'])
                selectedIon = spectrum['precursorList']['precursor'][0]['selectedIonList']['selectedIon'][0]

                assert len(mz_arr) == len(int_arr), "[ERR] Wrong data format: len(mz_arr) != len(int_arr)"

                print("BEGIN IONS", file=f)
                print("TITLE={0}.{1}".format(file_index, scan), file=f)
                print("PEPMASS={0}".format(selectedIon['selected ion m/z']), file=f)
                # sometimes they don't have a charge info. if so, we use the annotation file
                if 'charge state' in selectedIon:
                    print("CHARGE={0:d}+".format(int(selectedIon['charge state'])), file=f)
                else:
                    print("CHARGE={0:d}+".format(999), file=f)
                print("SCANS={0}:{1}".format(file_index, scan), file=f)
                print("RTINSECONDS={0}".format(rtsec), file=f)
                print("SEQ=UNKNOWN", file=f)
                for i in range(len(mz_arr)):
                    print("{0} {1}".format(mz_arr[i], int_arr[i]), file=f)
                print("END IONS", file=f)
                num_spectra += 1
            except:
                print('[ERR]', scan, spectrum)
                continue
                    
        return num_spectra

def generate_mgf_files(data_dir, dest_dir = './', dataset_pattern = '', gzipped = True):
    # collect mzML.gz files
    if dataset_pattern != '':
        if gzipped:
            mzML_files = glob.glob(data_dir + f"/{dataset_pattern}")
        else:
            mzML_files = glob.glob(data_dir + f"/{dataset_pattern}")
    else:
        if gzipped:
            mzML_files = glob.glob(data_dir + "/*.mzML.gz")
        else:
            mzML_files = glob.glob(data_dir + "/*.mzML")
        
    if not os.path.exists(dest_dir):
        os.makedirs(dest_dir)

    mzML_log_handler = open(dest_dir + '/mgf_list.log', 'w+')
    print("id\tmgf_file\tnum_scans\ttotal_scans", file=mzML_log_handler)
    
    start_time = time.time()
    num_mzML_files = len(mzML_files)
    total_scans = 0
    for i, mzML_file in enumerate(mzML_files):
        try:
            if gzipped:
                common_name = os.path.basename(mzML_file).rsplit('.mzML.gz')[0]
            else:
                common_name = os.path.basename(mzML_file).rsplit('.mzML')[0]
        
            if os.path.exists(dest_dir + '/' + common_name + '.mgf'):
                print('[{0:3d}/{1:3d}] {2}, Already exists' \
                      .format(i+1,
                              num_mzML_files,
                              common_name))
                continue

            num_spectra = 0
            mzml_spectra = inspect_mzML_file(mzML_file, gzipped)
            num_spectra = generate_mgf_without_annotation(mzml_spectra,
                                                            file_index=f'{dataset_pattern}--{i}',
                                                            out_file=dest_dir + '/' + common_name + '.mgf')
            num_scans = num_spectra
            total_scans += num_spectra
            msg = "SUCCESS"
            print('[{0:3d}/{1:3d}] {2}, {3:d}/{4:d}/{5:d}, {6:.2f}sec' \
                      .format(i+1,
                              num_mzML_files,
                              common_name,
                              num_spectra,
                              num_scans,
                              total_scans,
                              time.time()-start_time))
            print("{0}\t{1}\t{2}\t{3}".format(i, common_name, num_spectra, total_scans), file=mzML_log_handler)
            sys.stdout.flush()
        except Exception as e:
            print('[ERR] {}'.format(mzML_file))
            print('[ERR]', e)



