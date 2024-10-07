import os
import re
import sys
import pandas as pd

from math import floor
from pathlib import Path, PureWindowsPath

def prepare_denovo_command(mgf_file, denovout_dir, config):
    ## Step 1. Run Denovo using subprocess.
    if config['denovo']['method'] == "Deepnovo":
        denovo_args = [sys.executable, "src/kaiko_main.py", 
                        "--mgf_dir", mgf_file.resolve(), 
                        "--train_dir", "model/",
                        "--decode_dir", denovout_dir.resolve(),
                        "--profile", config['denovo']['profile']]

        if config['denovo']['topk']:
            denovo_args = denovo_args + ["--topk"]
        if config['denovo']['multi_decode']:
            denovo_args = denovo_args + ["--multi_decode"]
        if config['denovo']['beam_search']:
            denovo_args = denovo_args + ["--beam_search", "--beam_size", config['denovo']['beam_size']]     
        cwd_folder = "Kaiko_Deepnovo"    
    
    elif config['denovo']['method'] == "PNNL_Casanovo":
        model_path = "PNNL_casanovo.ckpt"
        config_path  = "PNNL_casanovo_config.yaml"
        mztab_path = denovout_dir / f'{str(mgf_file.stem)}_denovo.mztab'
        denovo_args = ["casanovo", "sequence", 
                        "--model", model_path, 
                        "--config", config_path,
                        "--output", mztab_path.resolve(),
                        mgf_file.resolve()]
        cwd_folder = "Casanovo"
    
    elif config['denovo']['method'] == "Casanovo_massivekb":
        model_path = "casanovo_massivekb.ckpt"
        config_path  = "casanovo_massivekb_config.yaml"
        mztab_path = denovout_dir / f'{str(mgf_file.stem)}_denovo.mztab'
        denovo_args = ["casanovo", "sequence", 
                        "--model", model_path, 
                        "--config", config_path,
                        "--output", mztab_path.resolve(),
                        mgf_file.resolve()]
        cwd_folder = "Casanovo"
    for i in range(len(denovo_args)):
            denovo_args[i] = str(denovo_args[i])
    return (denovo_args, cwd_folder)

# @profile
def combine_denovo_output(directory, prefix, denovo_method, selection = 0.25):

    if denovo_method == "Deepnovo":
        files = [f for f in os.listdir(directory) if bool(re.search(r'_out.txt', f))]
    else:
        files = [f for f in os.listdir(directory) if bool(re.search(r'_denovo.mztab', f))]
    samples = []

    if denovo_method == "Deepnovo":
        for file in files:
            xx = pd.read_csv(directory / file, sep = "\t", header = 0)

            xx['output_seq'] = [re.sub(",", "", str(peptide)) for peptide in xx['output_seq']]
            xx['output_seq'] = [re.sub("mod", "", str(peptide)) for peptide in xx['output_seq']]
            xx = xx.loc[xx['output_score'] != float("Inf")]

            xx = xx.sort_values('output_score', ascending = False)
            xx = xx.head(floor(selection * floor(len(xx.index))))
            #xx = xx[['scan', 'output_seq']]

            xx['pep_length'] = [len(peptide) for peptide in xx['output_seq']]
            xx = xx.loc[(xx['pep_length'] >= 10) & (xx['pep_length'] <= 17)]
            xx['rank'] = list(range(1, len(xx.index) + 1))

            grouped = xx.groupby('output_seq')

            summary = grouped.apply(summary_times).to_frame()
            # summary = grouped.apply(summary_times)
            summary['output_seq'] = summary.index
            summary.columns = ['times', 'output_seq']
            summary = summary[['output_seq', 'times']]
            summary['rank'] = grouped.apply(summary_rank)
            summary['scans'] = grouped.apply(summary_scans)

            samples += [summary]
    else:
        for file in files:
            filename = directory / file
            with filename.open() as file:
                skip_lines = 0
                line = file.readline()
                while bool(re.search(r'^MTD	', line)):
                    skip_lines = skip_lines + 1
                    line = file.readline()
            xx = pd.read_csv(filename, sep = "\t", header = 0, skiprows = skip_lines)
            xx['mass_to_charge_err'] = [abs(xx['calc_mass_to_charge'][i] - xx['exp_mass_to_charge'][i]) for i in range(len(xx))]
            xx['pep_length'] = [len(re.split(r'(?<=.)(?=[A-Z])', peptide)) for peptide in xx['sequence']]
            xx['scan'] = [f'{filename.stem.replace("_denovo", "")}_ID={str(ID)}' for ID in xx['PSM_ID']]
            xx = xx.loc[(xx['pep_length'] >= 10) & (xx['pep_length'] <= 30)]

            xx = xx.sort_values('search_engine_score[1]', ascending = False)
            xx = xx.head(floor(selection * floor(len(xx.index))))
            
            xx['rank'] = list(range(1, len(xx.index) + 1))
            grouped = xx.groupby('sequence')

            summary = grouped.apply(summary_times).to_frame()
            # summary = grouped.apply(summary_times)
            summary['output_seq'] = summary.index
            summary['output_seq'] = [peptide.replace('+15.995', '') for peptide in summary['output_seq']]
            summary.columns = ['times', 'output_seq']
            summary = summary[['output_seq', 'times']]
            summary['rank'] = grouped.apply(summary_rank)
            summary['scans'] = grouped.apply(summary_scans)

            samples += [summary]

    combined_fasta = directory / (prefix + '_combined_denovo.fasta')
    # Path('Kaiko_volume/Kaiko_intermediate/' + prefix + '_combined_denovo.fasta')
    with combined_fasta.open('w') as fasta_file:
        for index in range(len(samples)):
            summary = samples[index]
            nms = [">S" + str(index) + '_' + summary['scans'][i] + "_" + str(summary['times'][i]) for i in range(len(summary))]
            to_write = [None]*2*len(nms)
            to_write[::2] = nms
            to_write[1::2] = summary['output_seq']       
            for line in to_write:
                fasta_file.write(f"{line}\n")
                # fasta_file.write(line + "\n")


def summary_times(group):
    scans = group['scan']
    return len(scans)

def summary_scans(group):
    scans = group['scan']
    return "_AND_".join(scans)

def summary_rank(group):
    ranks = group['rank']
    return min(ranks)

