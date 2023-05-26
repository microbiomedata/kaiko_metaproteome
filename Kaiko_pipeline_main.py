import subprocess
import cProfile
import yaml
import os
import re
import argparse

from pathlib import Path, PureWindowsPath
from s3path import S3Path

from Kaiko_4 import aggregate_fasta
from Kaiko_3 import run_diamond_tally
from Kaiko_2 import combine_denovo_output

parser = argparse.ArgumentParser()

parser.add_argument("-config", "--config", help="YAML file with parameters for pipeline, and input files path. Any found values in config replace the default values.")
args = parser.parse_args()


## Parsing
kaiko_defaults_path = Path('kaiko_defaults.yaml')
config = yaml.safe_load(kaiko_defaults_path.open())

assert args.config is not None, "Please provide a config file with the spectrum path (mgf files). Use flag --config to pass the path. See the tutorial for more information."
user_config_path = Path(args.config)
assert user_config_path.exists(), "File " + str(user_config_path.absolute()) + " does not exist."
config_user = yaml.safe_load(user_config_path.open())

## Overriding defaults if values found in user config.
for section in config_user.keys():
    for param in config_user[section].keys():
        config[section][param] = config_user[section][param]


## handling any backslashes with this. All final paths used here have forward slashes, 
## as they are compatible in Windows, Linux, and Mac.
working_dir = Path(Path('.').resolve().as_posix())
mgf_dir = Path(PureWindowsPath(config['denovo']['mgf_dir']).as_posix())
ncbi_taxa_folder = Path(PureWindowsPath(config['diamond tally']['ncbi_taxa_folder']).as_posix())
ref_fasta = Path(PureWindowsPath(config['taxa to fasta']['ref_fasta']).as_posix())
diamond_folder = Path(PureWindowsPath(config['diamond tally']['diamond_folder']).as_posix())
prefix = mgf_dir.name

denovout_dir = Path('Kaiko_volume/Kaiko_intermediate/denovo_output/' + prefix)
if not denovout_dir.exists():
    denovout_dir.mkdir()

## Step 1. Run Denovo using subprocess.
kaiko_1_args = ["python", "src/kaiko_main.py", 
                "--mgf_dir", mgf_dir.resolve(), 
                "--train_dir", "model/",
                "--decode_dir", denovout_dir.resolve(),
                "--profile", config['denovo']['profile']]

if config['denovo']['topk']:
    kaiko_1_args = kaiko_1_args + ["--topk"]
if config['denovo']['multi_decode']:
    kaiko_1_args = kaiko_1_args + ["--multi_decode"]
if config['denovo']['beam_search']:
    kaiko_1_args = kaiko_1_args + ["--beam_search", "--beam_size", config['denovo']['beam_size']]     

print("DeNovo: Running the following command:\n")
for i in range(len(kaiko_1_args)):
    kaiko_1_args[i] = str(kaiko_1_args[i])

if not config['denovo']['cached']:
    print(" ".join(kaiko_1_args) + "\n")
    subprocess.run(kaiko_1_args, cwd = "Kaiko_denovo")
    # subprocess.call(kaiko_1_args, cwd = "Kaiko_denovo")

if (config['denovo'])['profile']:
    profiler = cProfile.Profile()
    profiler.enable()

## Step 2. Combine into fasta
print("\n Combinining denovo output\n")

combine_denovo_output(denovout_dir, prefix)

denovo_combined_fasta = Path("Kaiko_volume/Kaiko_intermediate/" + prefix + "_combined_denovo.fasta")
diamond_search_out = Path("Kaiko_volume/Kaiko_intermediate/" + prefix + "_diamond_search_output.dmd")

## Step 3. Passing to diamond
if os.name == 'posix':
    diamond_args = ["./diamond"]
else:
    diamond_args = ["diamond"]
    
diamond_args = diamond_args + ["blastp", "-d",
                "../uniref100", "--min-score", "1",
                "-q", denovo_combined_fasta.resolve().as_posix(), "-o",
                diamond_search_out.resolve().as_posix(), "-f", "6", "qseqid", 
                "stitle", "pident", "evalue", "mismatch"]

if not config['diamond tally']['cached']:
    print("DeNovo: Running the following command:\n")
    print(" ".join(diamond_args) + "\n")    
    # os.chdir("Kaiko_volume/Kaiko_stationary_files/diamond-2.0.15")
    os.chdir(diamond_folder)
    print(os.getcwd())
    os.system(" ".join(diamond_args))
    os.chdir(working_dir)

nprot = '{:.5e}'.format(int(config['diamond tally']['n_protein_cutoff']))
top_strains = str(config['taxa to fasta']['top_strains'])
kaiko_tally = Path("Kaiko_volume/Kaiko_intermediate/" + prefix + "_kaiko_prediction" + f'_top_taxa_nprot_{nprot}_top_{top_strains}_strains.csv')

# Step 4. Tallying the diamond results
run_diamond_tally(diamond_search_out, 
                  int(config['taxa to fasta']['top_strains']), 
                  ncbi_taxa_folder, 
                  config['diamond tally']['mode'], 
                  kaiko_tally, 
                  float(config['diamond tally']['pident']),
                  int(config['diamond tally']['n_protein_cutoff']))


## Step 5. Putting together the final fasta file.

if config['taxa to fasta']['kingdom_list'] != "":
    kingdom_list = config['taxa to fasta']['kingdom_list'].split(', ')
    print(kingdom_list)
else:
    kingdom_list = []

coverage_target = str(config['taxa to fasta']['coverage_target'])
kaiko_final_output = Path("Kaiko_volume/Kaiko_output/" + prefix + f'_kaiko_fasta_coverage_{coverage_target}_nprot_{nprot}_top{top_strains}_strains.fasta')

aggregate_fasta(ref_fasta,
                kaiko_tally,
                kaiko_final_output,
                config['taxa to fasta']['coverage_target'],
                int(config['taxa to fasta']['top_strains']),
                ncbi_taxa_folder,
                config['taxa to fasta']['taxa_key'],
                kingdom_list)


if (config['denovo'])['profile']:
    profiler.enable()
    profiler.dump_stats('Kaiko_volume/Kaiko_taxa.prof')
