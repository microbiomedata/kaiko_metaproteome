import subprocess
import configparser
import cProfile
import yaml
import os
import re

from Kaiko_4 import aggregate_fasta
from Kaiko_3 import run_diamond_tally
from Kaiko_2 import combine_denovo_output


## Parsing

config = yaml.safe_load(open('kaiko_defaults.yaml'))

if (os.path.isfile('config.yaml')):
    config_user = yaml.safe_load(open('config.yaml'))

    for section in config_user.keys():
        for param in config_user[section].keys():
            config[section][param] = config_user[section][param]

# config = configparser.ConfigParser()
# config.sections()

# ## Read the defaults
# config.read('kaiko_defaults.ini')

# ## Read any changes to defaults by user

#     config.read('config.ini')




prefix = config['denovo']['mgf_dir']
prefix = re.sub("^.*Kaiko_volume/Kaiko_input_files/", "", prefix)
if not os.path.exists('Kaiko_volume/Kaiko_intermediate/denovo_output/' + prefix):
    os.makedirs('Kaiko_volume/Kaiko_intermediate/denovo_output/' + prefix)

## Step 1. Run Denovo using subprocess.
kaiko_1_args = ["python", "src/kaiko_main.py", 
                "--mgf_dir", "../" + config['denovo']['mgf_dir'], 
                "--train_dir", config['denovo']['train_dir'],
                "--decode_dir", "../Kaiko_volume/Kaiko_intermediate/denovo_output/" + prefix,
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
print(" ".join(kaiko_1_args) + "\n")

subprocess.run(kaiko_1_args, cwd = "Kaiko_denovo")
# subprocess.call(kaiko_1_args, cwd = "Kaiko_denovo")

if (config['denovo'])['profile']:
    profiler = cProfile.Profile()
    profiler.enable()

## Step 2. Combine into fasta
print("\n Combinining denovo output\n")


combine_denovo_output("Kaiko_volume/Kaiko_intermediate/denovo_output/" + prefix + '/', prefix)


## Step 3. Passing to diamond
diamond_args = ["./diamond", "blastp", "-d",
                "../uniref100", "--min-score", "1",
                "-q", "../../Kaiko_intermediate/" + prefix + "_combined_denovo.fasta", "-o",
                "../../Kaiko_intermediate/" + prefix + "_diamond_search_output.dmd", "-f", "6", "qseqid", 
                "stitle", "pident", "evalue", "mismatch"]


print("DeNovo: Running the following command:\n")
print(" ".join(diamond_args) + "\n")

# subprocess.run(diamond_args, cwd = "diamond-2.0.15")
os.chdir("Kaiko_volume/Kaiko_stationary_files/diamond-linux")
os.system(" ".join(diamond_args))
os.chdir("../../../")

# Step 4. Tallying the diamond results
run_diamond_tally("Kaiko_volume/Kaiko_intermediate/" + prefix + "_diamond_search_output.dmd", 
                  int(config['diamond tally']['ntops']), 
                  config['diamond tally']['ncbi_taxa_folder'] , 
                  config['diamond tally']['mode'] , 
                  "Kaiko_volume/Kaiko_intermediate/" + prefix + "_kaiko_prediction_top_taxa.csv", 
                  float(config['diamond tally']['pident']))


## Step 5. Putting together the final fasta file.


if config['taxa to fasta']['kingdom_list'] != "":
    kingdom_list = config['taxa to fasta']['kingdom_list'].split(', ')
    print(kingdom_list)
else:
    kingdom_list = []


aggregate_fasta(config['taxa to fasta']['ref_fasta'],
                "Kaiko_volume/Kaiko_intermediate/" + prefix + "_kaiko_prediction_top_taxa.csv",
                "Kaiko_volume/Kaiko_output/" + prefix + "_kaiko_output.fasta",
                int(config['taxa to fasta']['ntops']),
                config['taxa to fasta']['taxa_key'],
                kingdom_list)


if (config['denovo'])['profile']:
    profiler.enable()
    profiler.dump_stats('Kaiko_volume/Kaiko_taxa.prof')
