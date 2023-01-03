import subprocess
import configparser
import os

from Kaiko_4 import aggregate_fasta
from Kaiko_3 import run_diamond_tally
from Kaiko_2 import combine_denovo_output


## Parsing

config = configparser.ConfigParser()
config.sections()

## Read the defaults
config.read('kaiko_defaults.ini')

## Read any changes to defaults by user
if (os.path.isfile('config.ini')):
    config.read('config.ini')



## Step by step pipeline.


## Step 1. Run Denovo using subprocess.
kaiko_1_args = ["python", "src/kaiko_main.py", 
                "--mgf_dir", "../" + config['denovo']['mgf_dir'], 
                "--train_dir", config['denovo']['train_dir'],
                "--decode_dir", "../" + config['denovo']['decode_dir']]

if config['denovo'].getboolean('topk'):
    kaiko_1_args = kaiko_1_args + ["--topk"]
if config['denovo'].getboolean('multi_decode'):
    kaiko_1_args = kaiko_1_args + ["--multi_decode"]
if config['denovo'].getboolean('beam_search'):
    kaiko_1_args = kaiko_1_args + ["--beam_search", "--beam_size", config['denovo']['beam_size']]     

print("DeNovo: Running the following command:\n")
print(" ".join(kaiko_1_args) + "\n")

# subprocess.run(kaiko_1_args, cwd = "Kaiko_denovo")
subprocess.call(kaiko_1_args, cwd = "Kaiko_denovo")



## Step 2. Combine into fasta
print("\n Combinining denovo output\n")
combine_denovo_output(config['denovo']['decode_dir'])




## Step 3. Passing to diamond
diamond_args = ["./diamond", "blastp", "-d",
                "../uniref100", "--min-score", "1",
                "-q", "../../Kaiko_intermediate/combined_denovo.fasta", "-o",
                "diamond_search_output.dmd", "-f", "6", "qseqid", 
                "stitle", "pident", "evalue", "mismatch"]


print("DeNovo: Running the following command:\n")
print(" ".join(diamond_args) + "\n")

# subprocess.run(diamond_args, cwd = "diamond-2.0.15")
os.chdir("Kaiko_volume/Kaiko_stationary_files/diamond-linux")
os.system(" ".join(diamond_args))
os.chdir("../../")




# Step 4. Tallying the diamond results
run_diamond_tally("Kaiko_volume/Kaiko_intermediate/diamond_search_output.dmd", 
                  int(config['diamond tally']['ntops']), 
                  config['diamond tally']['ncbi_taxa_folder'] , 
                  config['diamond tally']['mode'] , 
                  config['diamond tally']['fout'] , 
                  float(config['diamond tally']['pident']))




## Step 5. Putting together the final fasta file.


if config['taxa to fasta']['kingdom_list'] != "":
    kingdom_list = config['taxa to fasta']['kingdom_list'].split(', ')
    print(kingdom_list)
else:
    kingdom_list = []


aggregate_fasta(config['taxa to fasta']['ref_fasta'],
                config['taxa to fasta']['diamond_tally'],
                config['taxa to fasta']['fout'],
                int(config['taxa to fasta']['ntops']),
                config['taxa to fasta']['taxa_key'],
                kingdom_list)




