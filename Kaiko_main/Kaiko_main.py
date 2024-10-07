import subprocess
import cProfile
import yaml
import os
import re
import argparse


from pathlib import Path, PureWindowsPath
from s3path import S3Path

from .Kaiko_4 import aggregate_fasta
from .Kaiko_3 import run_diamond_tally
from .Kaiko_2 import prepare_denovo_command, combine_denovo_output

from Kaiko_dms.Kaiko_dms_functions import get_request_dataset_paths, generate_mgf_files

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
output_dir = Path(PureWindowsPath(config['general']['output_dir']).as_posix())
ncbi_taxa_folder = Path(PureWindowsPath(config['diamond tally']['ncbi_taxa_folder']).as_posix())
ref_fasta = Path(PureWindowsPath(config['taxa to fasta']['ref_fasta']).as_posix())
diamond_folder = Path(PureWindowsPath(config['diamond tally']['diamond_folder']).as_posix())
diamond_database = Path(PureWindowsPath(config['diamond tally']['diamond_database']).as_posix())
ref_fasta_igzip_index = Path(PureWindowsPath(config['taxa to fasta']['gz_index']).as_posix())
index_path = Path(PureWindowsPath(config['taxa to fasta']['proteome_index']).as_posix())
index_s_path = Path(PureWindowsPath(config['taxa to fasta']['proteome_index_s']).as_posix())

ref_proteome_log = Path(PureWindowsPath(config['taxa to fasta']['ref_proteome_log']).as_posix())
prefix = mgf_dir.name
if config['diamond tally']['db_pattern'] == 'TaxID':
    mode = 'uniref100'
elif config['diamond tally']['db_pattern'] == 'OX':
    mode = 'ref_prot'

## Creating drectories in output folder:
denovout_dir = output_dir / ('Kaiko_intermediate/' + prefix + '/denovo_output/')
intermediate_dir = output_dir / ("Kaiko_intermediate/" + prefix)
final_dir = output_dir / ('Kaiko_output/' + prefix)
denovout_dir.mkdir(parents=True, exist_ok=True)
denovo_completion_log = denovout_dir / 'denovo_completion_log.txt'
intermediate_dir.mkdir(parents=True, exist_ok=True)
final_dir.mkdir(parents=True, exist_ok=True)

if not config['denovo']['cached']:
    all_mgf = []
    
    if config['denovo']['dms_analysis_job'] != '':
        # print('dms access from docker not implemented')
        print('Gathering paths from dms')
        mzml_paths = get_request_dataset_paths(config['denovo']['dms_analysis_job'])
        for mzml_path in mzml_paths:
            individual_mgf_name = str(mzml_path.name).split('.mzML.gz')[0]
            all_mgf = all_mgf + [f'{individual_mgf_name}.mgf']
            if config['denovo']['method'] == "Deepnovo":  
                expected_output_path = denovout_dir / f'{individual_mgf_name}_out.txt'
            else:
                expected_output_path = denovout_dir / f'{individual_mgf_name}_denovo.mztab'
            individual_mgf_dir = mgf_dir / individual_mgf_name
            
            ## If denovo is done, don't need to do it again
            if not expected_output_path.exists():
                ## If another process already started working on this dataset, skip.
                if not individual_mgf_dir.exists():
                    print(f'Converting file {str(mzml_path.name)}')
                    generate_mgf_files(str(mzml_path.parent), dest_dir = str(individual_mgf_dir.resolve()), dataset_pattern = str(mzml_path.name), gzipped = True)

                    kaiko_1_args, cwd_folder = prepare_denovo_command(individual_mgf_dir, denovout_dir, config)
                    print(" ".join(kaiko_1_args) + "\n")
                    subprocess.run(kaiko_1_args, cwd = cwd_folder)

                    with denovo_completion_log.open('a') as completion_log:
                        completion_log.write(f'{individual_mgf_name}.mgf\t{expected_output_path.name}\t completed denovo sequencing\n')

                    if not config['denovo']['keep_dms_locally']:
                        print(f'Removing the locally stored coverted spectra {str(individual_mgf_dir)}.mgf')
                        try:
                            os.remove(individual_mgf_dir / f'{str(individual_mgf_name)}.mgf')
                        except:
                            print('woops')
                else:
                    print(f'{str(individual_mgf_dir)} already exists, another process is working on this dataset. Moving onto the next dataset')
            else:
                print(f'{str(expected_output_path.name)} already exists. Moving onto the next dataset')
    else:
        mgf_files = [x for x in mgf_dir.glob('*.mgf')]
        for mgf_file in mgf_files:
            mgf_name = str(mgf_file.name).split('.mgf')[0]
            all_mgf = all_mgf + [str(mgf_file.name)]
            if config['denovo']['method'] == "Deepnovo":  
                expected_output_path = denovout_dir / f'{mgf_name}_out.txt'
            else:
                expected_output_path = denovout_dir / f'{mgf_name}_denovo.mztab'
            ## If another process is working on this dataset, or if the denovo for that dataset is already done, skip
            if not expected_output_path.exists():
                ## Make the expected output handle immediately, so any other process knows NOT to start the same dataset
                with expected_output_path.open('w') as output_file:
                    pass
                kaiko_1_args, cwd_folder = prepare_denovo_command(mgf_file, denovout_dir, config)
                print("DeNovo: Running the following command:\n")
                print(" ".join(kaiko_1_args) + "\n")
                subprocess.run(kaiko_1_args, cwd = cwd_folder)
                
                with denovo_completion_log.open('a') as completion_log:
                    completion_log.write(f'{mgf_file.name}\t{expected_output_path.name}\t completed denovo sequencing\n')
            else:
                print(f'{str(expected_output_path)} already exists, another process is working on this or already finished this dataset. Moving onto the next dataset')

    with denovo_completion_log.open('r') as completion_log:
        completed_mgf = []
        for line in completion_log:
            completed_mgf = completed_mgf + [line.split('\t')[0]]
        ## We ONLY proceed whenever ALL datasets are DONE. If not, then the process fails and ONLY the process that finishes last will continue 
        ## with the DIAMOND search
        assert set(completed_mgf) == set(all_mgf)
        print('All the denovo sequencing has been finished. Moving to DIAMOND search')
        f = open(denovout_dir / f'{prefix}_config.yaml', 'a')
        yaml.safe_dump(config, f)   

if (config['denovo']['profile']):
    profiler = cProfile.Profile()
    profiler.enable()

nprot = '{:.5e}'.format(int(config['diamond tally']['n_protein_cutoff']))
top_strains = str(config['taxa to fasta']['top_strains'])
benchmark = config['diamond tally']['benchmark']
if config['diamond tally']['db_pattern'] == 'OX':
    db_name = 'ref_prot'
elif config['diamond tally']['db_pattern'] == 'TaxID':
    db_name = 'uniref100'
if benchmark:
    suffix = f'_benchmark_{benchmark}'
else:
    suffix = ''

denovo_combined_fasta = denovout_dir / (f'{prefix}_combined_denovo.fasta')
diamond_search_out = intermediate_dir / (f'{prefix}_diamond_search_output_{db_name}.dmd')

if not config['diamond tally']['cached']:
    ## Step 2. Combine into fasta
    print("\n Combinining denovo output\n")

    combine_denovo_output(denovout_dir, prefix, config['denovo']['method'])
    
    ## Step 3. Passing to diamond
    if os.name == 'posix':
        diamond_args = [f'{diamond_folder}/diamond']
    else:
        diamond_args = [f'{diamond_folder}\diamond']
        
    diamond_args = diamond_args + ["blastp", "-d",
                    diamond_database.resolve().as_posix(), "--min-score", "1",
                    "-q", denovo_combined_fasta.resolve().as_posix(), "-o",
                    diamond_search_out.resolve().as_posix(), "-f", "6", "qseqid", 
                    "stitle", "pident", "evalue", "mismatch"]

    print("DeNovo: Running the following command:\n")
    print(" ".join(diamond_args) + "\n")  
    # os.chdir(diamond_folder)
    # print(os.getcwd())
    ## Use subprocess
    os.system(" ".join(diamond_args))
    # os.chdir(working_dir)

species_tally_path = intermediate_dir / (prefix + f'_{db_name}{suffix}.xlsx')
detailed_fout = intermediate_dir / f'{prefix}_{db_name}_detailed.csv'
# Step 4. Tallying the diamond results
run_diamond_tally(diamond_search_out, 
                  int(config['taxa to fasta']['top_strains']), 
                  ncbi_taxa_folder, 
                  config['diamond tally']['mode'], 
                  species_tally_path, detailed_fout,
                  int(config['diamond tally']['n_protein_cutoff']),
                  config['diamond tally']['db_pattern'],
                  benchmark, config['diamond tally']['taxa_stats'])


## Step 5. Putting together the final fasta file.

if config['taxa to fasta']['kingdom_list'] != "":
    kingdom_list = config['taxa to fasta']['kingdom_list'].split(', ')
    print(kingdom_list)
else:
    kingdom_list = []

coverage_target = str(config['taxa to fasta']['coverage_target'])
output_fasta_path = final_dir / (prefix + f'_kaiko_fasta_coverage_{coverage_target}.fasta')
output_annotation_path = final_dir / (prefix + f'_kaiko_fasta_coverage_{coverage_target}_annotations.JSON')

aggregate_fasta(ref_fasta,
                ref_proteome_log,
                species_tally_path,
                output_fasta_path, output_annotation_path,
                config['taxa to fasta']['coverage_target'],
                int(config['taxa to fasta']['top_strains']),
                ref_fasta_igzip_index, index_path, index_s_path,
                kingdom_list, mode)

# aggregate_annotations()


if (config['denovo'])['profile']:
    profiler.enable()
    profiler.dump_stats('Kaiko_volume/Kaiko_taxa.prof')
