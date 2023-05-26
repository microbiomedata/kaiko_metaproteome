import configparser
import yaml


## This contains all the possible flags which can be used.
config = {}

# config['DEFAULT'] = {'prefix' : ""}

config['denovo'] = {'topk' : False, 
                    'multi_decode' : True,
                    'beam_search' : True,
                    'profile' : False,
                    'beam_size' : 5,
                    'mgf_dir' : 'Kaiko_volume/Kaiko_input_files/',
                    'cached' : False,
                    # 'train_dir' : 'model/'
                    # 'decode_dir' : 'Kaiko_volume/Kaiko_intermediate/denovo_output/',
                    }


config['diamond tally'] = {'diamond_folder' : 'Kaiko_volume/Kaiko_stationary_files/diamond',
                           'ncbi_taxa_folder' : "Kaiko_volume/Kaiko_stationary_files/ncbi_taxa",
                           'mode' : 'member',
                          #  'fout' : 'Kaiko_volume/Kaiko_intermediate/kaiko_prediction_top_taxa.csv',
                          #  'diamond_output' : "Kaiko_volume/Kaiko_intermediate/denovo_output/diamond_search_output.dmd",
                           'n_protein_cutoff' : 300000,
                           'pident' : 100,
                           'cached' : False}


config['taxa to fasta'] = {'ref_fasta' : "Kaiko_volume/Kaiko_stationary_files/uniref100.fasta.gz",
                          #  'diamond_tally' : "Kaiko_volume/Kaiko_intermediate/kaiko_prediction_top_taxa.csv",
                          #  'fout' : "Kaiko_volume/Kaiko_output/kaiko_output.fasta",
                           'coverage_target' : 0.66,
                           'top_strains' : 1,
                           'taxa_key' : "TaxID",
                           'kingdom_list' : ""}

# with open('kaiko_defaults.ini', 'w') as configfile:
#   config.write(configfile)

f = open("kaiko_defaults.yaml", 'w')
yaml.safe_dump(config, f)



