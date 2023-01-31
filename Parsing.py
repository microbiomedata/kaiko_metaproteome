import configparser
import yaml


## This contains all the possible flags which can be used.
config = {}

# config['DEFAULT'] = {'prefix' : ""}

config['denovo'] = {'topk' : False, 
                    'multi_decode' : True,
                    'beam_search' : False,
                    'profile' : False,
                    'beam_size' : 5,
                    'mgf_dir' : 'Kaiko_volume/Kaiko_input_files/',
                    # 'train_dir' : 'model/'
                    # 'decode_dir' : 'Kaiko_volume/Kaiko_intermediate/denovo_output/',
                    }


config['diamond tally'] = {'ntops' : -1,
                           'ncbi_taxa_folder' : "Kaiko_volume/Kaiko_stationary_files/ncbi_taxa_2022-10-20",
                           'mode' : 'member',
                          #  'fout' : 'Kaiko_volume/Kaiko_intermediate/kaiko_prediction_top_taxa.csv',
                          #  'diamond_output' : "Kaiko_volume/Kaiko_intermediate/denovo_output/diamond_search_output.dmd",
                           'pident' : 100}


config['taxa to fasta'] = {'ref_fasta' : "Kaiko_volume/Kaiko_stationary_files/uniref100.fasta.gz",
                          #  'diamond_tally' : "Kaiko_volume/Kaiko_intermediate/kaiko_prediction_top_taxa.csv",
                          #  'fout' : "Kaiko_volume/Kaiko_output/kaiko_output.fasta",
                           'ntops' : 5,
                           'taxa_key' : "TaxID",
                           'kingdom_list' : ""}

# with open('kaiko_defaults.ini', 'w') as configfile:
#   config.write(configfile)

f = open("kaiko_defaults.yaml", 'w')
yaml.safe_dump(config, f)



