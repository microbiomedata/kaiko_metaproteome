import configparser


## This contains all the possible flags which can be used.
config = configparser.ConfigParser()

# config['DEFAULT'] = {'prefix' : ""}

config['denovo'] = {'topk' : False, 
                    'multi_decode' : True,
                    'beam_search' : False,
                    'beam_size' : 5,
                    'mgf_dir' : 'pipeline_input/',
                    'train_dir' : 'model/',
                    'decode_dir' : 'pipeline_intermediary/denovo_output/'}

config['diamond tally'] = {'diamond_output' : "pipeline_intermediary/diamond_search_output.dmd",
                           'ntops' : -1,
                           'ncbi_taxa_folder' : "ncbi_taxa_2022-10-20",
                           'mode' : 'member',
                           'fout' : 'pipeline_intermediary/kaiko_prediction_top_taxa.csv',
                           'pident' : 100}



config['taxa to fasta'] = {'ref_fasta' : "uniref100.fasta.gz",
                           'diamond_tally' : "pipeline_intermediary/kaiko_prediction_top_taxa.csv",
                           'fout' : "pipeline_output/kaiko_output.fasta",
                           'ntops' : 5,
                           'taxa_key' : "TaxID",
                           'kingdom_list' : ""}

with open('kaiko_defaults.ini', 'w') as configfile:
  config.write(configfile)





