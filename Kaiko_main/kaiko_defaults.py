import configparser
import yaml


## This contains all the possible flags which can be used.
config = {}

# config['DEFAULT'] = {'prefix' : ""}

config['general'] = {'output_dir' : 'Kaiko_volume'}

config['denovo'] = {'method' : "PNNL_Casanovo",
                    'topk' : False, 
                    'multi_decode' : True,
                    'beam_search' : True,
                    'profile' : False,
                    'beam_size' : 5,
                    'mgf_dir' : 'Kaiko_volume/Kaiko_input_files/',
                    'dms_analysis_job' : '',
                    'keep_dms_locally' : False,
                    'cached' : False,
                    # 'train_dir' : 'model/'
                    # 'decode_dir' : 'Kaiko_volume/Kaiko_intermediate/denovo_output/',
                    }


config['diamond tally'] = {'diamond_folder' : 'Kaiko_volume/Kaiko_stationary_files/diamond',
                           'ncbi_taxa_folder' : 'E:/Kaiko_stationary_files/ncbi_taxa',
                           'mode' : 'member',
                          #  'fout' : 'Kaiko_volume/Kaiko_intermediate/kaiko_prediction_top_taxa.csv',
                          #  'diamond_output' : "Kaiko_volume/Kaiko_intermediate/denovo_output/diamond_search_output.dmd",
                           'diamond_database' : 'E:/Kaiko_stationary_files/Reference_Proteomes/reference_proteomes_db',
                           'n_protein_cutoff' : 3000000,  ## Not important when using ref proteomes (default)
                           'cached' : False,
                           'db_pattern' : 'OX',
                           'benchmark' : [],
                           'taxa_stats' : 'Kaiko_volume/Kaiko_stationary_files/ncbi_taxa/uniref100_member_stats_with_lineage.txt'}


config['taxa to fasta'] = {'ref_fasta' : 'C:/Kaiko_stationary_files/Reference_Proteomes',
                          #  'diamond_tally' : "Kaiko_volume/Kaiko_intermediate/kaiko_prediction_top_taxa.csv",
                          #  'fout' : "Kaiko_volume/Kaiko_output/kaiko_output.fasta",
                           'gz_index' : 'Kaiko_volume/Kaiko_stationary_files/uniref100_fasta_gzindex.gzidx',
                           'proteome_index' : 'Kaiko_volume/Kaiko_stationary_files/uniref100_index.txt',
                           'proteome_index_s' : 'Kaiko_volume/Kaiko_stationary_files/uniref100_index_s.txt',
                           'ref_proteome_log' : 'Kaiko_volume/Kaiko_stationary_files/Reference Proteomes/database_log.xlsx',
                           'coverage_target' : 0.48,
                           'top_strains' : 1,
                           'taxa_key' : "TaxID",
                           'kingdom_list' : "",
                          # Always output the primary accession and the ID of the entry. These others are optional, but by default included.
                           'annotation_cats' : ['GO', 'KEGG', 'Reactome', 'verstion_date', 'version', 'gene_name', 'description']
                           }

# with open('kaiko_defaults.ini', 'w') as configfile:
#   config.write(configfile)

f = open("kaiko_defaults.yaml", 'w')
yaml.safe_dump(config, f)



