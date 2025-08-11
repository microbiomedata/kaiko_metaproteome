# Kaiko Pipeline

## Introduction

Put simply, this tool takes a dataset (as .raw proteomic input) and outputs a FASTA file of those organisms most likely to be present in the input. This FASTA is meant to be used as a protein database to search the .raw input with tools such as MSGF+.

In short, the pipeline uses a neural network to denovo sequence peptides from the raw proteomic input. These peptides are aligned against all protein sequences from the Uniprot Reference proteomes using DIAMOND. From these alignments, we produce a ranking of species most likely to have proteomes matching the spectra, and aggregate a FASTA consisting the top matching proteomes. Additionally, we output annotation containing pfam, EC, ko and cog annotations for the proteins in GFF format.

This pipeline is dockerized, with image available here: [here](https://hub.docker.com/r/camiloposso15/kaiko_2.0-py3.10)


## Setup

Runs using python 3.10 and pytorch 2.5. The full list of requirements can be found in ```Kaiko_volume/setup_libraries.txt```.

To use, the pipeline requires database files to be downloaded. The bulk of this setup can be done nearly automatically, but cannot be completed from docker.

## Database setup

NOTE: The database files can be obtained from [here](fakewebsitekejfbkwjbfkf). Download all the files into a folder `reference_protomes_db`.

Kaiko 2.0 relies on the reference proteomes from UniProt, along with their annotations. Set up is easy, in three steps! (Note: Once complete, the size of the database files is about 270 Gb). 

1) Download the source code for this pipeline from the release section. Then extract the contents to a folder of your choice. For this tutorial, we will call it `./Kaiko_metaproteome`.

2) Create a folder in which to store the database files, in this tutorial we will refer to this folder's path as `database_folder`. In a command prompt, navigate to the folder `./Kaiko_metaproteome/Kaiko_main/database_setup/` and run the following command: `python -m kaiko_fetch_annotations --out_dir database_folder --N_process 4 --cache_size 15`. This will start 4 separate, parallel processes to download both the FASTA for each proteome, and a JSON containing all the annotations of the proteins (pfam, EC, ko, etc). 

The download script will log the number of proteins in each proteome, as well as any issues with downloads. This log can be found in the `database_folder`. If the processor has many more cores, the value of `N_processes` can be raised to complete the setup faster. Note: If any proteomes fail to download, or the main process is stopped mid way, running the same command again will fetch only proteomes and annotations which have not been downloaded.

To check if there's discrepencies between the proteomes and annotations after the download, search (ctrl+F) for 'Failed integrity check' in the log. To check for any Proteomes which could not be found, search (ctrl+F) for 'Proteome not found'. 

3) Finally, once all the proteomes have downloaded successfully, download the DIAMOND alignment tool from the official github [here](https://github.com/bbuchfink/diamond/releases/tag/v2.1.11). Extract the file into the `database_folder` from step 2. Then, from a powershell prompt, navigate to the `database_folder` and run `cat *.fasta | .\diamond makedb -d reference_proteomes_db`. This will make a DIAMOND compatible database using all the sequences from the Reference Proteomes, which will be used to map denovo peptides to proteins.


## Usage

Currently, only .mgf files are supported. To run Kaiko, it is required that the database files be downloaded first.

1) Place the input mgf files into a folder within the ```Kaiko_metaproteome/Kaiko_volume/Kaiko_input_files/``` directory. The folder name should be unique to the dataset.

2) Edit the ```config.yaml``` file within the ```Kaiko_volume``` directory to include the location of the folder with the input. An example can be found in the current file ```config.yaml```. Ensure the 

3) From a command prompt, navigate to the `Kaiko_metaproteome` folder, and run the command ```python -m Kaiko_pipeline_main.py --config your_config.yaml```. 

The ```Kaiko_volume/Kaiko_output/``` will contain a subfolder with the same name as the input in step 1. Inside, we can find the The FASTA and GFF output for the dataset, as well intermediate files created by the pipeline.


## Usage with Docker

To use the pipeline within Docker, follow steps 1-4 in Database Setup. 

1) Create a folder named ```Kaiko_volume``` somewhere in your PC. Place the input mgf files into a subfolder inside. The subfolder name should be unique to the dataset. In the following steps, we will call the path to this folder ```absolute_path_Kaiko_volume```. 

2) Download a copy of the ```config.yaml``` file and place it in the folder ```Kaiko_volume``` from step 1. Edit this yaml file to include the location of the input mgf folder from step 1.

3) From a command prompt, run the following: ```docker run -v absolute_path_Kaiko_volume:/Kaiko_metaproteome/Kaiko_volume camiloposso15/kaiko_2.0-py3.10 python -m Kaiko_main.Kaiko_main --config Kaiko_volume/config_kansas_soil.yaml```

3) (Docker) Run the command ```docker build -f Dockerfile_tensorflow2.12.0-py310 -t tensorflow2.12.0-py310 .``` to make the tensorflow image.

4) (Docker) Run the command ```docker build . -t kaiko-py310``` to build the Kaiko docker image using the tensorflow image from step 3)

5) (Docker) Run the command ```docker run --name Kaiko_container-py310 -v path_Kaiko_volume:/Kaiko_pipeline/Kaiko_volume kaiko-py310 python Kaiko_pipeline_main.py```, where path_Kaiko_volume is the absolute path to the Kaiko_volume folder. This allows Docker to store the outputs in Kaiko_volume. For example, such a command may look like ```docker run --name Kaiko_container-py310 -v C:/Users/memmys/Documents/GitHub/Kaiko_pipeline/Kaiko_volume/:/Kaiko_pipeline/Kaiko_volume kaiko-py310 python Kaiko_pipeline_main.py```

6) (Docker) Make sure to update the config file to point to the Linux version of diamond. See the setup for more details.

The ```Kaiko_volume/Kaiko_intermediate/``` folder will be populated with a few intermediate files. These are named using the ```mgf_input``` folder name. The final FASTA output can be found within ```Kaiko_volume/Kaiko_output/``` folder, again named using the folder name of the input.


## Unit Tests

After installing the files, we should ensure the denovo network is producing the expected output given the model. To do this, navigate to the main repo folder in a command prompt and run ```python kaiko_unit_test.py```. This runs the denovo model on a predetermined dataset and compares line by line to stored output.

