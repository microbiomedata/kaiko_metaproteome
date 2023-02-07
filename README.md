# Kaiko Pipeline

## Introduction

Put simply, this tool takes raw proteomic input and outputs a FASTA file of those organisms most likely to be present in the proteomic input.

The pipeline uses neural networks to identify peptide sequences from raw proteomic input, which are then aligned against all protein sequences using a diamond search. This offers us a view of those organisms most likely to be present in the proteomic samples, with which we make a FASTA file from the most likely organisms identified.


## Setup

Before first use, a few database files are needed.

### Downloading Files

Download the following files to the ```Kaiko_volume/Kaiko_intermediate``` folder.

1) [UniRef100 FASTA](https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.fasta.gz)

2) [UniRef100 XML](https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.xml.gz)

3) [NCBI Taxanomy dump](https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref100/uniref100.xml.gz)

4) [Diamond search](https://github.com/bbuchfink/diamond/releases), choosing the appropriate system. If using Docker, get the Linux version.

### Processing

1) Extract the diamond file from step 4 into its own folder within ```Kaiko_volume/Kaiko_stationary```, eg ```Kaiko_volume/Kaiko_stationary/diamond```.

2) Within a command prompt, navigate to the diamond folder created in the previous step and run ```diamond makedb --in ../uniref100.fasta.gz --db ../uniref100```.

3) Extract the contents of NCBI Taxonomy dump to its own folder within ```Kaiko_volume/Kaiko_stationary```, eg ```Kaiko_volume/Kaiko_stationary/ncbi_taxa```.

4) 


## Usage

Currently, only .mgf files are supported. To use, simply follow these steps.

1) Place the input into a separate folder WITHIN the ```Kaiko_volume\Kaiko_input_files\``` directory. This folder should have a descriptive name. 

2) Edit the ```config.yaml``` file within the root directory of this repo to include the location of the folder with the input. An example can be found in the current file ```config.yaml```.

3) Run the command ``` python Kaiko_pipeline_main.py ```. The ```kaiko_defaults.yaml``` file will fill in any necessary parameters not present in ```config.yaml```


The ```Kaiko_volume/Kaiko_intermediate/``` folder will be populated with a few intermediate files. These are named using the ```mgf_input``` folder name. The final FASTA output can be found within ```Kaiko_volume/Kaiko_output/``` folder, again named using the folder name of the input.

4) If you would like to profile the pipeline using cProfile, add the profile = True flag.

## Usage with Docker

To use the pipeline within Docker, follow steps 1-2 above, then jump here:

3) (Docker) Run the command ```docker build -f Dockerfile_tensorflow1.2.1-py36``` to make the tensorflow image.

4) (Docker) Run the command ```docker build . -t kaiko-py36``` to build the Kaiko docker image using the tensorflow image from step 3)

5) (Docker) Run the command ```docker run --name Kaiko_container-py36 -v path_Kaiko_volume:/Kaiko_pipeline/Kaiko_volume -t -d kaiko-py36```, where path_Kaiko_volume is the absolute path to the Kaiko_volume folder. This allows Docker to store the outputs in Kaiko_volume. For example, such a command may look like ```docker run --name Kaiko_container-py36 -v C:/Users/memmys/Documents/GitHub/Kaiko_pipeline/Kaiko_volume/:/Kaiko_pipeline/Kaiko_volume -t -d kaiko-py36```

The ```Kaiko_volume/Kaiko_intermediate/``` folder will be populated with a few intermediate files. These are named using the ```mgf_input``` folder name. The final FASTA output can be found within ```Kaiko_volume/Kaiko_output/``` folder, again named using the folder name of the input.

