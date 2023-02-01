# Kaiko Pipeline

## Introduction

Put simply, this tool takes raw proteomic input and outputs a FASTA file of those organisms most likely to be present in the proteomic input.

The pipeline uses neural networks to identify peptide sequences from raw proteomic input, which are then aligned against all protein sequences using a diamond search. This offers us a view of those organisms most likely to be present in the proteomic samples, with which we make a FASTA file from the most likely organisms identified.


## Setup


## Usage

Currently, only .mgf files are supported. To use, simply follow these steps.

1) Place the input into a separate folder WITHIN the ```Kaiko_volume\Kaiko_input_files\``` directory. This folder should have a descriptive name. 

2) Edit the ```config.yaml``` file within the root directory of this repo to include the location of the folder with the input. An example can be found in the current file ```config.yaml```.

3) Run the command ``` python Kaiko_pipeline_main.py ```. The ```kaiko_defaults.yaml``` file will fill in any necessary parameters not present in ```config.yaml```


The ```Kaiko_volume\Kaiko_intermediate\``` folder will be populated with a few intermediate files. These are named using the ```mgf_input``` folder name. The final FASTA output can be found within ```Kaiko_volume\Kaiko_output\``` folder, again named using the folder name of the input.

4) If you would like to profile the pipeline using cProfile, add the profile = True flag.

