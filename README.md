# Kaiko Pipeline

## Introduction

This is a tool which uses neural networks to denovo identify peptides from raw proteomic input, which are then aligned against all protein sequences using a diamond search. This offers us a view of those organisms most likely to be present in the proteomic samples. Finally, we make a .FASTA file from the most likely organisms identified.

Put simply, this tool takes raw proteomic input and outputs a FASTA file of those organisms most likely to be present in the proteomic input.


## Setup


## Usage

Currently, only .mgf files are supported. To use, simply follow these steps.

1) Place the input into a separate folder WITHIN the ```Kaiko_volume\Kaiko_input_files\``` directory. This folder should have a descriptive name. 

2) Edit the ```config.yaml``` file within the root directory of this repo to include the location of the folder with the input. An example can be found in the current file ```config.yaml```.

3) Run the command ``` python Kaiko_pipeline_main.py ```. The ```kaiko_defaults.yaml``` file will fill in any necessary parameters not present in ```config.yaml```