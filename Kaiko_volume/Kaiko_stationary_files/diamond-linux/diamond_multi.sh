#!/bin/bash

args=("$@")
for entry in "${args[0]}"/*.fasta
do
  echo "$entry"
  diamond blastp -d uniref100 --min-score 1 -q $entry -o ${entry/.fasta/.dmd} -f 6 qseqid stitle pident evalue mismatch
done