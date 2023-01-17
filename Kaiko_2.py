import os
import re
import pandas as pd

from math import floor

## This is a the same in function as the R script which combines the denovo output into a single fasta file.
## The only difference is that python and R handle ties in ordering differently. When a table is ordered,
## And the top 25% of the rows are taken, there can be slight differences between R and python.
directory = 'Kaiko_volume/Kaiko_intermediate/denovo_output/'
files = [f for f in os.listdir(directory) if bool(re.search(r'_out.txt', f))]
selection = 0.25
samples = []

def summary_times(group):
    scans = group['scan']
    return len(scans)

def summary_scans(group):
    scans = group['scan']
    return "_".join(scans)

def summary_rank(group):
    ranks = group['rank']
    return min(ranks)

def combine_denovo_output(directory, selection = 0.25):

    files = [f for f in os.listdir(directory) if bool(re.search(r'_out.txt', f))]
    samples = []

    for file in files:
        xx = pd.read_csv(directory + file, sep = "\t", header = 0)

        xx['output_seq'] = [re.sub(",", "", str(peptide)) for peptide in xx['output_seq']]
        xx['output_seq'] = [re.sub("mod", "", str(peptide)) for peptide in xx['output_seq']]
        xx = xx.loc[xx['output_score'] != float("Inf")]

        xx = xx.sort_values('output_score', ascending = False)
        xx = xx.head(floor(selection * floor(len(xx.index))))
        #xx = xx[['scan', 'output_seq']]

        xx['pep_length'] = [len(peptide) for peptide in xx['output_seq']]
        xx = xx.loc[(xx['pep_length'] >= 10) & (xx['pep_length'] <= 17)]
        xx['rank'] = list(range(1, len(xx.index) + 1))

        grouped = xx.groupby('output_seq')

        summary = grouped.apply(summary_times).to_frame()
        # summary = grouped.apply(summary_times)
        print(file)
        print(type(summary))
        print(summary)
        summary['output_seq'] = summary.index
        print(summary.columns)
        summary.columns = ['times', 'output_seq']
        summary = summary[['output_seq', 'times']]
        summary['rank'] = grouped.apply(summary_rank)
        summary['scans'] = grouped.apply(summary_scans)

        samples += [summary]


    for summary in samples:
        nms = [">S" + summary['scans'][i] + "_" + str(summary['times'][i]) for i in range(len(summary))]
        to_write = [None]*2*len(nms)
        to_write[::2] = nms
        to_write[1::2] = summary['output_seq']

        with open('Kaiko_volume/Kaiko_intermediate/combined_denovo.fasta', 'a') as fasta_file:
            for line in to_write:
                fasta_file.write(f"{line}\n")
                # fasta_file.write(line + "\n")


