#!/usr/bin/env python

import pysam
import pandas as pd
import csv


PROMETHION = pd.read_csv("/projects/rmorin/projects/gambl-repos/gambl-hshaalan/promethion_samples.tsv", sep = '\t')
paths = PROMETHION['link_name'].tolist()

chrz = list(range(1,23))
chromosomes = [str(int) for int in chrz]
chromosomes.append("X")
chr_prefixed = []
for chro in chromosomes:
    chr_prefix = "chr" + chro
    chr_prefixed.append(chr_prefix)

print(chr_prefixed)

def get_read_lengths(chrs, bam):

    readlengths = []
    for read in infile.fetch(chrs, until_eof = True, multiple_iterators = True):
        readlengths.append(read.query_length)
    return(readlengths)



with open("bam_readlengths_BL.tsv", "wt") as out_file:
    tsv_writer = csv.writer(out_file, delimiter = '\t')
    tsv_writer.writerow(['bam','read_length'])
    for bam in paths:
        infile = pysam.AlignmentFile(bam)
        lists = map(get_read_lengths, chr_prefixed, bam)
        all_lengths = [item for sublist in lists for item in sublist]
        mean = ((sum(all_lengths))/len(all_lengths))
        tsv_writer.writerow([bam,mean])

