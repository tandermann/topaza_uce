#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 10:30:13 2017

@author: tobias
"""
import os
import re
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


wd = '/Users/tobias/Desktop/topaza_review/3rd_submission/new_analyses/finding_uce_contigs/blast_contigs_uces'
lastz_iupac = 'lastz_results_iupac.txt'
contig_file_name = 'abyss_contigs_T_pella5.fa'
match_output_name = 'contigs_matching_uces.fasta'
sample_id = 'T_pella5'

# read and define input and output files
lastz_file = os.path.join(wd,lastz_iupac) 
lastz_df = pd.read_csv(lastz_file,sep='\t')
# load the fasta formatted contig file
contig_file = os.path.join(wd,contig_file_name)
contig_sequences = SeqIO.parse(open(contig_file),'fasta')
# define the output file
match_output_file = os.path.join(wd,match_output_name)


# make a dictionary with all contig names that match a uce locus
uce_contig_dict = {}
contig_uce_dict = {}
for row in lastz_df.iterrows():
    locus = row[1].name2
    locus_name = re.sub('\_p[0-9]* \|.*', '', locus)
    locus_name = re.sub('^>', '', locus_name)
    contig_header = row[1].name1
    #print(contig_header)
    contig_name = re.sub('^\>([0-9]*) .*', '\\1', contig_header)
    #print(contig_name)
    uce_contig_dict.setdefault(locus_name,[])
    uce_contig_dict[locus_name].append(contig_name)
    contig_uce_dict.setdefault(contig_name,[])
    contig_uce_dict[contig_name].append(locus_name)

# get uces that have multiple contigs matching them
invalid_uce_loci = []
uces_with_multiple_hits = []
for uce in uce_contig_dict.keys():
    if len(uce_contig_dict[uce]) > 1:
        uces_with_multiple_hits.append(uce)
        invalid_uce_loci.append(uce)

# get uces that match on multiple contigs
contigs_matching_multiple_uces = []
for contig in contig_uce_dict.keys():
    if len(contig_uce_dict[contig]) > 1:
        contigs_matching_multiple_uces.append(contig)
        for uce in contig_uce_dict[contig]:
            invalid_uce_loci.append(uce)

# summarize all uces that should be excluded form further processing (duplicates)
invalid_uces_unique = list(set(invalid_uce_loci))
print(len(invalid_uces_unique), 'duplicate UCE loci had to be excluded')

# get list of valid contig names
valid_contig_names = []
for uce in uce_contig_dict:
    if uce not in invalid_uces_unique:
        contig_name = uce_contig_dict[uce]
        valid_contig_names.append(contig_name[0])

# extract valid contigs form contig file and print to fasta file with uce-names+ sample_id as headers
with open(match_output_file, "w") as out_file:
    for fasta in contig_sequences:
        if fasta.id in valid_contig_names:
            seq = fasta.seq
            # get the corresponding uce locus name from the dictionary
            header = '%s_%s' %(contig_uce_dict[fasta.id][0],sample_id)
            new_fasta = SeqRecord(seq, id=header, name='', description='')
            out_file.write(new_fasta.format('fasta'))




# for all valid uce loci:
# get the matching contig sequence
# export as fasta sequence in joined contig file (header = uce_name+sample_name, sequence = contig_sequence)
