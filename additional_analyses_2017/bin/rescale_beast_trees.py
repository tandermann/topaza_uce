#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 09:32:04 2017

@author: tobias
"""
import os
import sys
import re
import io
import csv
import pandas as pd
import numpy as np
from Bio import Phylo

#_______________________________functions______________________________
def remove_burnin(log_file,burnin):
    no_burnin = log_file.iloc[int(len(log_file)*burnin):].copy()
    return no_burnin

def get_average_mutation_rate(rate_log_columns):
    mean_across_all_loci = np.mean([np.mean(rate_log_columns[x]) for x in rate_log_columns.columns])
    mcmc_array = np.array([np.mean(x[1][:]) for x in rate_log_columns.iterrows()])
    return mean_across_all_loci,mcmc_array

def read_beast_tree_file(treefile):
    reader = csv.reader(open(treefile, 'r'), delimiter='\n')
    reader = list(reader)
    header_index = ''
    for row in reader:
        if len(row) > 0 and row[0].startswith('tree STATE_0'):
            header_index = int(reader.index(row)-1)
            print('Length of header:', header_index, 'lines')
    header = reader[0:header_index]
    trees = reader[header_index+1:]
    return header,trees

def rescale_tree(tree_object,rescaling_factor):
    for clade in tree_object.find_clades():
        clade.branch_length = clade.branch_length*rescaling_factor
    return tree_object


#_____________________________input____________________________
#beast_tree_file = '/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/analyses/stacey/simulated/10_reps/rep1/allele_alignments/species.trees'
#log_file = '/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/analyses/stacey/simulated/10_reps/rep1/allele_alignments/sim_allele_alignment_2.log'
beast_tree_file = sys.argv[1]
log_file = sys.argv[2]
print('-----------')
print('Processing file',log_file)
workdir = '/'.join(os.path.abspath(beast_tree_file).split('/')[:-1])
#burnin=0.1

#____________________________workflow______________________________
# read the log file
log_content = pd.read_csv(log_file,sep='\t',comment='#')
# remove burnin
#log_content = remove_burnin(log_content,burnin)
# only extract columns containing logs of mutation rates
index_column = np.array(log_content.Sample)
rate_columns = log_content[log_content.columns[[x.startswith('clockRate') for x in log_content.columns]]]
# get the average mutation rate across all loci
average_mutation_rate,mcmc_array = get_average_mutation_rate(rate_columns)
# read in the trees file
header,trees = read_beast_tree_file(beast_tree_file)
# set up a new output file
rescaled_trees_file = os.path.join(workdir,'rescaled_species.trees')
tmp_file = os.path.join(workdir,'temp.txt')
rescaled_output = open(rescaled_trees_file, "w")
rescaled_trees_file_log=csv.writer(rescaled_output,delimiter='\n')
# print the header to the output file
for row in header:
    rescaled_trees_file_log.writerow(row)
rescaled_trees_file_log.writerow([';'])
# iterate through trees, rescale each by average mutation rate and print to file
print('Average clock rate across all mcmc iterations: %f. Rescaling species trees with individual clock rate for each mcmc iteration' %average_mutation_rate)
for row in trees:
    if not 'End;' in row:
        name = row[0].split('=')[0]
        state = int(name.split('_')[-1])
        # get the clock rate from this MCMC state
        rescaling_factor = mcmc_array[index_column==state]
        tree_object = Phylo.read(io.StringIO(row[0]), "newick")
        #Phylo.draw(tree_object)
        rescaled_tree_object = rescale_tree(tree_object,rescaling_factor)
        #Phylo.draw(rescaled_tree_object)
        Phylo.write(rescaled_tree_object,tmp_file,format='newick',format_branch_length='%1.15f')
        #tmp_file.flush()
        read_tmp = csv.reader(open(tmp_file, 'r'), delimiter='\n')
        tree = list(read_tmp)[0][0]
        tree = re.sub('=','',tree)
        new_tree = '= '.join([name,tree])
        rescaled_trees_file_log.writerow([new_tree])
rescaled_output.close()
