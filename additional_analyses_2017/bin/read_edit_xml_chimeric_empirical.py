#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 18:53:16 2017

@author: tobias
"""
import os
import re
import xml.etree.ElementTree as ET

def read_fasta(fasta):
    name, seq = None, []
    for line in fasta:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


# _______________Read the chimeric allele alignments___________________
work_dir = '/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/data/empirical/selection/chimeric_allele_alignments_selection_topaza'

# Create a dictionary with the locus names and their corresponding file names
locus_file_name_dict = {} 
for fasta in os.listdir(work_dir):
    if fasta.endswith(".fasta") or fasta.endswith(".fa"):
        locus_name = re.sub('_chimeric','',fasta)
        locus_name = re.sub('_allele','_phased',locus_name)
        locus_name = re.sub('.fasta','',locus_name)
        locus_file_name_dict.setdefault(locus_name,fasta)
        
# read the BEAST xml file
tree = ET.ElementTree(file='/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/analyses/stacey/empirical/chimeric_allele_alignments/template.xml')
root = tree.getroot()
for child in root:
    if child.tag == 'data':
        locus = child.attrib['id']
        # open the fasta file for that locus to retrieve the sequences
        file_name = locus_file_name_dict[locus]
        seq_dict = {}
        with open("%s/%s" %(work_dir,file_name)) as f:
            for name, seq in read_fasta(f):
                name = re.sub('>', '', name)
                seq_dict.setdefault(name,[])
                seq_dict[name].append(seq)
            for neighbor in child.iter('sequence'):
                taxon = neighbor.attrib['taxon']
                old_sequence = neighbor.attrib['value']
                new_sequence = seq_dict[taxon][0]
                neighbor.attrib['value']=new_sequence
                
tree.write("/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/analyses/stacey/empirical/chimeric_allele_alignments/modified_topaza-150-top-uces-chimeric_alleles-stacey_no_outgroup.xml")     
