#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 26 18:53:16 2017

@author: tobias
"""
import os
import re
import xml.etree.ElementTree as ET
import argparse

# Complete path function
class CompletePath(argparse.Action):
    """give the full path of an input file/folder"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))

# Get arguments
def get_args():
	parser = argparse.ArgumentParser(
		description="This script will replace the alignmnets in a BEAST xml file with alignmnets provided by the --alignments flag. Requires same number of alignments in xml file and in fasta folder",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	parser.add_argument(
		'--alignments',
		required=True,
		action=CompletePath,
		default=None,
		help='The directory containing selected fasta alignments'
	)
	parser.add_argument(
		'--xml_path',
		required=True,
		action=CompletePath,
		default=None,
		help='The directory containing the xml template file'
	)
	args = parser.parse_args()
	return args

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


# _______________Read the allele alignments___________________
#alignments = '/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/data/simulated/10_rep/selection/allele_alignments_simulated/rep1'
#xml_dict = '/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/analyses/stacey/simulated/10_reps/rep1/allele_alignments/'

args = get_args()
# Set working directory
alignments = args.alignments
xml_dict = args.xml_path

file_type = ''
path = alignments.split('/')
if 'rep' in path[-1]:
    file_type = path[-2]
else:
    file_type = path[-1]

seq_header_dict = {"f1_0":"f1^F","f1_1":"f2^F","d1_0":"d1^D","d1_1":"d2^D","d2_0":"d3^D","d2_1":"d4^D","e1_0":"e1^E","e1_1":"e2^E","e2_0":"e3^E","e2_1":"e4^E","x1_0":"x1^X","x1_1":"x2^X","x2_0":"x3^X","x2_1":"x4^X","y1_0":"y1^Y","y1_1":"y2^Y","z1_0":"z1^Z","z1_1":"z2^Z","z2_0":"z3^Z","z2_1":"z4^Z","f1":"F1","d1":"D1","d2":"D2","e1":"E1","e2":"E2","x1":"X1","x2":"X2","y1":"Y1","z1":"Z1","z2":"Z2"}

# Create a dictionary with the locus indeces and their corresponding file names
locus_file_name_dict = {}
counter = 0
for fasta in os.listdir(alignments):
    if fasta.endswith(".fasta") or fasta.endswith(".fa"):
        locus_file_name_dict.setdefault(counter,fasta)
        counter += 1

# read the BEAST xml file
tree = ET.ElementTree(file=os.path.join(xml_dict,'template.xml'))
root = tree.getroot()
counter = 0
for child in root:
    if child.tag == 'data':
        #locus = child.attrib['id']
        locus = counter
        # open the fasta file for that locus to retrieve the sequences
        file_name = locus_file_name_dict[locus]
        loc_name = re.sub('.fasta', '', file_name)
        #child.attrib['id'] = loc_name
        seq_dict = {}
        with open("%s/%s" %(alignments,file_name)) as f:
            for name, seq in read_fasta(f):
                name = re.sub('>', '', name)
                if name in seq_header_dict.keys():
                    name = seq_header_dict[name]
                seq_dict.setdefault(name,[])
                seq_dict[name].append(seq)
            for neighbor in child.iter('sequence'):
                taxon = neighbor.attrib['taxon']
                #taxon = re.sub('\^[A-Za-z]*', '', taxon)
                #taxon = re.sub('[1]', '1_0', taxon)
                #taxon = re.sub('[2]', '1_1', taxon)
                #taxon = re.sub('[3]', '2_0', taxon)
                #taxon = re.sub('[4]', '2_1', taxon)
                old_sequence = neighbor.attrib['value']
                #print(seq_dict.keys(),taxon)
                new_sequence = seq_dict[taxon][0]
                neighbor.attrib['value']=new_sequence
        counter += 1
                
tree.write(os.path.join(xml_dict,'%s.xml' %file_type))     
