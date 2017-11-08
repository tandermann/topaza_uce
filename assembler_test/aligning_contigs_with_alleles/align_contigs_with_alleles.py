#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  8 14:39:39 2017

@author: tobias
"""

import os
import re
import random
from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Alphabet import generic_dna

# get names of all allele alignments
allele_alignments_path = '/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/data/empirical/all_alignments/allele_alignments_topaza'
chimeric_alignments_path = '/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/data/empirical/all_alignments/chimeric_allele_alignments_topaza'
contig_consensus_path = '/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/data/empirical/all_alignments/contig_consensus_alignments_topaza'
iupac_consensus_path = '/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/data/empirical/all_alignments/iupac_consensus_alignments_topaza'
allele_alignment_files = os.listdir(allele_alignments_path)

taxon = 'T_pella5'

# read the target contig file
contig_file = '/Users/tobias/GitHub/topaza_uce/assembler_test/finding_uce_contigs/blast_contigs_uces/contigs_matching_uces.fasta'
contig_sequences = SeqIO.parse(open(contig_file),'fasta')

# pick n random loci for which to extract both allele sequences and the contig sequence
n = 10
# shuffle list randomely
random.shuffle(allele_alignment_files)
# take first n elements
random_pick = allele_alignment_files[0:n]

def get_seq_objects_for_taxon(taxon, fasta_path):
    file_name = fasta_path.split('/')[-1]
    suffix = '_' + file_name
    seq_list = []
    alignment = SeqIO.parse(open(fasta_path),'fasta')
    for seq in alignment:
        if seq.id.startswith(taxon):
            seq_list.append(seq)
            new_header = seq.id+suffix
            seq.id = re.sub('.fasta','',new_header)
            seq.name=''
            seq.description=''
    return seq_list

    
locus_fasta_dict = {}
for contig_fasta in contig_sequences:
    fasta_list = []
    locus = re.sub('^(uce-[0-9]*)_.*','\\1',contig_fasta.id)
    
    # compile the corresponding file names and check if locus is in our random selection
    allele_file_name = re.sub('^(uce-[0-9]*)','\\1_phased.fasta',locus)
    if allele_file_name in random_pick:
        chimeric_allele_file_name = re.sub('^(uce-[0-9]*)','\\1_chimeric_allele.fasta',locus)
        consensus_file_name = re.sub('^(uce-[0-9]*)','\\1.fasta',locus)
        iupac_cons_file_name = allele_file_name
    
        # get full path for each target alignment file
        final_allele_path = os.path.join(allele_alignments_path,allele_file_name)
        final_chimeric_path = os.path.join(chimeric_alignments_path,chimeric_allele_file_name)
        final_contig_cons_path = os.path.join(contig_consensus_path,consensus_file_name)
        final_iupac_cons_path = os.path.join(iupac_consensus_path,iupac_cons_file_name)
        
        # get fasta objects from file names        
        allele_fastas = get_seq_objects_for_taxon(taxon,final_allele_path)
        chimeric_fastas = get_seq_objects_for_taxon(taxon,final_chimeric_path)
        consensus_fastas = get_seq_objects_for_taxon(taxon,final_contig_cons_path)
        iupac_consensus_fastas = get_seq_objects_for_taxon(taxon,final_iupac_cons_path)
        
        fasta_list.append(contig_fasta)
        [fasta_list.append(item) for item in allele_fastas]
        #[fasta_list.append(item) for item in chimeric_fastas]
        #[fasta_list.append(item) for item in consensus_fastas]
        [fasta_list.append(item) for item in iupac_consensus_fastas]
        
        locus_fasta_dict.setdefault(locus,fasta_list)
        

out_path = '/Users/tobias/GitHub/topaza_uce/assembler_test/aligning_contigs_with_alleles/sequence_comparisons'
for locus in locus_fasta_dict:
    filename = '%s_sequences.fasta' %locus
    with open(os.path.join(out_path,filename), "w") as out_file:
        seq_list = locus_fasta_dict[locus]
        for sequence in seq_list:
            if '_0' in sequence.id:
                sequence.id = 'phased_allele0'
            elif '_1' in sequence.id:
                sequence.id = 'phased_allele1'
            elif sequence.id.startswith('T_pella5'):
                sequence.id = 'iupac_allele_consensus'
            else:
                sequence.id = 'abyss_contig'
                sequence.name=''
                sequence.description=''                
            out_file.write(sequence.format('fasta'))


"""
out_path = '/Users/tobias/GitHub/topaza_uce/assembler_test/aligning_contigs_with_alleles/sequence_comparisons'
for locus in locus_fasta_dict:
    filename = '%s_sequences.fasta' %locus
    seq_list = locus_fasta_dict[locus]
    align = MultipleSeqAlignment(seq_list, generic_dna)
    output_handle = open(os.path.join(out_path,filename), "w")
    AlignIO.write(align, output_handle, "fasta")
    output_handle.close()
"""
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        