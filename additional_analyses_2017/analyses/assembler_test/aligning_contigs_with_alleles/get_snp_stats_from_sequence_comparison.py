#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 10:23:54 2017

@author: tobias
"""

import os
import re
from Bio import AlignIO


alignment_path = '/Users/tobias/GitHub/topaza_uce/assembler_test/aligning_contigs_with_alleles/sequence_comparisons/alignments'
alignment_list = os.listdir(alignment_path)

valid_snp_contig_list = []
all_pos_contig_list = []
loci_with_heterozygous_alleles = []
for fasta in alignment_list:
    fasta_path = os.path.join(alignment_path,fasta)
    alignment = AlignIO.read(open(fasta_path), "fasta")
    #def variable_positions(alignment):
    for x in range(alignment.get_alignment_length()):
        # get the position in the two allele sequences
        allele_base_calls = alignment[1:3,x]
        # get the contig base call
        abyss_contig_call = alignment[0,x]
        all_pos_contig_list.append((allele_base_calls,abyss_contig_call,re.sub('_sequence_alignment.fasta','',fasta),x))

        # get all 3 sequences of interest (abyss, and both alleles)
        col = alignment[0:3,x]
        column = list(col)
        nucleotides = set(column)
        # get only positions with basecalls present for all sequences
        #if "-" not in nucleotides and "?" not in nucleotides:
        # if there is a difference in basecall within alleles, examine
        if len(set(allele_base_calls))>1:
            # excluse positions where one of the two alleles has an ambiguous basecall ('n')
            if 'n' in list(allele_base_calls):
                pass
                #print('this',allele_base_calls)
            else:
                locus_name = re.sub('_sequence_alignment.fasta','',fasta)
                valid_snp_contig_list.append((allele_base_calls,abyss_contig_call,locus_name,x))
                if locus_name not in loci_with_heterozygous_alleles:
                    loci_with_heterozygous_alleles.append(locus_name)

# get all positions where IUPAC ambiguity codes are in contigs
iupac_contig_bases = []
loci_with_iupac_ambiguities = []
for snp in all_pos_contig_list:
    if snp[1] not in ['a','c','t','g','-']:
        #if snp[0] != 'nn':
        iupac_contig_bases.append(snp)
        if snp[2] not in loci_with_iupac_ambiguities:
            loci_with_iupac_ambiguities.append(snp[2])
     
valid_iupac_contig_bases = []
for snp in all_pos_contig_list:
    if snp[1] not in ['a','c','t','g','-']:
        if snp[0] not in ['nn','n-','-n','--']:
            valid_iupac_contig_bases.append(snp)

# how many heterozygous positions were called by allele phasing?        
print('In total',len(valid_snp_contig_list),'heterozygous positions in %i loci were found among allele sequences for this sample' %len(loci_with_heterozygous_alleles))
# how many where called in abyss contigs
print('In total',len(iupac_contig_bases),'ambiguous positions in %i loci were found in abyss contigs' %len(loci_with_iupac_ambiguities))
# how many of the abyss positions are valid (supported by >2 reads)?
print('In total',len(valid_iupac_contig_bases),'of the abyss snps are at positions that are also covered by allele seqeunces (no missing data)')
# overlap between heterozygous positions called by abyss and by allele phasing
shared_snps = set(iupac_contig_bases).intersection(valid_snp_contig_list)
print('In total',len(shared_snps),'heterozygous positions shared between abyss contigs and alleles')


