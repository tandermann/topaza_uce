cd all_alignments/allele_alignments_simulated
for rep in rep*; do phyluce_align_convert_one_align_to_another --alignments $rep --input-format fasta --output nexus_$rep --output-format nexus; done
for rep in nexus_rep*; do phyluce_align_get_informative_sites --alignments $rep --input-format nexus > stats_$rep; done
for rep in nexus_rep*; do python2.7 ../../../../../bin/extract_uces_from_get_inf_sites_output.py --input stats_$rep --output selection_$rep --mode top --threshold 150; done
for rep in rep*; do mkdir -p ../../selection/allele_alignments_simulated/$rep; for locus in $(cat top_150_selection_nexus_$rep | sed 's/.nexus/.fasta/g'); do cp $rep/$(echo $locus| tr -d '\r') ../../selection/allele_alignments_simulated/$rep; done; done
for rep in rep*; do mkdir -p ../../selection/chimeric_allele_alignments_simulated/$rep; for locus in $(cat top_150_selection_nexus_$rep | sed 's/.nexus/.fasta/g'); do cp ../chimeric_allele_alignments_simulated/$rep/$(echo $locus| tr -d '\r'| sed 's/sim_allele_alignment/*/g') ../../selection/chimeric_allele_alignments_simulated/$rep; done; done
for rep in rep*; do mkdir -p ../../selection/contig_consensus_alignments_simulated/$rep; for locus in $(cat top_150_selection_nexus_$rep | sed 's/.nexus/.fasta/g'); do cp ../contig_consensus_alignments_simulated/$rep/$(echo $locus| tr -d '\r'| sed 's/sim_allele_alignment/*/g') ../../selection/contig_consensus_alignments_simulated/$rep; done; done
for rep in rep*; do mkdir -p ../../selection/iupac_consensus_alignments_simulated/$rep; for locus in $(cat top_150_selection_nexus_$rep | sed 's/.nexus/.fasta/g'); do cp ../iupac_consensus_alignments_simulated/$rep/$(echo $locus| tr -d '\r'| sed 's/sim_allele_alignment/*/g') ../../selection/iupac_consensus_alignments_simulated/$rep; done; done
rm -r *nexus*
rm phyluce*
cd ../../
