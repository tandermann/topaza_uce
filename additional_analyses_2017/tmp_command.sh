
for sim in analyses/species_da/simulated/10_reps/*;

do sim_rep=$(echo $sim | sed 's/.*\///g' | sed 's/\\r//g');
sed -i -e 's/:0.0;/;/g' analyses/treeannotator/$sim_rep/simulated_allele_alignments_stacey.tre
sed -i -e 's/:0.0;/;/g' analyses/treeannotator/$sim_rep/simulated_chimeric_allele_alignments_stacey.tre
sed -i -e 's/:0.0;/;/g' analyses/treeannotator/$sim_rep/simulated_consensus_contig_alignments_stacey.tre
sed -i -e 's/:0.0;/;/g' analyses/treeannotator/$sim_rep/simulated_iupac_consensus_alignments_stacey.tre
rm analyses/treeannotator/$sim_rep/*-e
done