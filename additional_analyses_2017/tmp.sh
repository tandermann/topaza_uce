python bin/rescale_beast_trees.py analyses/stacey/empirical/allele_alignments/species.trees analyses/stacey/empirical/allele_alignments/*.log
python bin/rescale_beast_trees.py analyses/stacey/empirical/consensus_contig_alignments/species.trees analyses/stacey/empirical/consensus_contig_alignments/*.log
python bin/rescale_beast_trees.py analyses/stacey/empirical/chimeric_allele_alignments/species.trees analyses/stacey/empirical/chimeric_allele_alignments/*.log
python bin/rescale_beast_trees.py analyses/stacey/empirical/iupac_consensus_alignments/species.trees analyses/stacey/empirical/iupac_consensus_alignments/*.log
python bin/rescale_beast_trees.py analyses/stacey/simulated/allele_alignments/species.trees analyses/stacey/simulated/allele_alignments/*.log
python bin/rescale_beast_trees.py analyses/stacey/simulated/consensus_contig_alignments/species.trees analyses/stacey/simulated/consensus_contig_alignments/*.log
python bin/rescale_beast_trees.py analyses/stacey/simulated/chimeric_allele_alignments/species.trees analyses/stacey/simulated/chimeric_allele_alignments/*.log
python bin/rescale_beast_trees.py analyses/stacey/simulated/iupac_consensus_alignments/species.trees analyses/stacey/simulated/iupac_consensus_alignments/*.log
for sim in analyses/stacey/simulated/10_reps/*;
do echo $sim;
python bin/rescale_beast_trees.py $sim/allele_alignments/species.trees $sim/allele_alignments/*.log;
python bin/rescale_beast_trees.py $sim/consensus_contig_alignments/species.trees $sim/consensus_contig_alignments/*.log;
python bin/rescale_beast_trees.py $sim/chimeric_allele_alignments/species.trees $sim/chimeric_allele_alignments/*.log;
python bin/rescale_beast_trees.py $sim/iupac_consensus_alignments/species.trees $sim/iupac_consensus_alignments/*.log;
done;

for rep_result in $(ls analyses/stacey/simulated/10_reps);
do rep=$(basename $rep_result);
cp bin/extract_node_heights.py bin/tmp_extract_node_heights.py;
sed -i -e "s/xxxxx/10_reps\/$rep\/allele_alignments/g" bin/tmp_extract_node_heights.py;
sed -i -e "s/yyyyy/10_reps/g" bin/tmp_extract_node_heights.py;
python bin/tmp_extract_node_heights.py;
rm bin/tmp_extract_node_heights.py;
cp bin/extract_node_heights.py bin/tmp_extract_node_heights.py;
sed -i -e "s/xxxxx/10_reps\/$rep\/chimeric_allele_alignments/g" bin/tmp_extract_node_heights.py;
sed -i -e "s/yyyyy/10_reps/g" bin/tmp_extract_node_heights.py;
python bin/tmp_extract_node_heights.py;
rm bin/tmp_extract_node_heights.py;
cp bin/extract_node_heights.py bin/tmp_extract_node_heights.py;
sed -i -e "s/xxxxx/10_reps\/$rep\/consensus_contig_alignments/g" bin/tmp_extract_node_heights.py;
sed -i -e "s/yyyyy/10_reps/g" bin/tmp_extract_node_heights.py;
python bin/tmp_extract_node_heights.py;
rm bin/tmp_extract_node_heights.py;
cp bin/extract_node_heights.py bin/tmp_extract_node_heights.py;
sed -i -e "s/xxxxx/10_reps\/$rep\/iupac_consensus_alignments/g" bin/tmp_extract_node_heights.py;
sed -i -e "s/yyyyy/10_reps/g" bin/tmp_extract_node_heights.py;
python bin/tmp_extract_node_heights.py;
rm bin/tmp_extract_node_heights.py;
rm bin/*-e;
done

# plot the node height distribution for each simulation replicate
for rep_result in $(ls analyses/stacey/simulated/10_reps);
do rep=$(basename $rep_result);
repnum=$(echo $rep | sed 's/rep//g')
cp bin/plot_node_heigh_distribution_simreps.r bin/tmp_plot_node_heigh_distribution_simreps.r;
sed -i -e "s/xxxxx/$rep/g" bin/tmp_plot_node_heigh_distribution_simreps.r;
sed -i -e "s/yyyyy/$repnum/g" bin/tmp_plot_node_heigh_distribution_simreps.r;
Rscript bin/tmp_plot_node_heigh_distribution_simreps.r;
rm bin/tmp_plot_node_heigh_distribution_simreps.r;
rm bin/*-e;
done

pdfjoin results/node_depth_distribution/10_reps/*.pdf --outfile results/node_depth_distribution/10_reps/all_node_depths_simreps_combined.pdf
