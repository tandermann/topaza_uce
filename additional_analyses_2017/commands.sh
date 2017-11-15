# The following workflow creates plots and sim-matrices from the results stored in the analyses folder

# _______________run speciesDA on all stacey results____________________
mkdir -p analyses/species_da/empirical/allele_alignments/
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1500 -collapseheight 1e-5 analyses/stacey/empirical/allele_alignments/species.trees analyses/species_da/empirical/allele_alignments/species_da_results_1e-5.txt
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1500 -collapseheight 1e-4 analyses/stacey/empirical/allele_alignments/species.trees analyses/species_da/empirical/allele_alignments/species_da_results_1e-4.txt
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1500 -collapseheight 3e-4 analyses/stacey/empirical/allele_alignments/species.trees analyses/species_da/empirical/allele_alignments/species_da_results_3e-4.txt
mkdir -p analyses/species_da/empirical/consensus_contig_alignments/
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1500 -collapseheight 1e-5 analyses/stacey/empirical/consensus_contig_alignments/species.trees analyses/species_da/empirical/consensus_contig_alignments/species_da_results_1e-5.txt
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1500 -collapseheight 1e-4 analyses/stacey/empirical/consensus_contig_alignments/species.trees analyses/species_da/empirical/consensus_contig_alignments/species_da_results_1e-4.txt
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1500 -collapseheight 3e-4 analyses/stacey/empirical/consensus_contig_alignments/species.trees analyses/species_da/empirical/consensus_contig_alignments/species_da_results_3e-4.txt

mkdir -p analyses/species_da/empirical/chimeric_allele_alignments/
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1500 -collapseheight 1e-5 analyses/stacey/empirical/chimeric_allele_alignments/species.trees analyses/species_da/empirical/chimeric_allele_alignments/species_da_results_1e-5.txt
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1500 -collapseheight 1e-4 analyses/stacey/empirical/chimeric_allele_alignments/species.trees analyses/species_da/empirical/chimeric_allele_alignments/species_da_results_1e-4.txt
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1500 -collapseheight 3e-4 analyses/stacey/empirical/chimeric_allele_alignments/species.trees analyses/species_da/empirical/chimeric_allele_alignments/species_da_results_3e-4.txt
mkdir -p analyses/species_da/empirical/iupac_consensus_alignments/
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1500 -collapseheight 1e-5 analyses/stacey/empirical/iupac_consensus_alignments/species.trees analyses/species_da/empirical/iupac_consensus_alignments/species_da_results_1e-5.txt
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1500 -collapseheight 1e-4 analyses/stacey/empirical/iupac_consensus_alignments/species.trees analyses/species_da/empirical/iupac_consensus_alignments/species_da_results_1e-4.txt
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1500 -collapseheight 3e-4 analyses/stacey/empirical/iupac_consensus_alignments/species.trees analyses/species_da/empirical/iupac_consensus_alignments/species_da_results_3e-4.txt


mkdir -p analyses/species_da/simulated/allele_alignments/
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1000 -collapseheight 1e-5 analyses/stacey/simulated/allele_alignments/species.trees analyses/species_da/simulated/allele_alignments/species_da_results_1e-5.txt
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1000 -collapseheight 1e-4 analyses/stacey/simulated/allele_alignments/species.trees analyses/species_da/simulated/allele_alignments/species_da_results_1e-4.txt
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1000 -collapseheight 3e-4 analyses/stacey/simulated/allele_alignments/species.trees analyses/species_da/simulated/allele_alignments/species_da_results_3e-4.txt
mkdir -p analyses/species_da/simulated/consensus_contig_alignments/
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1000 -collapseheight 1e-5 analyses/stacey/simulated/consensus_contig_alignments/species.trees analyses/species_da/simulated/consensus_contig_alignments/species_da_results_1e-5.txt
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1000 -collapseheight 1e-4 analyses/stacey/simulated/consensus_contig_alignments/species.trees analyses/species_da/simulated/consensus_contig_alignments/species_da_results_1e-4.txt
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1000 -collapseheight 3e-4 analyses/stacey/simulated/consensus_contig_alignments/species.trees analyses/species_da/simulated/consensus_contig_alignments/species_da_results_3e-4.txt

mkdir -p analyses/species_da/simulated/chimeric_allele_alignments/
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1000 -collapseheight 1e-5 analyses/stacey/simulated/chimeric_allele_alignments/species.trees analyses/species_da/simulated/chimeric_allele_alignments/species_da_results_1e-5.txt
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1000 -collapseheight 1e-4 analyses/stacey/simulated/chimeric_allele_alignments/species.trees analyses/species_da/simulated/chimeric_allele_alignments/species_da_results_1e-4.txt
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1000 -collapseheight 3e-4 analyses/stacey/simulated/chimeric_allele_alignments/species.trees analyses/species_da/simulated/chimeric_allele_alignments/species_da_results_3e-4.txt
mkdir -p analyses/species_da/simulated/iupac_consensus_alignments/
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1000 -collapseheight 1e-5 analyses/stacey/simulated/iupac_consensus_alignments/species.trees analyses/species_da/simulated/iupac_consensus_alignments/species_da_results_1e-5.txt
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1000 -collapseheight 1e-4 analyses/stacey/simulated/iupac_consensus_alignments/species.trees analyses/species_da/simulated/iupac_consensus_alignments/species_da_results_1e-4.txt
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1000 -collapseheight 3e-4 analyses/stacey/simulated/iupac_consensus_alignments/species.trees analyses/species_da/simulated/iupac_consensus_alignments/species_da_results_3e-4.txt


for sim in analyses/stacey/simulated/10_reps/*;
do mkdir -p analyses/species_da/simulated/10_reps/$(echo $(echo $sim | sed 's/.*\///g') | sed 's/.*\///g')/allele_alignments/;
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1000 -collapseheight 1e-5 analyses/stacey/simulated/10_reps/$(echo $sim | sed 's/.*\///g')/allele_alignments/species.trees analyses/species_da/simulated/10_reps/$(echo $sim | sed 's/.*\///g')/allele_alignments/species_da_results_1e-5.txt;
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1000 -collapseheight 1e-4 analyses/stacey/simulated/10_reps/$(echo $sim | sed 's/.*\///g')/allele_alignments/species.trees analyses/species_da/simulated/10_reps/$(echo $sim | sed 's/.*\///g')/allele_alignments/species_da_results_1e-4.txt;
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1000 -collapseheight 3e-4 analyses/stacey/simulated/10_reps/$(echo $sim | sed 's/.*\///g')/allele_alignments/species.trees analyses/species_da/simulated/10_reps/$(echo $sim | sed 's/.*\///g')/allele_alignments/species_da_results_3e-4.txt;
mkdir -p analyses/species_da/simulated/10_reps/$(echo $sim | sed 's/.*\///g')/consensus_contig_alignments/;
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1000 -collapseheight 1e-5 analyses/stacey/simulated/10_reps/$(echo $sim | sed 's/.*\///g')/consensus_contig_alignments/species.trees analyses/species_da/simulated/10_reps/$(echo $sim | sed 's/.*\///g')/consensus_contig_alignments/species_da_results_1e-5.txt;
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1000 -collapseheight 1e-4 analyses/stacey/simulated/10_reps/$(echo $sim | sed 's/.*\///g')/consensus_contig_alignments/species.trees analyses/species_da/simulated/10_reps/$(echo $sim | sed 's/.*\///g')/consensus_contig_alignments/species_da_results_1e-4.txt;
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1000 -collapseheight 3e-4 analyses/stacey/simulated/10_reps/$(echo $sim | sed 's/.*\///g')/consensus_contig_alignments/species.trees analyses/species_da/simulated/10_reps/$(echo $sim | sed 's/.*\///g')/consensus_contig_alignments/species_da_results_3e-4.txt;

mkdir -p analyses/species_da/simulated/10_reps/$(echo $sim | sed 's/.*\///g')/chimeric_allele_alignments/;
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1000 -collapseheight 1e-5 analyses/stacey/simulated/10_reps/$(echo $sim | sed 's/.*\///g')/chimeric_allele_alignments/species.trees analyses/species_da/simulated/10_reps/$(echo $sim | sed 's/.*\///g')/chimeric_allele_alignments/species_da_results_1e-5.txt;
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1000 -collapseheight 1e-4 analyses/stacey/simulated/10_reps/$(echo $sim | sed 's/.*\///g')/chimeric_allele_alignments/species.trees analyses/species_da/simulated/10_reps/$(echo $sim | sed 's/.*\///g')/chimeric_allele_alignments/species_da_results_1e-4.txt;
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1000 -collapseheight 3e-4 analyses/stacey/simulated/10_reps/$(echo $sim | sed 's/.*\///g')/chimeric_allele_alignments/species.trees analyses/species_da/simulated/10_reps/$(echo $sim | sed 's/.*\///g')/chimeric_allele_alignments/species_da_results_3e-4.txt;
mkdir -p analyses/species_da/simulated/10_reps/$(echo $sim | sed 's/.*\///g')/iupac_consensus_alignments/;
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1000 -collapseheight 1e-5 analyses/stacey/simulated/10_reps/$(echo $sim | sed 's/.*\///g')/iupac_consensus_alignments/species.trees analyses/species_da/simulated/10_reps/$(echo $sim | sed 's/.*\///g')/iupac_consensus_alignments/species_da_results_1e-5.txt;
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1000 -collapseheight 1e-4 analyses/stacey/simulated/10_reps/$(echo $sim | sed 's/.*\///g')/iupac_consensus_alignments/species.trees analyses/species_da/simulated/10_reps/$(echo $sim | sed 's/.*\///g')/iupac_consensus_alignments/species_da_results_1e-4.txt;
java -jar /Users/tobias/bin/speciesDA.jar -burnin 1000 -collapseheight 3e-4 analyses/stacey/simulated/10_reps/$(echo $sim | sed 's/.*\///g')/iupac_consensus_alignments/species.trees analyses/species_da/simulated/10_reps/$(echo $sim | sed 's/.*\///g')/iupac_consensus_alignments/species_da_results_3e-4.txt;

done;


# _________________plot similarity matrices_____________________

mkdir -p results/simmatrix_plots

# EMPIRICAL DATA_______________________
# empirical allele sequences
cp bin/template_similarity_matrix_empirical_allele.R bin/temp_similarity_matrix_empirical_allele.R
sed -i -e 's/xxx/allele_alignments/g' bin/temp_similarity_matrix_empirical_allele.R
sed -i -e 's/xe-x/1e-5/g' bin/temp_similarity_matrix_empirical_allele.R
Rscript bin/temp_similarity_matrix_empirical_allele.R
rm bin/temp_*
cp bin/template_similarity_matrix_empirical_allele.R bin/temp_similarity_matrix_empirical_allele.R
sed -i -e 's/xxx/allele_alignments/g' bin/temp_similarity_matrix_empirical_allele.R
sed -i -e 's/xe-x/1e-4/g' bin/temp_similarity_matrix_empirical_allele.R
Rscript bin/temp_similarity_matrix_empirical_allele.R
rm bin/temp_*
cp bin/template_similarity_matrix_empirical_allele.R bin/temp_similarity_matrix_empirical_allele.R
sed -i -e 's/xxx/allele_alignments/g' bin/temp_similarity_matrix_empirical_allele.R
sed -i -e 's/xe-x/3e-4/g' bin/temp_similarity_matrix_empirical_allele.R
Rscript bin/temp_similarity_matrix_empirical_allele.R
rm bin/temp_*

# empirical chimeric allele sequences
cp bin/template_similarity_matrix_empirical_allele.R bin/temp_similarity_matrix_empirical_allele.R
sed -i -e 's/xxx/chimeric_allele_alignments/g' bin/temp_similarity_matrix_empirical_allele.R
sed -i -e 's/xe-x/1e-5/g' bin/temp_similarity_matrix_empirical_allele.R
Rscript bin/temp_similarity_matrix_empirical_allele.R
rm bin/temp_*
cp bin/template_similarity_matrix_empirical_allele.R bin/temp_similarity_matrix_empirical_allele.R
sed -i -e 's/xxx/chimeric_allele_alignments/g' bin/temp_similarity_matrix_empirical_allele.R
sed -i -e 's/xe-x/1e-4/g' bin/temp_similarity_matrix_empirical_allele.R
Rscript bin/temp_similarity_matrix_empirical_allele.R
rm bin/temp_*
cp bin/template_similarity_matrix_empirical_allele.R bin/temp_similarity_matrix_empirical_allele.R
sed -i -e 's/xxx/chimeric_allele_alignments/g' bin/temp_similarity_matrix_empirical_allele.R
sed -i -e 's/xe-x/3e-4/g' bin/temp_similarity_matrix_empirical_allele.R
Rscript bin/temp_similarity_matrix_empirical_allele.R
rm bin/temp_*

# empirical consensus contig sequences
cp bin/template_similarity_matrix_empirical_consensus.R bin/temp_similarity_matrix_empirical_consensus.R
sed -i -e 's/xxx/consensus_contig_alignments/g' bin/temp_similarity_matrix_empirical_consensus.R
sed -i -e 's/xe-x/1e-5/g' bin/temp_similarity_matrix_empirical_consensus.R
Rscript bin/temp_similarity_matrix_empirical_consensus.R
rm bin/temp_*
cp bin/template_similarity_matrix_empirical_consensus.R bin/temp_similarity_matrix_empirical_consensus.R
sed -i -e 's/xxx/consensus_contig_alignments/g' bin/temp_similarity_matrix_empirical_consensus.R
sed -i -e 's/xe-x/1e-4/g' bin/temp_similarity_matrix_empirical_consensus.R
Rscript bin/temp_similarity_matrix_empirical_consensus.R
rm bin/temp_*
cp bin/template_similarity_matrix_empirical_consensus.R bin/temp_similarity_matrix_empirical_consensus.R
sed -i -e 's/xxx/consensus_contig_alignments/g' bin/temp_similarity_matrix_empirical_consensus.R
sed -i -e 's/xe-x/3e-4/g' bin/temp_similarity_matrix_empirical_consensus.R
Rscript bin/temp_similarity_matrix_empirical_consensus.R
rm bin/temp_*

# empirical consensus iupac sequences
cp bin/template_similarity_matrix_empirical_consensus.R bin/temp_similarity_matrix_empirical_consensus.R
sed -i -e 's/xxx/iupac_consensus_alignments/g' bin/temp_similarity_matrix_empirical_consensus.R
sed -i -e 's/xe-x/1e-5/g' bin/temp_similarity_matrix_empirical_consensus.R
Rscript bin/temp_similarity_matrix_empirical_consensus.R
rm bin/temp_*
cp bin/template_similarity_matrix_empirical_consensus.R bin/temp_similarity_matrix_empirical_consensus.R
sed -i -e 's/xxx/iupac_consensus_alignments/g' bin/temp_similarity_matrix_empirical_consensus.R
sed -i -e 's/xe-x/1e-4/g' bin/temp_similarity_matrix_empirical_consensus.R
Rscript bin/temp_similarity_matrix_empirical_consensus.R
rm bin/temp_*
cp bin/template_similarity_matrix_empirical_consensus.R bin/temp_similarity_matrix_empirical_consensus.R
sed -i -e 's/xxx/iupac_consensus_alignments/g' bin/temp_similarity_matrix_empirical_consensus.R
sed -i -e 's/xe-x/3e-4/g' bin/temp_similarity_matrix_empirical_consensus.R
Rscript bin/temp_similarity_matrix_empirical_consensus.R
rm bin/temp_*


# SIMULATED DATA_______________________
# simulated allele sequences
cp bin/template_similarity_matrix_simulated_allele.R bin/temp_similarity_matrix_simulated_allele.R
sed -i -e 's/xxx/allele_alignments/g' bin/temp_similarity_matrix_simulated_allele.R
sed -i -e 's/xe-x/1e-5/g' bin/temp_similarity_matrix_simulated_allele.R
Rscript bin/temp_similarity_matrix_simulated_allele.R
rm bin/temp_*
cp bin/template_similarity_matrix_simulated_allele.R bin/temp_similarity_matrix_simulated_allele.R
sed -i -e 's/xxx/allele_alignments/g' bin/temp_similarity_matrix_simulated_allele.R
sed -i -e 's/xe-x/1e-4/g' bin/temp_similarity_matrix_simulated_allele.R
Rscript bin/temp_similarity_matrix_simulated_allele.R
rm bin/temp_*
cp bin/template_similarity_matrix_simulated_allele.R bin/temp_similarity_matrix_simulated_allele.R
sed -i -e 's/xxx/allele_alignments/g' bin/temp_similarity_matrix_simulated_allele.R
sed -i -e 's/xe-x/3e-4/g' bin/temp_similarity_matrix_simulated_allele.R
Rscript bin/temp_similarity_matrix_simulated_allele.R
rm bin/temp_*

# simulated chimeric allele sequences
cp bin/template_similarity_matrix_simulated_allele.R bin/temp_similarity_matrix_simulated_allele.R
sed -i -e 's/xxx/chimeric_allele_alignments/g' bin/temp_similarity_matrix_simulated_allele.R
sed -i -e 's/xe-x/1e-5/g' bin/temp_similarity_matrix_simulated_allele.R
Rscript bin/temp_similarity_matrix_simulated_allele.R
rm bin/temp_*
cp bin/template_similarity_matrix_simulated_allele.R bin/temp_similarity_matrix_simulated_allele.R
sed -i -e 's/xxx/chimeric_allele_alignments/g' bin/temp_similarity_matrix_simulated_allele.R
sed -i -e 's/xe-x/1e-4/g' bin/temp_similarity_matrix_simulated_allele.R
Rscript bin/temp_similarity_matrix_simulated_allele.R
rm bin/temp_*
cp bin/template_similarity_matrix_simulated_allele.R bin/temp_similarity_matrix_simulated_allele.R
sed -i -e 's/xxx/chimeric_allele_alignments/g' bin/temp_similarity_matrix_simulated_allele.R
sed -i -e 's/xe-x/3e-4/g' bin/temp_similarity_matrix_simulated_allele.R
Rscript bin/temp_similarity_matrix_simulated_allele.R
rm bin/temp_*

# simulated consensus contig sequences
cp bin/template_similarity_matrix_simulated_consensus.R bin/temp_similarity_matrix_simulated_consensus.R
sed -i -e 's/xxx/consensus_contig_alignments/g' bin/temp_similarity_matrix_simulated_consensus.R
sed -i -e 's/xe-x/1e-5/g' bin/temp_similarity_matrix_simulated_consensus.R
Rscript bin/temp_similarity_matrix_simulated_consensus.R
rm bin/temp_*
cp bin/template_similarity_matrix_simulated_consensus.R bin/temp_similarity_matrix_simulated_consensus.R
sed -i -e 's/xxx/consensus_contig_alignments/g' bin/temp_similarity_matrix_simulated_consensus.R
sed -i -e 's/xe-x/1e-4/g' bin/temp_similarity_matrix_simulated_consensus.R
Rscript bin/temp_similarity_matrix_simulated_consensus.R
rm bin/temp_*
cp bin/template_similarity_matrix_simulated_consensus.R bin/temp_similarity_matrix_simulated_consensus.R
sed -i -e 's/xxx/consensus_contig_alignments/g' bin/temp_similarity_matrix_simulated_consensus.R
sed -i -e 's/xe-x/3e-4/g' bin/temp_similarity_matrix_simulated_consensus.R
Rscript bin/temp_similarity_matrix_simulated_consensus.R
rm bin/temp_*

# simulated consensus iupac sequences
cp bin/template_similarity_matrix_simulated_consensus.R bin/temp_similarity_matrix_simulated_consensus.R
sed -i -e 's/xxx/iupac_consensus_alignments/g' bin/temp_similarity_matrix_simulated_consensus.R
sed -i -e 's/xe-x/1e-5/g' bin/temp_similarity_matrix_simulated_consensus.R
Rscript bin/temp_similarity_matrix_simulated_consensus.R
rm bin/temp_*
cp bin/template_similarity_matrix_simulated_consensus.R bin/temp_similarity_matrix_simulated_consensus.R
sed -i -e 's/xxx/iupac_consensus_alignments/g' bin/temp_similarity_matrix_simulated_consensus.R
sed -i -e 's/xe-x/1e-4/g' bin/temp_similarity_matrix_simulated_consensus.R
Rscript bin/temp_similarity_matrix_simulated_consensus.R
rm bin/temp_*
cp bin/template_similarity_matrix_simulated_consensus.R bin/temp_similarity_matrix_simulated_consensus.R
sed -i -e 's/xxx/iupac_consensus_alignments/g' bin/temp_similarity_matrix_simulated_consensus.R
sed -i -e 's/xe-x/3e-4/g' bin/temp_similarity_matrix_simulated_consensus.R
Rscript bin/temp_similarity_matrix_simulated_consensus.R
rm bin/temp_*


# calculate simmatrices for the 10 simulation replicates
for sim in analyses/species_da/simulated/10_reps/*;
do sim_rep=$(echo $sim | sed 's/.*\///g' | sed 's/\\r//g');
mkdir results/simmatrix_plots/$sim_rep

# simulated allele sequences
cp bin/template_similarity_matrix_simulated_allele.R bin/temp_similarity_matrix_simulated_allele.R
sed -i -e "s/simulated\/xxx/simulated\/10_reps\/$sim_rep\/allele_alignments/g" bin/temp_similarity_matrix_simulated_allele.R
sed -i -e 's/xe-x/1e-5/g' bin/temp_similarity_matrix_simulated_allele.R
sed -i -e "s/\"\/simmatrix_simulated_xxx/\"$sim_rep\/simmatrix_simulated_allele_alignments/g" bin/temp_similarity_matrix_simulated_allele.R
Rscript bin/temp_similarity_matrix_simulated_allele.R
rm bin/temp_*
cp bin/template_similarity_matrix_simulated_allele.R bin/temp_similarity_matrix_simulated_allele.R
sed -i -e "s/simulated\/xxx/simulated\/10_reps\/$sim_rep\/allele_alignments/g" bin/temp_similarity_matrix_simulated_allele.R
sed -i -e 's/xe-x/1e-4/g' bin/temp_similarity_matrix_simulated_allele.R
sed -i -e "s/\"\/simmatrix_simulated_xxx/\"$sim_rep\/simmatrix_simulated_allele_alignments/g" bin/temp_similarity_matrix_simulated_allele.R
Rscript bin/temp_similarity_matrix_simulated_allele.R
rm bin/temp_*
cp bin/template_similarity_matrix_simulated_allele.R bin/temp_similarity_matrix_simulated_allele.R
sed -i -e "s/simulated\/xxx/simulated\/10_reps\/$sim_rep\/allele_alignments/g" bin/temp_similarity_matrix_simulated_allele.R
sed -i -e 's/xe-x/3e-4/g' bin/temp_similarity_matrix_simulated_allele.R
sed -i -e "s/\"\/simmatrix_simulated_xxx/\"$sim_rep\/simmatrix_simulated_allele_alignments/g" bin/temp_similarity_matrix_simulated_allele.R
Rscript bin/temp_similarity_matrix_simulated_allele.R
rm bin/temp_*


# simulated chimeric allele sequences
cp bin/template_similarity_matrix_simulated_allele.R bin/temp_similarity_matrix_simulated_allele.R
sed -i -e "s/simulated\/xxx/simulated\/10_reps\/$sim_rep\/chimeric_allele_alignments/g" bin/temp_similarity_matrix_simulated_allele.R
sed -i -e 's/xe-x/1e-5/g' bin/temp_similarity_matrix_simulated_allele.R
sed -i -e "s/\"\/simmatrix_simulated_xxx/\"$sim_rep\/simmatrix_simulated_chimeric_allele_alignments/g" bin/temp_similarity_matrix_simulated_allele.R
Rscript bin/temp_similarity_matrix_simulated_allele.R
rm bin/temp_*
cp bin/template_similarity_matrix_simulated_allele.R bin/temp_similarity_matrix_simulated_allele.R
sed -i -e "s/simulated\/xxx/simulated\/10_reps\/$sim_rep\/chimeric_allele_alignments/g" bin/temp_similarity_matrix_simulated_allele.R
sed -i -e 's/xe-x/1e-4/g' bin/temp_similarity_matrix_simulated_allele.R
sed -i -e "s/\"\/simmatrix_simulated_xxx/\"$sim_rep\/simmatrix_simulated_chimeric_allele_alignments/g" bin/temp_similarity_matrix_simulated_allele.R
Rscript bin/temp_similarity_matrix_simulated_allele.R
rm bin/temp_*
cp bin/template_similarity_matrix_simulated_allele.R bin/temp_similarity_matrix_simulated_allele.R
sed -i -e "s/simulated\/xxx/simulated\/10_reps\/$sim_rep\/chimeric_allele_alignments/g" bin/temp_similarity_matrix_simulated_allele.R
sed -i -e 's/xe-x/3e-4/g' bin/temp_similarity_matrix_simulated_allele.R
sed -i -e "s/\"\/simmatrix_simulated_xxx/\"$sim_rep\/simmatrix_simulated_chimeric_allele_alignments/g" bin/temp_similarity_matrix_simulated_allele.R
Rscript bin/temp_similarity_matrix_simulated_allele.R
rm bin/temp_*

# simulated consensus contig sequences
cp bin/template_similarity_matrix_simulated_consensus.R bin/temp_similarity_matrix_simulated_consensus.R
sed -i -e "s/simulated\/xxx/simulated\/10_reps\/$sim_rep\/consensus_contig_alignments/g" bin/temp_similarity_matrix_simulated_consensus.R
sed -i -e 's/xe-x/1e-5/g' bin/temp_similarity_matrix_simulated_consensus.R
sed -i -e "s/\"\/simmatrix_simulated_xxx/\"$sim_rep\/simmatrix_simulated_consensus_contig_alignments/g" bin/temp_similarity_matrix_simulated_consensus.R
Rscript bin/temp_similarity_matrix_simulated_consensus.R
rm bin/temp_*
cp bin/template_similarity_matrix_simulated_consensus.R bin/temp_similarity_matrix_simulated_consensus.R
sed -i -e "s/simulated\/xxx/simulated\/10_reps\/$sim_rep\/consensus_contig_alignments/g" bin/temp_similarity_matrix_simulated_consensus.R
sed -i -e 's/xe-x/1e-4/g' bin/temp_similarity_matrix_simulated_consensus.R
sed -i -e "s/\"\/simmatrix_simulated_xxx/\"$sim_rep\/simmatrix_simulated_consensus_contig_alignments/g" bin/temp_similarity_matrix_simulated_consensus.R
Rscript bin/temp_similarity_matrix_simulated_consensus.R
rm bin/temp_*
cp bin/template_similarity_matrix_simulated_consensus.R bin/temp_similarity_matrix_simulated_consensus.R
sed -i -e "s/simulated\/xxx/simulated\/10_reps\/$sim_rep\/consensus_contig_alignments/g" bin/temp_similarity_matrix_simulated_consensus.R
sed -i -e 's/xe-x/3e-4/g' bin/temp_similarity_matrix_simulated_consensus.R
sed -i -e "s/\"\/simmatrix_simulated_xxx/\"$sim_rep\/simmatrix_simulated_consensus_contig_alignments/g" bin/temp_similarity_matrix_simulated_consensus.R
Rscript bin/temp_similarity_matrix_simulated_consensus.R
rm bin/temp_*


# simulated consensus iupac sequences
cp bin/template_similarity_matrix_simulated_consensus.R bin/temp_similarity_matrix_simulated_consensus.R
sed -i -e "s/simulated\/xxx/simulated\/10_reps\/$sim_rep\/iupac_consensus_alignments/g" bin/temp_similarity_matrix_simulated_consensus.R
sed -i -e 's/xe-x/1e-5/g' bin/temp_similarity_matrix_simulated_consensus.R
sed -i -e "s/\"\/simmatrix_simulated_xxx/\"$sim_rep\/simmatrix_simulated_iupac_consensus_alignments/g" bin/temp_similarity_matrix_simulated_consensus.R
Rscript bin/temp_similarity_matrix_simulated_consensus.R
rm bin/temp_*
cp bin/template_similarity_matrix_simulated_consensus.R bin/temp_similarity_matrix_simulated_consensus.R
sed -i -e "s/simulated\/xxx/simulated\/10_reps\/$sim_rep\/iupac_consensus_alignments/g" bin/temp_similarity_matrix_simulated_consensus.R
sed -i -e 's/xe-x/1e-4/g' bin/temp_similarity_matrix_simulated_consensus.R
sed -i -e "s/\"\/simmatrix_simulated_xxx/\"$sim_rep\/simmatrix_simulated_iupac_consensus_alignments/g" bin/temp_similarity_matrix_simulated_consensus.R
Rscript bin/temp_similarity_matrix_simulated_consensus.R
rm bin/temp_*
cp bin/template_similarity_matrix_simulated_consensus.R bin/temp_similarity_matrix_simulated_consensus.R
sed -i -e "s/simulated\/xxx/simulated\/10_reps\/$sim_rep\/iupac_consensus_alignments/g" bin/temp_similarity_matrix_simulated_consensus.R
sed -i -e 's/xe-x/3e-4/g' bin/temp_similarity_matrix_simulated_consensus.R
sed -i -e "s/\"\/simmatrix_simulated_xxx/\"$sim_rep\/simmatrix_simulated_iupac_consensus_alignments/g" bin/temp_similarity_matrix_simulated_consensus.R
Rscript bin/temp_similarity_matrix_simulated_consensus.R
rm bin/temp_*

done



# _________________summarize posterior tree distribution (treeannotator)_____________________
mkdir analyses/treeannotator

# empirical allele sequences
/Applications/BEAST\ 2.4.4/bin/treeannotator -burnin 10 -heights mean analyses/stacey/empirical/allele_alignments/species.trees analyses/treeannotator/empirical_allele_alignments_stacey.tre
# empirical chimeric allele sequences
/Applications/BEAST\ 2.4.4/bin/treeannotator -burnin 10 -heights mean analyses/stacey/empirical/chimeric_allele_alignments/species.trees analyses/treeannotator/empirical_chimeric_allele_alignments_stacey.tre
# empirical consensus contig sequences
/Applications/BEAST\ 2.4.4/bin/treeannotator -burnin 10 -heights mean analyses/stacey/empirical/consensus_contig_alignments/species.trees analyses/treeannotator/empirical_consensus_contig_alignments_stacey.tre
# empirical iupac consensus sequences
/Applications/BEAST\ 2.4.4/bin/treeannotator -burnin 10 -heights mean analyses/stacey/empirical/iupac_consensus_alignments/species.trees analyses/treeannotator/empirical_iupac_consensus_alignments_stacey.tre
# simulated allele sequences
/Applications/BEAST\ 2.4.4/bin/treeannotator -burnin 10 -heights mean analyses/stacey/simulated/allele_alignments/species.trees analyses/treeannotator/simulated_allele_alignments_stacey.tre
# simulated chimeric allele sequences
/Applications/BEAST\ 2.4.4/bin/treeannotator -burnin 10 -heights mean analyses/stacey/simulated/chimeric_allele_alignments/species.trees analyses/treeannotator/simulated_chimeric_allele_alignments_stacey.tre
# simulated consensus contig sequences
/Applications/BEAST\ 2.4.4/bin/treeannotator -burnin 10 -heights mean analyses/stacey/simulated/consensus_contig_alignments/species.trees analyses/treeannotator/simulated_consensus_contig_alignments_stacey.tre
# simulated iupac consensus sequences
/Applications/BEAST\ 2.4.4/bin/treeannotator -burnin 10 -heights mean analyses/stacey/simulated/iupac_consensus_alignments/species.trees analyses/treeannotator/simulated_iupac_consensus_alignments_stacey.tre

# make trees for 10 simreps
for sim in analyses/species_da/simulated/10_reps/*;
do sim_rep=$(echo $sim | sed 's/.*\///g' | sed 's/\\r//g');
mkdir analyses/treeannotator/$sim_rep/
# simulated allele sequences
/Applications/BEAST\ 2.4.4/bin/treeannotator -burnin 10 -heights mean analyses/stacey/simulated/10_reps/$sim_rep/allele_alignments/species.trees analyses/treeannotator/$sim_rep/simulated_allele_alignments_stacey.tre
# simulated chimeric allele sequences
/Applications/BEAST\ 2.4.4/bin/treeannotator -burnin 10 -heights mean analyses/stacey/simulated/10_reps/$sim_rep/chimeric_allele_alignments/species.trees analyses/treeannotator/$sim_rep/simulated_chimeric_allele_alignments_stacey.tre
# simulated consensus contig sequences
/Applications/BEAST\ 2.4.4/bin/treeannotator -burnin 10 -heights mean analyses/stacey/simulated/10_reps/$sim_rep/consensus_contig_alignments/species.trees analyses/treeannotator/$sim_rep/simulated_consensus_contig_alignments_stacey.tre
# simulated iupac consensus sequences
/Applications/BEAST\ 2.4.4/bin/treeannotator -burnin 10 -heights mean analyses/stacey/simulated/10_reps/$sim_rep/iupac_consensus_alignments/species.trees analyses/treeannotator/$sim_rep/simulated_iupac_consensus_alignments_stacey.tre
done


# the following is necessary in order to avoid a bug in the phyloch::read.beast() funtion in R. We need to remove the very root of the tree which is printed by treeannotator 2
sed -i -e 's/:0.0;/;/g' analyses/treeannotator/empirical_allele_alignments_stacey.tre
sed -i -e 's/:0.0;/;/g' analyses/treeannotator/empirical_chimeric_allele_alignments_stacey.tre
sed -i -e 's/:0.0;/;/g' analyses/treeannotator/empirical_consensus_contig_alignments_stacey.tre
sed -i -e 's/:0.0;/;/g' analyses/treeannotator/empirical_iupac_consensus_alignments_stacey.tre
sed -i -e 's/:0.0;/;/g' analyses/treeannotator/simulated_allele_alignments_stacey.tre
sed -i -e 's/:0.0;/;/g' analyses/treeannotator/simulated_chimeric_allele_alignments_stacey.tre
sed -i -e 's/:0.0;/;/g' analyses/treeannotator/simulated_consensus_contig_alignments_stacey.tre
sed -i -e 's/:0.0;/;/g' analyses/treeannotator/simulated_iupac_consensus_alignments_stacey.tre

# also take care of some naming issues in the empirical trees
for tree in analyses/treeannotator/*.tre;
do sed -i -e 's/--/_/g' $tree;
done

for tree in analyses/treeannotator/*.tre;
do sed -i -e 's/T_//g' $tree;
done

for tree in analyses/treeannotator/*.tre;
do sed -i -e 's/pyra/T_pyra/g' $tree;
done

for tree in analyses/treeannotator/*.tre;
do sed -i -e 's/pella/T_pella/g' $tree;
done

rm analyses/treeannotator/*-e


for sim in analyses/species_da/simulated/10_reps/*;
do sim_rep=$(echo $sim | sed 's/.*\///g' | sed 's/\\r//g');
sed -i -e 's/:0.0;/;/g' analyses/treeannotator/$sim_rep/simulated_allele_alignments_stacey.tre
sed -i -e 's/:0.0;/;/g' analyses/treeannotator/$sim_rep/simulated_chimeric_allele_alignments_stacey.tre
sed -i -e 's/:0.0;/;/g' analyses/treeannotator/$sim_rep/simulated_consensus_contig_alignments_stacey.tre
sed -i -e 's/:0.0;/;/g' analyses/treeannotator/$sim_rep/simulated_iupac_consensus_alignments_stacey.tre
rm analyses/treeannotator/$sim_rep/*-e
done

# _________________________plot trees______________________________

mkdir results/trees
for sim in analyses/species_da/simulated/10_reps/*;
do sim_rep=$(echo $sim | sed 's/.*\///g' | sed 's/\\r//g');
mkdir results/trees/$sim_rep
# open and run the r-script at bin/plot_beast_trees.r




# _________________________extract and plot node height distribution______________________________
# run bin/extract_node_heights.py for every stacey analysis and then plot with bin/plot_node_heigh_distribution.r

