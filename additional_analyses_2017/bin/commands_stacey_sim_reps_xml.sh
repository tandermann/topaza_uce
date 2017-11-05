cd ../analyses/stacey/simulated
mkdir 10_reps
cd 10_reps
mkdir rep1 rep2 rep3 rep4 rep5 rep6 rep7 rep8 rep9 rep10
for rep in rep*; do cd $rep; mkdir allele_alignments chimeric_allele_alignments consensus_contig_alignments iupac_consensus_alignments; cd ..; done
for rep in rep*; do cd $rep; cd allele_alignments; cp ../../../allele_alignments/*.xml .; cd ..;cd chimeric_allele_alignments;cp ../../../chimeric_allele_alignments/template.xml .; cd ..; cd consensus_contig_alignments;cp ../../../consensus_contig_alignments/*.xml .; cd ..; cd iupac_consensus_alignments;cp ../../../iupac_consensus_alignments/template.xml .; cd ..; cd ..; done
for rep in rep*; do cd $rep; cd allele_alignments;rm modified*.xml; cd ..;cd chimeric_allele_alignments;rm modified*.xml; cd ..; cd consensus_contig_alignments;rm modified*.xml; cd ..; cd iupac_consensus_alignments;rm modified*.xml; cd ..; cd ..; done
for rep in rep*; do cd $rep; cd allele_alignments;mv simulated*.xml template.xml; cd ..; cd consensus_contig_alignments;mv simulated*.xml template.xml; cd ..; cd ..; done

for rep in rep*;
do cd $rep;
cd allele_alignments;
python ../../../../../../bin/read_edit_xml_sim_master.py --alignments ../../../../../../data/simulated/10_rep/selection/allele_alignments_simulated/$rep/ --xml_path .;
cd ..;
cd chimeric_allele_alignments;
python ../../../../../../bin/read_edit_xml_sim_master.py --alignments ../../../../../../data/simulated/10_rep/selection/chimeric_allele_alignments_simulated/$rep/ --xml_path .;
cd ..;
cd consensus_contig_alignments;
python ../../../../../../bin/read_edit_xml_sim_master.py --alignments ../../../../../../data/simulated/10_rep/selection/contig_consensus_alignments_simulated/$rep/ --xml_path .;
cd ..;
cd iupac_consensus_alignments;
python ../../../../../../bin/read_edit_xml_sim_master.py --alignments ../../../../../../data/simulated/10_rep/selection/iupac_consensus_alignments_simulated/$rep/ --xml_path .;
cd ..;
cd ..;
done


