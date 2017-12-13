for folder in *;
do sed -i -e 's/<branchRateModel .*\/>//g' $folder/*.xml
sed -i -e 's/spec="TreeLikelihood"/spec="TreeLikelihood" branchRateModel="@StrictClock.c:sim_allele_alignment_2"/g' $folder/*.xml
sed -i -e 's/<distribution id="treeLikelihood.sim_allele_alignment_2" spec="TreeLikelihood" branchRateModel="@StrictClock.c:sim_allele_alignment_2"/<distribution id="treeLikelihood.sim_allele_alignment_2" spec="TreeLikelihood"/g' $folder/*.xml
#sed -i -e 's/branchRateModel="@StrictClock.c:sim_allele_alignment_2" branchRateModel="@StrictClock.c:sim_allele_alignment_2"/branchRateModel="@StrictClock.c:sim_allele_alignment_2"/g' $folder/*.xml
rm $folder/*-e
done

cd 10_reps

for rep in *;
do cd $rep;
for folder in *;
do sed -i -e 's/<branchRateModel .*\/>//g' $folder/*.xml
sed -i -e 's/spec="TreeLikelihood"/spec="TreeLikelihood" branchRateModel="@StrictClock.c:sim_allele_alignment_2"/g' $folder/*.xml
sed -i -e 's/id="treeLikelihood.sim_consensus_alignment_2" spec="TreeLikelihood" branchRateModel="@StrictClock.c:sim_allele_alignment_2"/id="treeLikelihood.sim_consensus_alignment_2" spec="TreeLikelihood"/g' $folder/*.xml
#sed -i -e 's/branchRateModel="@StrictClock.c:sim_allele_alignment_2" branchRateModel="@StrictClock.c:sim_allele_alignment_2"/branchRateModel="@StrictClock.c:sim_allele_alignment_2"/g' $folder/*.xml
rm $folder/*-e
done
cd ..
done