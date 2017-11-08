# I ran the following on albiorix, because both makeblastdb and last are installed there
# get the probe file
wget https://raw.githubusercontent.com/faircloth-lab/uce-probe-sets/master/uce-2.5k-probe-set/uce-2.5k-probes.fasta
# make blast database
makeblastdb -in uce-2.5k-probes.fasta -parse_seqids -dbtype nucl
# make copy of contig file in which we will mask ambiguities as 'N' (required by lastz)
cp abyss_contigs_T_pella5.fa abyss_contigs_T_pella5_no_ambiguities.fa

######### WITHOUT AMBIGUITIES
# replace all IUPAC ambiguity codes with 'N'
sed -i -e 's/[RYKMSWBDHVrykmswbdhv]/N/g' abyss_contigs_T_pella5_no_ambiguities.fa
# run the lastz blast algorithm
lastz abyss_contigs_T_pella5_no_ambiguities.fa[multiple,nameparse=full] uce-2.5k-probes.fasta[nameparse=full] --strand=both --seed=12of19 --transition --nogfextend --nochain --gap=400,30 --xdrop=910 --ydrop=8370 --hspthresh=3000 --gappedthresh=3000 --noentropy --coverage=80 --identity=80 --format=general:score,name1,strand1,zstart1,end1,length1,name2,strand2,zstart2,end2,length2,diff,cigar,identity,continuity > lastz_results.txt

######### WITH AMBIGUITIES
# run lastz with ambiguity flag
lastz abyss_contigs_T_pella5.fa[multiple,nameparse=full] uce-2.5k-probes.fasta[nameparse=full] --strand=both --seed=12of19 --transition --nogfextend --nochain --gap=400,30 --xdrop=910 --ydrop=8370 --hspthresh=3000 --gappedthresh=3000 --noentropy --coverage=80 --identity=80 --ambiguous=iupac --format=general:score,name1,strand1,zstart1,end1,length1,name2,strand2,zstart2,end2,length2,diff,cigar,identity,continuity > lastz_results_iupac.txt

# run python script to extract target contig sequences and rename headers accordingly
python ./read_lastz_output.py
