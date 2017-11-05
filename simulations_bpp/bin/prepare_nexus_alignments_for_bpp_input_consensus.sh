#this scripts needs a folder name as argument, that contains alignments in nexus format. The alignments will be converted into a bpp specific phylip format
alignments=$(basename $1);
echo processing the folder $alignments;
mkdir $alignments-edited;
mkdir $alignments-edited/modified_nex;
for nex in $1/*.nexus;
do cp $nex $alignments-edited/modified_nex;
done;
#cp -r $1 ./$alignments-backup;
for nex1 in $alignments-edited/modified_nex/*.nexus;
do sed -i '' 's/Flori//g' $nex1;
sed -i '' 's/_//g' $nex1;
done;
phyluce_align_convert_one_align_to_another --alignments $alignments-edited/modified_nex --input-format nexus --output $alignments-edited/phylip --output-format phylip;
mkdir $alignments-edited/bpp-formatted
for phy in $alignments-edited/phylip/*.phylip;
do echo processing $phy;
sed -i '' 's/suga/Florisuga/' $phy;

sed -i '' 's/9 /9^10 /g' $phy;
sed -i '' 's/8 /8^9 /g' $phy;
sed -i '' 's/7 /7^8 /g' $phy;
sed -i '' 's/6 /6^7 /g' $phy;
sed -i '' 's/5 /5^6 /g' $phy;
sed -i '' 's/4 /4^5 /g' $phy;
sed -i '' 's/3 /3^4 /g' $phy;
sed -i '' 's/2 /2^3 /g' $phy;
sed -i '' 's/1 /1^2 /g' $phy;
sed -i '' 's/suga/suga^1/g' $phy;

split -l 11 $phy;
paste xa* > $phy.bpp;
perl -p -i -e 's/\t//g' $phy.bpp;
perl -p -i -e 's/          //g' $phy.bpp;
rm xa*;
mv $phy.bpp $alignments-edited/bpp-formatted

done;