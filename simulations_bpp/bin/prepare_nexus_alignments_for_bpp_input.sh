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
sed -i '' 's/suga0/suga0^1/g' $phy;
sed -i '' 's/suga1/suga1^2/g' $phy;
sed -i '' 's/10 /10^3 /g' $phy;
sed -i '' 's/11 /11^4 /g' $phy;
sed -i '' 's/20 /20^5 /g' $phy;
sed -i '' 's/21 /21^6 /g' $phy;
sed -i '' 's/30 /30^7 /g' $phy;
sed -i '' 's/31 /31^8 /g' $phy;
sed -i '' 's/40 /40^9 /g' $phy;
sed -i '' 's/41 /41^10 /g' $phy;
sed -i '' 's/50 /50^11 /g' $phy;
sed -i '' 's/51 /51^12 /g' $phy;
sed -i '' 's/60 /60^13 /g' $phy;
sed -i '' 's/61 /61^14 /g' $phy;
sed -i '' 's/70 /70^15 /g' $phy;
sed -i '' 's/71 /71^16 /g' $phy;
sed -i '' 's/80 /80^17 /g' $phy;
sed -i '' 's/81 /81^18 /g' $phy;
sed -i '' 's/90 /90^19 /g' $phy;
sed -i '' 's/91 /91^20 /g' $phy;
#this is for the header of the resulting phylip file, which accidentally also received a '^5 at the end, which needs to be removed'
sed -i '' 's/ 20^5 / 20 /g' $phy;

split -l 21 $phy;
paste xa* > $phy.bpp;
perl -p -i -e 's/\t//g' $phy.bpp;
perl -p -i -e 's/          //g' $phy.bpp;
rm xa*;
mv $phy.bpp $alignments-edited/bpp-formatted

done;