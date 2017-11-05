setwd("/Users/tobias/GitHub/topaza_uce/summary_stats/snps/")
workdir<-getwd()
workdir
dir()

#read in input file
dist_count<-read.table("distance_to_center_count.txt")
dist_count<-data.frame(dist_count)


#make sure that the first column has no missing values, and if so, insert them and add value 0 in the second column
new_dist_count <- dist_count
#make a list with all whole numbers between 0 and 441(the longest dist in dataframe)
full_values <- c(0:441)
for (i in full_values) {
  if (i %in% dist_count[,1]) {
    print(i)
  }
  else {
    newrow = c(i,0)
    new_dist_count = rbind(new_dist_count,newrow)
    print (newrow)
  }
}
#sort the new dict that contains all missing values by its first column
new_dist_count <- new_dist_count[order(new_dist_count$V1),]

#call the two different columns 'dist' and 'count'
dist<-new_dist_count[,1]
count<-new_dist_count[,2]
count
barplot(count)

#group the dist column into equal bins (each size 10), to visualize the data better
dist_new<-cut(dist,breaks=c(seq(0, 440, 10)))
#sum up the counts within each bin
sums<-aggregate(count~dist_new,FUN=sum)
#transform the rather ugly values in the first column in a better format for the final plot
for(i in sums[,1]) {
  label<-toString(i)
  label_new<-strsplit(label,',')
  label_a<-as.numeric(gsub('\\(', '', label_new[[1]][1]))
  label_b<-as.numeric(gsub('\\]', '', label_new[[1]][2]))
  final_label<-(label_b-5)
  sums<-lapply(sums, gsub, pattern = i, replacement = final_label, fixed = TRUE)
}
sums
pdf('barplot_snp_distance_to_center.pdf')
barplot(as.numeric(sums$count),names.arg=sums$dist_new,xlab = "Distance from center",ylab = "SNPs extracted from position x across all alignments", main= "Positions of extracted SNPs in alignments")
dev.off()
?barplot
table


