setwd("/Users/tobias/Desktop/figures_topaza_review")
workdir<-getwd()
workdir
dir()

ref_DE<-0.00020020000000000000
ref_YZ<-0.00025410000000000000
ref_XYZ<-0.00035330000000000000
ref_DEXYZ<-0.00083790000000000000

#get all data files ********************************************************************

#read in stacey_sim_allele files
sim_phased_node_heights_DE<-read.csv("stacey/stacey_simulated_alleles_no_outgroup/collapse_height_1e-5/node_depths/DE_node_depths.txt")
sim_phased_node_heights_YZ<-read.table("stacey/stacey_simulated_alleles_no_outgroup/collapse_height_1e-5/node_depths/YZ_node_depths.txt")
sim_phased_node_heights_XYZ<-read.csv("stacey/stacey_simulated_alleles_no_outgroup/collapse_height_1e-5/node_depths/XYZ_node_depths.txt")
sim_phased_node_heights_DEXYZ<-read.csv("stacey/stacey_simulated_alleles_no_outgroup/collapse_height_1e-5/node_depths/DEXYZ_node_depths.txt")

#read in stacey_sim_contig files
sim_contig_node_heights_DE<-read.csv("stacey/stacey_simulated_contigs_no_outgroup/collapse_height_1e-5/node_depths/DE_node_depths.txt")
sim_contig_node_heights_YZ<-read.table("stacey/stacey_simulated_contigs_no_outgroup/collapse_height_1e-5/node_depths/YZ_node_depths.txt")
sim_contig_node_heights_XYZ<-read.csv("stacey/stacey_simulated_contigs_no_outgroup/collapse_height_1e-5/node_depths/XYZ_node_depths.txt")
sim_contig_node_heights_DEXYZ<-read.csv("stacey/stacey_simulated_contigs_no_outgroup/collapse_height_1e-5/node_depths/DEXYZ_node_depths.txt")

#read in snapp_top150 files
snapp_150_node_heights_DE<-read.csv("snapp/simulated_allele_snps/top_150/node_depths/DE_node_depths.txt")
snapp_150_node_heights_YZ<-read.table("snapp/simulated_allele_snps/top_150/node_depths/YZ_node_depths.txt")
snapp_150_node_heights_XYZ<-read.csv("snapp/simulated_allele_snps/top_150/node_depths/XYZ_node_depths.txt")
snapp_150_node_heights_DEXYZ<-read.csv("snapp/simulated_allele_snps/top_150/node_depths/DEXYZ_node_depths.txt")

# read in snapp_all files
snapp_all_node_heights_DE<-read.csv("snapp/simulated_allele_snps/rep2/node_depths/DE_node_depths.txt")
snapp_all_node_heights_YZ<-read.table("snapp/simulated_allele_snps/rep2/node_depths/YZ_node_depths.txt")
snapp_all_node_heights_XYZ<-read.csv("snapp/simulated_allele_snps/rep2/node_depths/XYZ_node_depths.txt")
snapp_all_node_heights_DEXYZ<-read.csv("snapp/simulated_allele_snps/rep2/node_depths/DEXYZ_node_depths.txt")


#correct all datasets (align by root)*************************************************

# simulated phased
reference_DEXYZ <- 0.00083790000000000000
mean_DEXYZ <- mean(sim_phased_node_heights_DEXYZ[,1])
factor <- reference_DEXYZ/mean_DEXYZ
corrected_sim_phased_node_heights_DEXYZ <- sim_phased_node_heights_DEXYZ*factor
corrected_sim_phased_node_heights_XYZ <- sim_phased_node_heights_XYZ*factor
corrected_sim_phased_node_heights_YZ <- sim_phased_node_heights_YZ*factor
corrected_sim_phased_node_heights_DE <- sim_phased_node_heights_DE*factor

# simulated consensus
reference_DEXYZ <- 0.00083790000000000000
mean_DEXYZ <- mean(sim_contig_node_heights_DEXYZ[,1])
factor <- reference_DEXYZ/mean_DEXYZ
corrected_sim_contig_node_heights_DEXYZ <- sim_contig_node_heights_DEXYZ*factor
corrected_sim_contig_node_heights_XYZ <- sim_contig_node_heights_XYZ*factor
corrected_sim_contig_node_heights_YZ <- sim_contig_node_heights_YZ*factor
corrected_sim_contig_node_heights_DE <- sim_contig_node_heights_DE*factor


# snps_150
reference_DEXYZ <- 0.00083790000000000000
mean_DEXYZ <- mean(snapp_150_node_heights_DEXYZ[,1])
factor <- reference_DEXYZ/mean_DEXYZ
corrected_snapp_150_node_heights_DEXYZ <- snapp_150_node_heights_DEXYZ*factor
corrected_snapp_150_node_heights_XYZ <- snapp_150_node_heights_XYZ*factor
corrected_snapp_150_node_heights_YZ <- snapp_150_node_heights_YZ*factor
corrected_snapp_150_node_heights_DE <- snapp_150_node_heights_DE*factor

# snps_all
reference_DEXYZ <- 0.00083790000000000000
mean_DEXYZ <- mean(snapp_all_node_heights_DEXYZ[,1])
factor <- reference_DEXYZ/mean_DEXYZ
corrected_snapp_all_node_heights_DEXYZ <- snapp_all_node_heights_DEXYZ*factor
corrected_snapp_all_node_heights_XYZ <- snapp_all_node_heights_XYZ*factor
corrected_snapp_all_node_heights_YZ <- snapp_all_node_heights_YZ*factor
corrected_snapp_all_node_heights_DE <- snapp_all_node_heights_DE*factor


#___________________________________________________________________________
#corrected plots of phased vs unphased data

pdf(file="stacey/phased_unphased_comparison_corrected_DE.pdf",width = 7)

plot(density(corrected_sim_contig_node_heights_DE[,1]),xlim=c(0,0.0006),ylim=c(0,20000),xlab="Node-height",ylab="Density",main='Node-height of clade (D,E)', col = "grey")
polygon(density(corrected_sim_contig_node_heights_DE[,1]),density = 70, col = "grey")
abline(v=mean(corrected_sim_contig_node_heights_DE[,1]),lty = 5,lwd = 2, col = "grey")

points(density(corrected_sim_phased_node_heights_DE[,1]), type = 'l', col = "black")
polygon(density(corrected_sim_phased_node_heights_DE[,1]),density = 70, col = "black")
abline(v=mean(corrected_sim_phased_node_heights_DE[,1]),lty = 5,lwd = 2, col = "black")

points(density(corrected_snapp_150_node_heights_DE[,1]), type = 'l', col = "yellow3")
polygon(density(corrected_snapp_150_node_heights_DE[,1]),density = 70, col = "yellow3")
abline(v=mean(corrected_snapp_150_node_heights_DE[,1]),lty = 5,lwd = 2, col = "yellow3")

points(density(corrected_snapp_all_node_heights_DE[,1]), type = 'l', col = "orange")
polygon(density(corrected_snapp_all_node_heights_DE[,1]),density = 70, col = "orange")
abline(v=mean(corrected_snapp_all_node_heights_DE[,1]),lty = 5,lwd = 2, col = "orange")

abline(v=0.00020020000000000000,col = "green",lwd = 2)
plot_colors <- c("grey","black", "yellow3", "orange","green")
legend(4.5e-04,19000,legend = c("consensus","allele","150 SNPs","all SNPs","'true' value"),col=plot_colors, bg="white", lwd=6, cex=.9, box.lty = 0)

dev.off()



pdf(file="stacey/phased_unphased_comparison_corrected_YZ.pdf",width = 7)

plot(density(corrected_sim_contig_node_heights_YZ[,1]),xlim=c(0,0.0006),ylim=c(0,20000),xlab="Node-height",ylab="Density",main='Node-height of clade (Y,Z)', col = "grey")
polygon(density(corrected_sim_contig_node_heights_YZ[,1]),density = 70, col = "grey")
abline(v=mean(corrected_sim_contig_node_heights_YZ[,1]),lty = 5,lwd = 2, col = "grey")

points(density(corrected_sim_phased_node_heights_YZ[,1]), type = 'l', col = "black")
polygon(density(corrected_sim_phased_node_heights_YZ[,1]),density = 70, col = "black")
abline(v=mean(corrected_sim_phased_node_heights_YZ[,1]),lty = 5,lwd = 2, col = "black")

points(density(corrected_snapp_150_node_heights_YZ[,1]), type = 'l', col = "yellow3")
polygon(density(corrected_snapp_150_node_heights_YZ[,1]),density = 70, col = "yellow3")
abline(v=mean(corrected_snapp_150_node_heights_YZ[,1]),lty = 5,lwd = 2, col = "yellow3")

points(density(corrected_snapp_all_node_heights_YZ[,1]), type = 'l', col = "orange")
polygon(density(corrected_snapp_all_node_heights_YZ[,1]),density = 70, col = "orange")
abline(v=mean(corrected_snapp_all_node_heights_YZ[,1]),lty = 5,lwd = 2, col = "orange")

abline(v=0.00025410000000000000,col = "green",lwd = 2)
plot_colors <- c("grey","black", "yellow3", "orange","green")
legend(4.5e-04,19000,legend = c("consensus","allele","150 SNPs","all SNPs","'true' value"),col=plot_colors, bg="white", lwd=6, cex=.9, box.lty = 0)

dev.off()



pdf(file="stacey/phased_unphased_comparison_corrected_XYZ.pdf",width = 7)

plot(density(corrected_sim_contig_node_heights_XYZ[,1]),xlim=c(0,0.0006),ylim=c(0,20000),xlab="Node-height",ylab="Density",main='Node-height of clade (X,(Y,Z))', col = "grey")
polygon(density(corrected_sim_contig_node_heights_XYZ[,1]),density = 70, col = "grey")
abline(v=mean(corrected_sim_contig_node_heights_XYZ[,1]),lty = 5,lwd = 2, col = "grey")

points(density(corrected_sim_phased_node_heights_XYZ[,1]), type = 'l', col = "black")
polygon(density(corrected_sim_phased_node_heights_XYZ[,1]),density = 70, col = "black")
abline(v=mean(corrected_sim_phased_node_heights_XYZ[,1]),lty = 5,lwd = 2, col = "black")

points(density(corrected_snapp_150_node_heights_XYZ[,1]), type = 'l', col = "yellow3")
polygon(density(corrected_snapp_150_node_heights_XYZ[,1]),density = 70, col = "yellow3")
abline(v=mean(corrected_snapp_150_node_heights_XYZ[,1]),lty = 5,lwd = 2, col = "yellow3")

points(density(corrected_snapp_all_node_heights_XYZ[,1]), type = 'l', col = "orange")
polygon(density(corrected_snapp_all_node_heights_XYZ[,1]),density = 70, col = "orange")
abline(v=mean(corrected_snapp_all_node_heights_XYZ[,1]),lty = 5,lwd = 2, col = "orange")

abline(v=0.00035330000000000000,col = "green",lwd = 2)
plot_colors <- c("grey","black", "yellow3", "orange","green")
legend(4.5e-04,19000,legend = c("consensus","allele","150 SNPs","all SNPs","'true' value"),col=plot_colors, bg="white", lwd=6, cex=.9, box.lty = 0)

dev.off()



pdf(file="stacey/phased_unphased_comparison_corrected_DEXYZ.pdf",width = 7)

plot(density(corrected_sim_contig_node_heights_DEXYZ[,1]),xlim=c(0.0005,0.0011),ylim=c(0,20000),xlab="Node-height",ylab="Density",main='Node-height of clade ((D,E),(X,(Y,Z)))', col = "black")
polygon(density(corrected_sim_contig_node_heights_DEXYZ[,1]),density = 70, col = "grey")
abline(v=mean(corrected_sim_contig_node_heights_DEXYZ[,1]),lty = 5,lwd = 2, col = "grey")

points(density(corrected_sim_phased_node_heights_DEXYZ[,1]), type = 'l', col = "grey")
polygon(density(corrected_sim_phased_node_heights_DEXYZ[,1]),density = 70, col = "black")
abline(v=mean(corrected_sim_phased_node_heights_DEXYZ[,1]),lty = 5,lwd = 2, col = "black")

points(density(corrected_snapp_150_node_heights_DEXYZ[,1]), type = 'l', col = "yellow3")
polygon(density(corrected_snapp_150_node_heights_DEXYZ[,1]),density = 70, col = "yellow3")
abline(v=mean(corrected_snapp_150_node_heights_DEXYZ[,1]),lty = 5,lwd = 2, col = "yellow3")

points(density(corrected_snapp_all_node_heights_DEXYZ[,1]), type = 'l', col = "orange")
polygon(density(corrected_snapp_all_node_heights_DEXYZ[,1]),density = 70, col = "orange")
abline(v=mean(corrected_snapp_all_node_heights_DEXYZ[,1]),lty = 5,lwd = 2, col = "orange")

abline(v=0.00083790000000000000,col = "green",lwd = 2)
plot_colors <- c("grey","black", "yellow3", "orange","green")
legend(9.5e-04,19000,legend = c("consensus","allele","150 SNPs","all SNPs","'true' value"),col=plot_colors, bg="white", lwd=6, cex=.9, box.lty = 0)

dev.off()


#________________________________________________________________________________________
#raw plots of phased vs unphased data

pdf(file="stacey/uncorrected_phased_unphased_comparison_DE.pdf",width = 7)

plot(density(sim_contig_node_heights_DE[,1]),xlim=c(0,0.0006),ylim=c(0,20000),xlab="Node-height",ylab="Density",main='Node-height of clade (D,E)', col = "grey")
polygon(density(sim_contig_node_heights_DE[,1]),density = 70, col = "grey")
abline(v=mean(sim_contig_node_heights_DE[,1]),lty = 5,lwd = 2, col = "grey")

points(density(sim_phased_node_heights_DE[,1]), type = 'l', col = "black")
polygon(density(sim_phased_node_heights_DE[,1]),density = 70, col = "black")
abline(v=mean(sim_phased_node_heights_DE[,1]),lty = 5,lwd = 2, col = "black")

abline(v=0.00020020000000000000,col = "green")
plot_colors <- c("grey","black", "green")
legend(4.5e-04,19000,legend = c("consensus","allele","'true' value"),col=plot_colors, bg="white", lwd=6, cex=.9, box.lty = 0)

dev.off()



pdf(file="stacey/uncorrected_phased_unphased_comparison_YZ.pdf",width = 7)

plot(density(sim_contig_node_heights_YZ[,1]),xlim=c(0,0.0006),ylim=c(0,20000),xlab="Node-height",ylab="Density",main='Node-height of clade (Y,Z)', col = "grey")
polygon(density(sim_contig_node_heights_YZ[,1]),density = 70, col = "grey")
abline(v=mean(sim_contig_node_heights_YZ[,1]),lty = 5,lwd = 2, col = "grey")

points(density(sim_phased_node_heights_YZ[,1]), type = 'l', col = "black")
polygon(density(sim_phased_node_heights_YZ[,1]),density = 70, col = "black")
abline(v=mean(sim_phased_node_heights_YZ[,1]),lty = 5,lwd = 2, col = "black")

abline(v=0.00025410000000000000,col = "green")
plot_colors <- c("grey","black", "green")
legend(4.5e-04,19000,legend = c("consensus","allele","'true' value"),col=plot_colors, bg="white", lwd=6, cex=.9, box.lty = 0)

dev.off()



pdf(file="stacey/uncorrected_phased_unphased_comparison_XYZ.pdf",width = 7)

plot(density(sim_contig_node_heights_XYZ[,1]),xlim=c(0.00015,0.00075),ylim=c(0,20000),xlab="Node-height",ylab="Density",main='Node-height of clade (X,(Y,Z))', col = "grey")
polygon(density(sim_contig_node_heights_XYZ[,1]),density = 70, col = "grey")
abline(v=mean(sim_contig_node_heights_XYZ[,1]),lty = 5,lwd = 2, col = "grey")

points(density(sim_phased_node_heights_XYZ[,1]), type = 'l', col = "black")
polygon(density(sim_phased_node_heights_XYZ[,1]),density = 70, col = "black")
abline(v=mean(sim_phased_node_heights_XYZ[,1]),lty = 5,lwd = 2, col = "black")

abline(v=0.00035330000000000000,col = "green")
plot_colors <- c("grey","black", "green")
legend(6.0e-04,19000,legend = c("consensus","allele","'true' value"),col=plot_colors, bg="white", lwd=6, cex=.9, box.lty = 0)

dev.off()




pdf(file="stacey/uncorrected_phased_unphased_comparison_DEXYZ.pdf",width = 7)

plot(density(sim_contig_node_heights_DEXYZ[,1]),xlim=c(0.0006,0.0012),ylim=c(0,20000),xlab="Node-height",ylab="Density",main='Node-height of clade ((D,E),(X,(Y,Z)))', col = "grey")
polygon(density(sim_contig_node_heights_DEXYZ[,1]),density = 70, col = "grey")
abline(v=mean(sim_contig_node_heights_DEXYZ[,1]),lty = 5,lwd = 2, col = "grey")

points(density(sim_phased_node_heights_DEXYZ[,1]), type = 'l', col = "black")
polygon(density(sim_phased_node_heights_DEXYZ[,1]),density = 70, col = "black")
abline(v=mean(sim_phased_node_heights_DEXYZ[,1]),lty = 5,lwd = 2, col = "black")

abline(v=0.00083790000000000000,col = "green")
plot_colors <- c("grey","black", "green")
legend(10.5e-04,19000,legend = c("consensus","allele","'true' value"),col=plot_colors, bg="white", lwd=6, cex=.9, box.lty = 0)

dev.off()


#____________________________________________________________________________________

# separate corrected plots for each dataset
#simulated_alleles_________________
pdf(file="stacey/stacey_simulated_alleles_no_outgroup/collapse_height_1e-5/node_depths/corrected_density_plot_node_depths_sim_alleles.pdf",width = 14)
plot(density(corrected_sim_phased_node_heights_DE[,1]),xlim=c(0,0.0012), col = "green")
points(density(corrected_sim_phased_node_heights_YZ[,1]),type = 'l',col = "blue")
points(density(corrected_sim_phased_node_heights_XYZ[,1]),type = 'l', col = "red")
points(density(corrected_sim_phased_node_heights_DEXYZ[,1]),type = 'l')

mean_DE <- mean(corrected_sim_phased_node_heights_DE[,1])
mean_YZ <- mean(corrected_sim_phased_node_heights_YZ[,1])
mean_XYZ <- mean(corrected_sim_phased_node_heights_XYZ[,1])
mean_DEXYZ <- mean(corrected_sim_phased_node_heights_DEXYZ[,1])

abline(v=0.00020020000000000000,col = "green")
#abline(v=mean_DE,col = "green")
abline(v=0.00025410000000000000,col = "blue")
#abline(v=mean_YZ,col = "blue")
abline(v=0.00035330000000000000,col = "red")
#abline(v=mean_XYZ,col = "red")
abline(v=0.00083790000000000000)
#abline(v=mean_DEXYZ)

dev.off()


#simulated consensus________________
pdf(file="stacey/stacey_simulated_contigs_no_outgroup/collapse_height_1e-5/node_depths/corrected_density_plot_node_depths_sim_contigs.pdf",width = 14)
plot(density(corrected_sim_contig_node_heights_DE[,1]),xlim=c(0,0.0012), col = "green")
points(density(corrected_sim_contig_node_heights_YZ[,1]),type = 'l',col = "blue")
points(density(corrected_sim_contig_node_heights_XYZ[,1]),type = 'l', col = "red")
points(density(corrected_sim_contig_node_heights_DEXYZ[,1]),type = 'l')

mean_DE <- mean(corrected_sim_contig_node_heights_DE[,1])
mean_YZ <- mean(corrected_sim_contig_node_heights_YZ[,1])
mean_XYZ <- mean(corrected_sim_contig_node_heights_XYZ[,1])
mean_DEXYZ <- mean(corrected_sim_contig_node_heights_DEXYZ[,1])

abline(v=0.00020020000000000000,col = "green")
#abline(v=mean_DE,col = "green")
abline(v=0.00025410000000000000,col = "blue")
#abline(v=mean_YZ,col = "blue")
abline(v=0.00035330000000000000,col = "red")
#abline(v=mean_XYZ,col = "red")
abline(v=0.00083790000000000000)
#abline(v=mean_DEXYZ)

dev.off()


#snps_150______________________________
pdf(file="snapp/simulated_allele_snps/top_150/node_depths/corrected_density_plot_node_depths.pdf", width = 14)

plot(density(corrected_snapp_150_node_heights_DE[,1]),xlim=c(0,0.0012), col = "green")
points(density(corrected_snapp_150_node_heights_YZ[,1]),type = 'l',col = "blue")
points(density(corrected_snapp_150_node_heights_XYZ[,1]),type = 'l', col = "red")
points(density(corrected_snapp_150_node_heights_DEXYZ[,1]),type = 'l')

mean_DE <- mean(corrected_snapp_150_node_heights_DE[,1])
mean_YZ <- mean(corrected_snapp_150_node_heights_YZ[,1])
mean_XYZ <- mean(corrected_snapp_150_node_heights_XYZ[,1])
mean_DEXYZ <- mean(corrected_snapp_150_node_heights_DEXYZ[,1])

abline(v=0.00020020000000000000,col = "green")
abline(v=0.00025410000000000000,col = "blue")
abline(v=0.00035330000000000000,col = "red")
abline(v=0.00083790000000000000)

dev.off()


#snps_all____________________________________
pdf(file="snapp/simulated_allele_snps/rep2/node_depths/corrected_density_plot_node_depths.pdf", width = 14)

plot(density(corrected_snapp_all_node_heights_DE[,1]),xlim=c(0,0.0012), col = "green")
points(density(corrected_snapp_all_node_heights_YZ[,1]),type = 'l',col = "blue")
points(density(corrected_snapp_all_node_heights_XYZ[,1]),type = 'l', col = "red")
points(density(corrected_snapp_all_node_heights_DEXYZ[,1]),type = 'l')

mean_DE <- mean(corrected_snapp_all_node_heights_DE[,1])
mean_YZ <- mean(corrected_snapp_all_node_heights_YZ[,1])
mean_XYZ <- mean(corrected_snapp_all_node_heights_XYZ[,1])
mean_DEXYZ <- mean(corrected_snapp_all_node_heights_DEXYZ[,1])

abline(v=0.00020020000000000000,col = "green")
abline(v=0.00025410000000000000,col = "blue")
abline(v=0.00035330000000000000,col = "red")
abline(v=0.00083790000000000000)

dev.off()

#________________________________________________________________________________________

#separate raw plots for each dataset:

#sim_alleles____________________________________
pdf(file="stacey/stacey_simulated_alleles_no_outgroup/collapse_height_1e-5/node_depths/density_plot_node_depths.pdf",width = 14)

plot(density(sim_phased_node_heights_DE[,1]),xlim=c(0,0.0012), col = "green")
points(density(sim_phased_node_heights_YZ[,1]),type = 'l',col = "blue")
points(density(sim_phased_node_heights_XYZ[,1]),type = 'l', col = "red")
points(density(sim_phased_node_heights_DEXYZ[,1]),type = 'l')

mean_DE <- mean(sim_phased_node_heights_DE[,1])
mean_YZ <- mean(sim_phased_node_heights_YZ[,1])
mean_XYZ <- mean(sim_phased_node_heights_XYZ[,1])
mean_DEXYZ <- mean(sim_phased_node_heights_DEXYZ[,1])

abline(v=0.00020020000000000000,col = "green")
#abline(v=mean_DE,col = "green")
abline(v=0.00025410000000000000,col = "blue")
#abline(v=mean_YZ,col = "blue")
abline(v=0.00035330000000000000,col = "red")
#abline(v=mean_XYZ,col = "red")
abline(v=0.00083790000000000000)
#abline(v=mean_DEXYZ)

dev.off()


#sim_contigs____________________________
pdf(file="stacey/stacey_simulated_contigs_no_outgroup/collapse_height_1e-5/node_depths/density_plot_node_depths.pdf", width = 14)

plot(density(sim_contig_node_heights_DE[,1]),xlim=c(0,0.0012), col = "green")
points(density(sim_contig_node_heights_YZ[,1]),type = 'l',col = "blue")
points(density(sim_contig_node_heights_XYZ[,1]),type = 'l', col = "red")
points(density(sim_contig_node_heights_DEXYZ[,1]),type = 'l')

mean_DE <- mean(sim_contig_node_heights_DE[,1])
mean_YZ <- mean(sim_contig_node_heights_YZ[,1])
mean_XYZ <- mean(sim_contig_node_heights_XYZ[,1])
mean_DEXYZ <- mean(sim_contig_node_heights_DEXYZ[,1])

abline(v=0.00020020000000000000,col = "green")
abline(v=0.00025410000000000000,col = "blue")
abline(v=0.00035330000000000000,col = "red")
abline(v=0.00083790000000000000)

dev.off()


#snps_150______________________________
pdf(file="snapp/simulated_allele_snps/top_150/node_depths/density_plot_node_depths.pdf", width = 14)

plot(density(snapp_150_node_heights_DE[,1]),xlim=c(0,0.5), col = "green")
points(density(snapp_150_node_heights_YZ[,1]),type = 'l',col = "blue")
points(density(snapp_150_node_heights_XYZ[,1]),type = 'l', col = "red")
points(density(snapp_150_node_heights_DEXYZ[,1]),type = 'l')

mean_DE <- mean(snapp_150_node_heights_DE[,1])
mean_YZ <- mean(snapp_150_node_heights_YZ[,1])
mean_XYZ <- mean(snapp_150_node_heights_XYZ[,1])
mean_DEXYZ <- mean(snapp_150_node_heights_DEXYZ[,1])

abline(v=0.00020020000000000000,col = "green")
abline(v=0.00025410000000000000,col = "blue")
abline(v=0.00035330000000000000,col = "red")
abline(v=0.00083790000000000000)

dev.off()


#snps_all____________________________________
pdf(file="snapp/simulated_allele_snps/rep2/node_depths/density_plot_node_depths.pdf", width = 14)

plot(density(snapp_all_node_heights_DE[,1]),xlim=c(0,0.06), col = "green")
points(density(snapp_all_node_heights_YZ[,1]),type = 'l',col = "blue")
points(density(snapp_all_node_heights_XYZ[,1]),type = 'l', col = "red")
points(density(snapp_all_node_heights_DEXYZ[,1]),type = 'l')

mean_DE <- mean(snapp_all_node_heights_DE[,1])
mean_YZ <- mean(snapp_all_node_heights_YZ[,1])
mean_XYZ <- mean(snapp_all_node_heights_XYZ[,1])
mean_DEXYZ <- mean(snapp_all_node_heights_DEXYZ[,1])

abline(v=0.00020020000000000000,col = "green")
abline(v=0.00025410000000000000,col = "blue")
abline(v=0.00035330000000000000,col = "red")
abline(v=0.00083790000000000000)

dev.off()
