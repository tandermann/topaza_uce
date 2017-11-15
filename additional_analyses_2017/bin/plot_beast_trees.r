#install.packages("/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/bin/phyloch_1.5-3.tgz",repos=NULL)
#install.packages("XML",repos="http://cran.us.r-project.org")

# install the required packages in the correct folder
#install.packages("phytools",dependencies = TRUE,lib="/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/bin/r_packages/")
#install.packages("/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/bin/phyloch_1.5-3.tgz",repos=NULL,lib="/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/bin/r_packages/")
#install.packages("ape",dependencies = TRUE,lib="/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/bin/r_packages/")
#install.packages("Hmisc",dependencies = TRUE,lib="/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/bin/r_packages/")

#library(XML,lib.loc="/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/bin/r_packages/")
#library(Formula,lib.loc="/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/bin/r_packages/")
#library(acepack,lib.loc="/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/bin/r_packages/")
#library(gridExtra,lib.loc="/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/bin/r_packages/")
#library(checkmate,lib.loc="/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/bin/r_packages/")
#library(htmlTable,lib.loc="/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/bin/r_packages/")


library(ape,lib.loc="/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/bin/r_packages/")
library(phyloch,lib.loc="/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/bin/r_packages/")
library(Hmisc,lib.loc="/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/bin/r_packages/")
library(phytools,lib.loc="/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/bin/r_packages/")


setwd("/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/")

# read the trees
tree1 <- phyloch::read.beast(file = "analyses/treeannotator/empirical_allele_alignments_stacey.tre")
tree2 <- phyloch::read.beast(file = "analyses/treeannotator/empirical_chimeric_allele_alignments_stacey.tre")
tree3 <- phyloch::read.beast(file = "analyses/treeannotator/empirical_consensus_contig_alignments_stacey.tre")
tree4 <- phyloch::read.beast(file = "analyses/treeannotator/empirical_iupac_consensus_alignments_stacey.tre")
tree5 <- phyloch::read.beast(file = "analyses/treeannotator/rep10/simulated_allele_alignments_stacey.tre")
tree6 <- phyloch::read.beast(file = "analyses/treeannotator/rep10/simulated_chimeric_allele_alignments_stacey.tre")
tree7 <- phyloch::read.beast(file = "analyses/treeannotator/rep10/simulated_consensus_contig_alignments_stacey.tre")
tree8 <- phyloch::read.beast(file = "analyses/treeannotator/rep10/simulated_iupac_consensus_alignments_stacey.tre")

# captialize tip labels in simulated trees
tree5$tip.label <- Hmisc::capitalize(tree5$tip.label)
tree6$tip.label <- Hmisc::capitalize(tree6$tip.label)
tree7$tip.label <- Hmisc::capitalize(tree7$tip.label)
tree8$tip.label <- Hmisc::capitalize(tree8$tip.label)

# check the node numbers and see where we want to turn the tree nodes
phytools::plotTree(tree5,pts=F,node.numbers=T)
tree1<-phytools::rotateNodes(tree1,c(29,30,33,32,31,35,34,20,21,22,23,24,26,25,27,28))
tree2<-phytools::rotateNodes(tree2,c(29,30,33,32,31,35,34,20,21,22,23,24,26,25,27,28))
tree3<-phytools::rotateNodes(tree3,c(15,16,17,11,12,13,14))
tree4<-phytools::rotateNodes(tree4,c(15,16,17,11,12,13,14))
tree6<-phytools::rotateNodes(tree6,c(19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35))
tree5<-rotateConstr(tree5, tree6$tip.label)
tree7<-phytools::rotateNodes(tree7,c(10,11,12,13,14,15,16,17))
tree8<-phytools::rotateNodes(tree8,c(10,11,12,13,14,15,16,17))

# define the pink that I want to use for the northern T pyra clade (to match coloration in mitochondrial tree)
pink_cust = rgb(255,0,238,maxColorValue = 255)
# define the order of the bracnh colorations for allele and consensus trees
color_allele <- c("black","black","blue","blue","blue","blue","blue","blue","blue","black","black","black","red","red","red","red","red","red","red","black","green","green","green","green","green","green","green",pink_cust,pink_cust,pink_cust,pink_cust,pink_cust,pink_cust,pink_cust)
color_consensus <- c("black","black","blue","blue","blue","black","red","red","red","black","green","green","green",pink_cust,pink_cust,pink_cust)
color_consensus_2 <- c("black","black","blue","blue","blue","black","red","red","red","black","black","green","green","green",pink_cust,pink_cust)

#pdf("results/trees/combined.pdf")
#par(mfrow=c(2,4))

pdf("results/trees/empirical_allele_alignments_stacey.pdf")
plot(tree1, edge.color = color_allele, edge.width = 2,x.lim = c(-0.00026, 0.00196))
add.scale.bar(length = 1e-04,lwd = 2)
node.support(tree1$posterior,cutoff=0,digits=2,cex=1,mode="numbers",pos = "above")
HPDbars(tree1, label = "height_95%_HPD",col = "skyblue",lwd=4)
dev.off()

pdf("results/trees/empirical_chimeric_allele_alignments_stacey.pdf")
plot(tree2, edge.color = color_allele, edge.width = 2,x.lim= c(-3.1e-04, 0.00237))
add.scale.bar(length = 1e-04,lwd = 2)
node.support(tree2$posterior,cutoff=0,digits=2,cex=1,mode="numbers",pos = "above")
HPDbars(tree2, label = "height_95%_HPD",col = "skyblue",lwd=4)
dev.off()

pdf("results/trees/empirical_consensus_contig_alignments_stacey.pdf")
plot(tree3, edge.color = color_consensus_2, edge.width = 2,x.lim=c(-0.00029, 0.002))
add.scale.bar(length = 1e-04,lwd = 2)
node.support(tree3$posterior,cutoff=0,digits=2,cex=1,mode="numbers",pos = "above")
HPDbars(tree3, label = "height_95%_HPD",col = "skyblue",lwd=4)
dev.off()

pdf("results/trees/empirical_iupac_consensus_alignments_stacey.pdf")
plot(tree4, edge.color = color_consensus, edge.width = 2,x.lim=c(-0.00024, 0.00119))
add.scale.bar(length = 1e-04,lwd = 2)
node.support(tree4$posterior,cutoff=0,digits=2,cex=1,mode="numbers",pos = "above")
HPDbars(tree4, label = "height_95%_HPD",col = "skyblue",lwd=4)
dev.off()

pdf("results/trees/rep10/simulated_allele_alignments_stacey.pdf")
plot(tree5, edge.color = color_allele, edge.width = 2,x.lim=c(-0.00017, 0.00106))
add.scale.bar(length = 1e-04,lwd = 2)
node.support(tree5$posterior,cutoff=0,digits=2,cex=1,mode="numbers",pos = "above")
HPDbars(tree5, label = "height_95%_HPD",col = "skyblue",lwd=4)
dev.off()

pdf("results/trees/rep10/simulated_chimeric_allele_alignments_stacey.pdf")
plot(tree6, edge.color = color_allele, edge.width = 2,x.lim=c(-0.00020, 0.00136))
add.scale.bar(length = 1e-04,lwd = 2)
node.support(tree6$posterior,cutoff=0,digits=2,cex=1,mode="numbers",pos = "above")
HPDbars(tree6, label = "height_95%_HPD",col = "skyblue",lwd=4)
dev.off()

pdf("results/trees/rep10/simulated_consensus_contig_alignments_stacey.pdf")
plot(tree7, edge.color = color_consensus, edge.width = 2,x.lim=c(-2.1e-04, 0.00134))
add.scale.bar(length = 1e-04,lwd = 2)
node.support(tree7$posterior,cutoff=0,digits=2,cex=1,mode="numbers",pos = "above")
HPDbars(tree7, label = "height_95%_HPD",col = "skyblue",lwd=4)
dev.off()

pdf("results/trees/rep10/simulated_iupac_consensus_alignments_stacey.pdf")
plot(tree8, edge.color = color_consensus, edge.width = 2,x.lim=c(-7.5e-05, 0.00028))
add.scale.bar(length = 1e-04,lwd = 2)
node.support(tree8$posterior,cutoff=0,digits=2,cex=1,mode="numbers",pos = "above")
HPDbars(tree8, label = "height_95%_HPD",col = "skyblue",lwd=4)
dev.off()


# this plots the reference tree
#ref_tree <- read.nexus(file = "simulation_input_tree.tre")
#
#new_tree<-ref_tree[[2]]
#new_tree<-rotateNodes(new_tree, c(6,7,9))
#tre$edge.color <- c("black","black","blue","black","red","black","green",rgb(255,0,238,maxColorValue = 255))
#
#pdf("simulation_input_tree.pdf")
#plot(new_tree, edge.color = tre$edge.color, edge.width = 2)
#add.scale.bar(length = 1e-04,lwd = 2)
#dev.off()


