require(phytools)
require(ape)
library(phyloch)
library(diversitree)
#library(ggtree)
library(phangorn)
library(Hmisc)
?read.beast()
setwd("/Users/tobias/Desktop/figures_topaza_review/snapp/simulated_allele_snps/rep2/")
tree1 <- phyloch::read.beast(file = "topaza_one_snp_per_locus_snapp_raw.tre")
#tree2 <- read.beast(file = "stacey_sim_allele_collapse_height_1e5.tre")
#tree3 <- phyloch::read.beast(file = "sim_allele_snps_rep2.tre")
#tree4 <- read.beast(file = "sim_top_150_allele_snps_snapp.tre")
plot(tree1)
tree1$tip.label <- capitalize(tree1$tip.label)
plotTree(tree1,pts=F,node.numbers=T)
tree1<-rotateNodes(tree1, c(15,16,17,11,12,13,14))

color <- c("black","black","blue","blue","blue","black","red","red","red","black","green","green","green",rgb(255,0,238,maxColorValue = 255),rgb(255,0,238,maxColorValue = 255),rgb(255,0,238,maxColorValue = 255))

pdf("snapp_snps_topaza.pdf")
plot(tree1, edge.color = color, edge.width = 2, x.lim=c(-0.00853, 0.04509))
add.scale.bar(length = 1e-03,lwd = 2)
node.support(tree1$posterior,cutoff=0,digits=2,cex=1,mode="numbers",pos = "above")
HPDbars(tree1, label = "height_95%_HPD",col = "skyblue",lwd=4)
dev.off()




# this plots the reference tree
tree5 <- read.nexus(file = "simulation_input_tree.tre")

new_tree<-tree5[[2]]
new_tree<-rotateNodes(new_tree, c(6,7,9))
tre$edge.color <- c("black","black","blue","black","red","black","green",rgb(255,0,238,maxColorValue = 255))


pdf("simulation_input_tree.pdf")
plot(new_tree, edge.color = tre$edge.color, edge.width = 2)
add.scale.bar(length = 1e-04,lwd = 2)
dev.off()

