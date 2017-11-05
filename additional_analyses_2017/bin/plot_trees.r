library(ape)
library(diversitree)
library(ggtree)
library(phytools)
source("/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/bin/color.terminal.branches.R")


setwd("/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/analyses/treeannotator/")
dir()




tre <- read.nexus("empirical_allele_alignments_stacey.tre")

new_nodelabel <- list()
i <- 1
for (bootstrap in tre$node.label) {
  value <- round(as.numeric(bootstrap),digits=2)
  new_nodelabel[i] <- value
  i <- i+1
}
tre$node.label <- rapply(new_nodelabel, c)

#pdf("mpest_sim_allele_10_reps.pdf")
plot(tre,main='j) MPEST species tree, simulated allele data, rep10',type = "phylogram",use.edge.length=FALSE,show.node.label = TRUE, edge.width = 3, cex=1.5, font=1)
add.scale.bar(0.1,10)
#dev.off()






#this function was used for coloring the branches for the ms figures

pdf("mpest/clade_assignments/all_joined_simulated_bootstrap_mpest_clade_assignments_no_edge_length.pdf")
par(mfrow=c(2,2))

tre$edge.color <- c("black","black","black","green",rgb(255,0,238,maxColorValue = 255),"black","red","black","blue","black")
tre2$edge.color <- c("black","black","green","black","red","black","black","black",rgb(255,0,238,maxColorValue = 255),"blue")

plot(tre2,type = "phylogram",use.edge.length=TRUE,show.node.label = TRUE,edge.color = tre2$edge.color, edge.width = 3, cex=1.5, font=1, main = "Simulated consensus sequences")

plot(tre,type = "phylogram",use.edge.length=TRUE,show.node.label = TRUE,edge.color = tre$edge.color, edge.width = 3, cex=1.5, font=1, main = "Simulated allele sequences")

plot(tre2,type = "phylogram",use.edge.length=FALSE,show.node.label = TRUE,edge.color = tre2$edge.color, edge.width = 3, cex=1.5, font=1, main = "consensus, no branch lengths")

plot(tre,type = "phylogram",use.edge.length=FALSE,show.node.label = TRUE,edge.color = tre$edge.color, edge.width = 3, cex=1.5, font=1, main = "allele, no branch lengths")

dev.off()

#plotTree(tre,pts=F,node.numbers=T)




