setwd("/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/analyses/species_da/empirical/xxx")
workdir<-getwd()
workdir
dir()
x<-read.table("species_da_results_xe-x.txt", sep="", header=TRUE)
x
x$Florisuga <- NULL
y<-names(x)
y
renames <- matrix(c(
"pyra3", "T_pyra3",
"pella5", "T_pella5",
"pyra2", "T_pyra2",
"pella6", "T_pella6",
"pyra1", "T_pyra1",
"pella7", "T_pella7",
"pella8", "T_pella8",
"pella9", "T_pella9",
"pyra4", "T_pyra4"),
nrow=9, ncol=2, byrow=TRUE)
renames
# define the columns that should be pursued (remove the 'count' 'fraction' 'similarity' and 'nclusters' column)
mincl.names<-colnames(x)[-(1:4)]
for (i in 1:length(mincl.names)) {
  stopifnot(mincl.names[i] == renames[i,1])
}
mincl.names[1]
renames[1,1]
#make similarity matrix
displaynames <- renames[,2]
nmincls <- length(displaynames)
sim <- matrix(0, ncol=nmincls, nrow=nmincls, dimnames=list(displaynames, displaynames))
for (i in 1:nmincls) {
  for (j in 1:nmincls) {
    coli <- x[,mincl.names[i]]
    colj <- x[,mincl.names[j]]
    w <- coli == colj
    sim[i,j] <- sum(x[w,"fraction"])
  }
}
sim <- pmin(sim,1)
neworder <- c(5,3,1,9, 2,4,6,7,8)
dividers<-c(0,4,9)
plot.rectangle <- function(v1,v2,...)
{
polygon(c(v1[1],v2[1],v2[1],v1[1]), c(v1[2],v1[2],v2[2],v2[2]), ...)
}
plot.simmatrix <- function() {
  par(mar= c(0,5,5,0)+.1)
  plot(NULL, xlim=c(0,nmincls), ylim=c(nmincls,0), axes=FALSE, ylab="", xlab="")
  axis(3, at=(1:nmincls)-.5, displaynames[neworder], tick=FALSE, las=2, line=-1)
  axis(2, at=(1:nmincls)-.5, displaynames[neworder], tick=FALSE, las=2, line=-1)
  for (i in 1:nmincls) {
    for (j in 1:nmincls) {
      d <- 1 - sim[neworder[i],neworder[j]]
      plot.rectangle(c(i-1,j-1), c(i,j), col=rgb(d,d,d), border="white")
    }
  }
  for (b in dividers) {
    lines(x=c(-.5,nmincls), y=c(b,b))
    lines(x=c(b,b), y=c(-.5,nmincls))
  }
  #legend(nmincls-4,0,legend = c("0%","25%","50%","75%","100%"),col=c(rgb(1,1,1),rgb(.75,.75,.75),rgb(.5,.5,.5),rgb(.25,.25,.25),rgb(0,0,0)),bg="white", lwd=15, cex=2, box.lty = 1)
}
print(sim[neworder,neworder], digits=2)
#plot.simmatrix()
pdf(file=paste('/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/results/simmatrix_plots/', "simmatrix_empirical_xxx_stacey_xe-x.pdf", sep=""))
plot.simmatrix()
dev.off()

