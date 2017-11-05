setwd("/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/analyses/species_da/simulated/xxx")
workdir<-getwd()
workdir
dir()
x<-read.table("species_da_results_xe-x.txt", sep="", header=TRUE)
x
x$Florisuga <- NULL
y<-names(x)
y
renames <- matrix(c(
"d1", "D1",
"d2", "D2",
"e1", "E1",
"e2", "E2",
"x1", "X1",
"x2", "X2",
"y1", "Y1",
"z1", "Z1",
"z2", "Z2"),
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
neworder <- c(1,2,3,4, 5,6,7,8,9)
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
pdf(file=paste('/Users/tobias/GitHub/topaza_uce/additional_analyses_2017/results/simmatrix_plots/', "/simmatrix_simulated_xxx_stacey_xe-x.pdf", sep=""))
plot.simmatrix()
dev.off()

