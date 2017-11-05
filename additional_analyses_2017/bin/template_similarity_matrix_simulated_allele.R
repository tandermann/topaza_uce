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
"e1..0", "E1_0",
"d1..0", "D1_0",
"z1..0", "Z1_0",
"y1..0", "Y1_0",
"x1..0", "X1_0",
"d2..0", "D2_0",
"e2..0", "E2_0",
"x2..0", "X2_0",
"z2..0", "Z2_0",
"e1..1", "E1_1",
"d1..1", "D1_1",
"z1..1", "Z1_1",
"y1..1", "Y1_1",
"x1..1", "X1_1",
"d2..1", "D2_1",
"e2..1", "E2_1",
"x2..1", "X2_1",
"z2..1", "Z2_1"),
nrow=18, ncol=2, byrow=TRUE)
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
neworder <- c(2,11,6,15,1,10,7,16, 5,14,8,17,4,13,3,12,9,18)
dividers<-c(0,8,18)
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

