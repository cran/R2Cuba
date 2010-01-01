#dyn.load("../src/cuba.so");source("../R/divonne.R");source("../R/cuhre.R");source("../R/vegas.R");source("../R/suave.R");source("../R/commoncuba.R")
library("R2Cuba")

#----------------------------------------
# Generation of the configuration of the centers:
set.seed(5)
ndim <- 2
ncomp <- 1
n.pts <- 50
pattern <- matrix(runif(ndim*n.pts),nrow=ndim)
sigma <- 0.01
rownames(pattern) <- c("x","y")

#----------------------------------------
# Exact computation of the integral:
int.val <- function(pattern,sigma) {
  fs <- array(pnorm(c(1,0),rep(pattern,each=2),sigma),c(2,2,ncol(pattern)))
  fs.diff <- apply(fs,c(2,3),diff)
  sum(apply(fs.diff,2,prod))
}
cat("Expected result:\n")
print(int.val(pattern,sigma)) # 49.66977
#----------------------------------------
#Function integrand: 
pat.dist <- function(x, sigma) {
  sum(apply(matrix(dnorm(x,pattern,sigma),2),2,prod))
}
#----------------------------------------
# Common arguments:
flags <- list(verbose=0,final=1)
rel.tol <- 0.1 # Not so good, for timing purpose
#ptm <- proc.time() # total timing 
#----------------------------------------
cat("\ncuhre\n")
 print(cuhre(ndim,ncomp,pat.dist,sigma,
    rel.tol=rel.tol,flags=flags))

#----------------------------------------
cat("\nvegas\n")
print(vegas(ndim,ncomp,pat.dist,sigma,
       rel.tol=rel.tol,flags=flags))
#----------------------------------------
cat("\nsuave\n")
print(suave(ndim,ncomp,pat.dist,sigma,
      rel.tol=rel.tol,flags=flags))
#----------------------------------------
cat("\ndivonne\n")
print(divonne(ndim,ncomp,pat.dist,sigma,
          rel.tol=rel.tol,flags=flags))
#----------------------------------------
cat("\ndivonne with  xgiven = the centers\n")
print(divonne(ndim,ncomp,pat.dist,sigma,xgiven=pattern,
                rel.tol=rel.tol,flags=flags))            
#----------------------------------------
cat("\ndivonne with  xgiven = the centers and the peaks\n")
# Vertices creation:
require("deldir")
voronoi <- deldir(pattern["x",],pattern["y",],rw=c(0,1,0,1))
vertices <- unique(rbind(as.matrix(voronoi$dirsgs[,1:2]),
                         as.matrix(voronoi$dirsgs[,3:4])))
# The maxima:
peaks <- rbind(vertices,t(pattern))
print(divonne(ndim,ncomp,pat.dist,sigma,xgiven=t(peaks),
                  rel.tol=rel.tol,flags=flags))             
#----------------------------------------
cat("\ndivonne with  peakfinder\n")
peakfinder <- function(bounds) {
  idx <- (peaks[,1]>=bounds[1,1]) & (peaks[,1]<=bounds[1,2])
  idx <- idx & (peaks[,2]>=bounds[2,1]) & (peaks[,2]<=bounds[2,2])
  retour <- t(peaks[idx,])
  if (ncol(retour)==0) {
    retour <- bounds
  }
  return(retour)
}
print(divonne(ndim,ncomp,pat.dist,sigma,
                          peakfinder=peakfinder,nextra=nrow(peaks),
                rel.tol=rel.tol,flags=flags))            
#----------------------------------------
#cat("Total time", (proc.time() - ptm),"\n")
#----------------------------------------

