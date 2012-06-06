pkgname <- "R2Cuba"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('R2Cuba')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("R2Cuba-package")
### * R2Cuba-package

flush(stderr()); flush(stdout())

### Name: R2Cuba-package
### Title: Multidimensional Numerical Integration
### Aliases: R2Cuba-package R2Cuba
### Keywords: package

### ** Examples

integrand <- function(arg, weight) {
  x <- arg[1]
  y <- arg[2]
  z <- arg[3]
  ff <- sin(x)*cos(y)*exp(z);
return(ff)
} # end integrand
NDIM <-3
NCOMP <- 1
vegas(NDIM, NCOMP, integrand, rel.tol=1e-3,  abs.tol=1e-12)



cleanEx()
nameEx("cuhre")
### * cuhre

flush(stderr()); flush(stdout())

### Name: cuhre
### Title: Multidimensional numerical integration with a deterministic
###   iterative adaptive algorithm
### Aliases: cuhre
### Keywords: math

### ** Examples

integrand <- function(arg) {
  x <- arg[1]
  y <- arg[2]
  z <- arg[3]
  ff <- sin(x)*cos(y)*exp(z);
return(ff)
} # End integrand

NDIM <- 3
NCOMP <- 1
cuhre(NDIM, NCOMP,  integrand,
             rel.tol= 1e-3, abs.tol= 1e-12,
             flags= list(verbose=2, final=0))



cleanEx()
nameEx("divonne")
### * divonne

flush(stderr()); flush(stdout())

### Name: divonne
### Title: Multidimensional numerical integration by stratified sampling
###   for variance reduction
### Aliases: divonne
### Keywords: math

### ** Examples

NDIM <- 3
NCOMP <- 1
integrand <- function(arg, phase) {
  x <- arg[1]
  y <- arg[2]
  z <- arg[3]
  ff <- sin(x)*cos(y)*exp(z);
return(ff)
}
divonne(NDIM, NCOMP, integrand, rel.tol=1e-3,  abs.tol=1e-12,
        flags=list(verbose=2),  key1= 47)

# Example with a peak-finder function
NMAX <- 4

peakf <- function(bounds) {
#  print(bounds) # matrix (ndim,2)
  x <- matrix(0, ncol=NMAX, nrow=NDIM)
   pas <- 1/(NMAX-1)
   # 1ier point
   x[,1] <- rep(0, NDIM)
   # Les autres points
   for (i in 2:NMAX) {
      x[,i] <- x[,(i-1)] + pas
    }
  return(x)
} #end peakf

divonne(NDIM, NCOMP, integrand,
               flags=list(verbose=0) ,
                peakfinder=peakf, nextra=NMAX)



cleanEx()
nameEx("suave")
### * suave

flush(stderr()); flush(stdout())

### Name: suave
### Title: Multidimensional numerical integration with SUbregion-Adaptive
###   Vegas algorithm
### Aliases: suave
### Keywords: math

### ** Examples

integrand <- function(arg, weight) {
  x <- arg[1]
  y <- arg[2]
  z <- arg[3]
  ff <- sin(x)*cos(y)*exp(z);
return(ff)
} # end integrand
suave(3, 1, integrand, rel.tol=1e-3,  abs.tol=1e-12,
             flags=list(verbose=2, final=0))



cleanEx()
nameEx("vegas")
### * vegas

flush(stderr()); flush(stdout())

### Name: vegas
### Title: Multidimensional numerical integration with a Monte Carlo
###   algorithm
### Aliases: vegas
### Keywords: math

### ** Examples

integrand <- function(arg, weight) {
  x <- arg[1]
  y <- arg[2]
  z <- arg[3]
  ff <- sin(x)*cos(y)*exp(z);
return(ff)
} # end integrand
vegas(3, 1, integrand, rel.tol=1e-3,  abs.tol=1e-12, flags=list(verbose=2))



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
