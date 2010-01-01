#dyn.load("../src/R2Cuba.so");source("../R/cuhre.R"); source("../R/commoncuba.R")
library("R2Cuba")
# DEMO EXAMPLES
# +++++++++++++++++++++++++++++++
integrand <- function(arg) {
  x <- arg[1]
  y <- arg[2]
  z <- arg[3]
  ff <- sin(x)*cos(y)*exp(z);
return(ff)
} # End integrand



NDIM <- 3
NCOMP <- 1
EPSREL <- 1e-3
EPSABS <- 1e-12
VERBOSE <-  2
LAST <- 1


 cuhre(NDIM, NCOMP,  integrand,  flags= list(verbose=VERBOSE))
 cuhre(NDIM, NCOMP,  integrand,  flags= list(pseudo.random=1))

# Essai en déplacant les bornes
integrand2 <- function(arg) {
  x <- arg[1]
  y <- arg[2]
  z <- arg[3]
  ff <- sin(x-3)*cos(y-2)*exp(z-1);
return(ff)
} # End integrand2
 cuhre(NDIM, NCOMP,  integrand2, lower=c(3,2,1), upper=c(4,3,2), flags=list(verbose=3))

