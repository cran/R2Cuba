#dyn.load("../src/cuba.so")
#source("../R/cuhre.R"); source("../R/methods.R")
library("Cuba")
source("integrand.R")

NDIM <- 3
NCOMP <- 1
EPSREL <- 1e-3
EPSABS <- 1e-12
VERBOSE <-  2
LAST <- 4
MINEVAL <- 0
MAXEVAL <- 50000
KEY <- 0
ret <- cuhre(NDIM, NCOMP, integrand,
                  EPSREL,  EPSABS, VERBOSE| LAST,
                  MINEVAL,  MAXEVAL, KEY)
print(ret)
