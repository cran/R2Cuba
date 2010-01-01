#dyn.load("../src/cuba.so")
#source("../R/suave.R")
library("Cuba")
source("integrand.R")
# VOIR: tester integrand avec l'argu weight (après ncomp); il est ici ignoré

NDIM <- 3
NCOMP <- 1
EPSREL <- 1e-3
EPSABS <- 1e-12
VERBOSE <- 2
LAST <- 4
MINEVAL <- 0
MAXEVAL <- 50000
NNEW <- 1000
FLATNESS <- 25.

ret <- suave(NDIM, NCOMP, integrand,
                  EPSREL,  EPSABS, VERBOSE | LAST,
                  MINEVAL,  MAXEVAL, NNEW, FLATNESS)

print(ret)
