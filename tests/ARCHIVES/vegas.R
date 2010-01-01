#dyn.load("../src/cuba.so")
#source("../R/vegas.R")
library("Cuba")
source("integrand.R")
# VOIR: tester integrand avec l'argu weight (après ncomp); il est ici ignoré


NDIM <- 3
NCOMP <- 1
EPSREL <- 1e-3
EPSABS <- 1e-12
VERBOSE <- 2
MINEVAL <- 0
MAXEVAL <- 50000
NSTART <- 1000
NINCREASE <- 500

ret <- vegas(NDIM, NCOMP, integrand,
                  EPSREL,  EPSABS, VERBOSE,
                  MINEVAL,  MAXEVAL,
                  NSTART,  NINCREASE)
print(ret)
