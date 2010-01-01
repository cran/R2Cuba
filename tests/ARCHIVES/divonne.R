#dyn.load("../src/cuba.so")
#source("../R/divonne.R")
library("Cuba")
source("integrand.R")
# VOIR: tester integrand avec l'argu  phase  (après ncomp); il est ici ignoré

NDIM <- 3
NCOMP <- 1
EPSREL <- 1e-3
EPSABS <- 1e-12
VERBOSE <- 2
MINEVAL <- 0
MAXEVAL <- 50000
KEY1 <- 47
KEY2 <- 1
KEY3 <- 1
MAXPASS <- 5
BORDER <- 0.
MAXCHISQ <- 10.
MINDEVIATION <- .25
NGIVEN <- 0
LDXGIVEN <- NDIM
NEXTRA <- 0
ret <- divonne(NDIM, NCOMP, integrand,
                  EPSREL,  EPSABS, VERBOSE ,
                  MINEVAL,  MAXEVAL, KEY1, KEY2, KEY3, MAXPASS, BORDER, MAXCHISQ, MINDEVIATION,
    NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL)
# VOIR: les argu qui sont ici NULL
print(ret)
