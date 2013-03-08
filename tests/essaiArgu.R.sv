#dyn.load("../src/R2Cuba.so");source("../R/cuhre.R"); source("../R/commoncuba.R"); source("../R/divonne.R");source("../R/suave.R");source("../R/vegas.R")
library("R2Cuba")
# Essai d'integrandes avec diverses configurations d'arguments
# +++++++++++++++++++++++++++++++
integrand <- function(arg, a, sigma, toto) {
  x <- arg[1]
  y <- arg[2]
  z <- arg[3]
  ff <- sin(x)*cos(y)*exp(z);
#  print(sigma)
return(ff)
} # End integrand

integrand2 <- function(arg) {
    x <- arg[1]
  y <- arg[2]
  z <- arg[3]
  ff <- sin(x)*cos(y)*exp(z);
 return(ff)
} # End integrand2

integrand3 <- function(arg, phase) {
    x <- arg[1]
  y <- arg[2]
  z <- arg[3]
  ff <- sin(x)*cos(y)*exp(z);
  #  print(phase)
 return(ff)
} # End integrand3

integrand4 <- function(arg, phase, sigma=0) {
    x <- arg[1]
  y <- arg[2]
  z <- arg[3]
#  print(sigma)
    ff <- sin(x)*cos(y)*exp(z);
 return(ff)
} # End integrand4

NDIM <- 3
NCOMP <- 1

# --------------------divonne ------------------------------
# plante: sigma manquant dans l'appel : print(divonne(NDIM, NCOMP,  integrand,  flags= list(verbose=0)))
print(divonne(NDIM, NCOMP,  integrand2,  flags= list(verbose=0)))
print(divonne(NDIM, NCOMP,  integrand3,  flags= list(verbose=0)))
# plante: sigma manquant dans l'integrande: print(divonne(NDIM, NCOMP,  integrand2,  flags= list(verbose=0), sigma=1))
print(divonne(NDIM, NCOMP,  integrand4,  flags= list(verbose=0), sigma=1))

# --------------------cuhre ------------------------------
# plante: sigma manquant dans l'appel : print(cuhre(NDIM, NCOMP,  integrand,  flags= list(verbose=0)))
print(cuhre(NDIM, NCOMP,  integrand2,  flags= list(verbose=0)))
print(cuhre(NDIM, NCOMP,  integrand3,  flags= list(verbose=0)))
# plante: sigma manquant dans l'integrande: print(cuhre(NDIM, NCOMP,  integrand2,  flags= list(verbose=0), sigma=1))
print(cuhre(NDIM, NCOMP,  integrand4,  flags= list(verbose=0), sigma=1))
print(cuhre(3,1,  integrand4,  flags= list(verbose=0), sigma=1))

# --------------------suave ------------------------------
print(suave(NDIM, NCOMP,  integrand2,  flags= list(verbose=0)))
# --------------------vegas ------------------------------
print(vegas(NDIM, NCOMP,  integrand2,  flags= list(verbose=0)))
integrandp <- function(arg, weight, ...) {

  x <- arg[1]
  y <- arg[2]
  z <- arg[3]
 # print( weight, ...)
    ff <- sin(x)*cos(y)*exp(z);
return(ff)
} # fin integrandp
EPSREL <- 1e-3
EPSABS <- 1e-12
print(vegas(NDIM, NCOMP, integrandp, digits=1,
                  rel.tol=EPSREL,  abs.tol=EPSABS,
             flags=list(verbose=0)))
# ----------------------------------------
# Essai avec state
print(vegas(NDIM, NCOMP, integrand4,
                  rel.tol=1e-50,  abs.tol=1e-50,
             flags=list(verbose=0), state.file="toto"))
print(vegas(NDIM, NCOMP, integrand,
                  rel.tol=EPSREL,  abs.tol=EPSABS,
             flags=list(verbose=3, pseudo.random=1, mersenne.seed=10,
               vegas.gridno=4), sigma=5))
