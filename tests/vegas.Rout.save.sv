
R version 2.15.0 (2012-03-30)
Copyright (C) 2012 The R Foundation for Statistical Computing
ISBN 3-900051-07-0
Platform: x86_64-unknown-linux-gnu (64-bit)

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

> #dyn.load("../src/R2Cuba.so");source("../R/vegas.R"); source("../R/commoncuba.R")
> library("R2Cuba")
> 
> # DEMO EXAMPLES
> # +++++++++++++++++++++++++++++++
> integrand <- function(arg, weight) {
+ 
+   x <- arg[1]
+   y <- arg[2]
+   z <- arg[3]
+ #  print( weight)
+     ff <- sin(x)*cos(y)*exp(z);
+ return(ff)
+ } # fin integrand
> # ----------------------------------------
> 
> NDIM <- 3
> NCOMP <- 1
> EPSREL <- 1e-3
> EPSABS <- 1e-12
> VERBOSE <- 0
> 
> vegas(NDIM, NCOMP, integrand,
+                   rel.tol=EPSREL,  abs.tol=EPSABS,
+              flags=list(verbose=VERBOSE))
integral: 0.6648107 (+-0.00049)
number of evaluations:  10000; probability:  0.08870493 
> # ----------------------------------------
> # Essai avec ...
> integrandp <- function(arg, weight, ...) {
+ 
+   x <- arg[1]
+   y <- arg[2]
+   z <- arg[3]
+  # print( weight, ...)
+     ff <- sin(x)*cos(y)*exp(z);
+ return(ff)
+ } # fin integrandp
> vegas(NDIM, NCOMP, integrandp, digits=1,
+                   rel.tol=EPSREL,  abs.tol=EPSABS,
+              flags=list(verbose=0))
integral: 0.6648107 (+-0.00049)
number of evaluations:  10000; probability:  0.08870493 
> # ----------------------------------------
> # Essai avec state
> vegas(NDIM, NCOMP, integrand,
+                   rel.tol=1e-50,  abs.tol=1e-50,
+              flags=list(verbose=VERBOSE), state.file="toto")
integral: 0.6646931 (+-0.00013)
number of evaluations:  67500; probability:  0.005078484 
failed with message  'Dimension out of range' 
> 
> vegas(NDIM, NCOMP, integrand,
+                    rel.tol=EPSREL,  abs.tol=EPSABS,
+              flags=list(verbose=VERBOSE), state.file="toto")
integral: 0.6646589 (+-0.00012)
number of evaluations:  76000; probability:  0.01009976 
> 
> # ----------------------------------------
> # Essai avec mersenne
> vegas(NDIM, NCOMP, integrand,
+   rel.tol=EPSREL,  abs.tol=EPSABS,
+       flags=list(verbose=0, pseudo.random=1, mersenne.seed=10))
integral: 0.6640216 (+-0.00065)
number of evaluations:  13500; probability:  0.07234709 
> 
> # ----------------------------------------
> # Essai avec nbatch
> vegas(NDIM, NCOMP, integrand,
+                   rel.tol=EPSREL,  abs.tol=EPSABS,
+              flags=list(verbose=0, vegas.nbatch=500))
integral: 0.6648107 (+-0.00049)
number of evaluations:  10000; probability:  0.08870493 
> # ----------------------------------------
> # Essai avec gridno
> #dyn.load("../src/cuba.so");source("../R/vegas.R"); source("../R/commoncuba.R")
> vegas(NDIM, NCOMP, integrand,
+                   rel.tol=EPSREL,  abs.tol=EPSABS,
+              flags=list(verbose=3, pseudo.random=1, mersenne.seed=10,
+                vegas.gridno=4))
Vegas input parameters:
  ndim 3
  ncomp 1
  rel.tol 0.001
  abs.tol 1e-12
  smooth 0
  pseudo.random  1
  final 0
  verbose 3
  min.eval 0
  max.eval 50000
  nstart 1000
  nincrease 500
  vegas.gridno 0
  vegas.state ""
Iteration 1:  1000 integrand evaluations so far
[1] 0.666187 +- 0.014126  	chisq 0 (0 df)
Iteration 2:  2500 integrand evaluations so far
[1] 0.663884 +- 0.0048976  	chisq 0.0302145 (1 df)
Iteration 3:  4500 integrand evaluations so far
[1] 0.6631 +- 0.00224811  	chisq 0.0626448 (2 df)
Iteration 4:  7000 integrand evaluations so far
[1] 0.662786 +- 0.00127313  	chisq 0.0914002 (3 df)
Iteration 5:  10000 integrand evaluations so far
[1] 0.663624 +- 0.000891107  	chisq 0.942088 (4 df)
Iteration 6:  13500 integrand evaluations so far
[1] 0.664022 +- 0.000651656  	chisq 1.3692 (5 df)
integral: 0.6640216 (+-0.00065)
number of evaluations:  13500; probability:  0.07234709 
> vegas(NDIM, NCOMP, integrand,
+                   rel.tol=EPSREL,  abs.tol=EPSABS,
+              flags=list(verbose=0, pseudo.random=1, mersenne.seed=10,
+                vegas.gridno=4))
integral: 0.6640216 (+-0.00065)
number of evaluations:  13500; probability:  0.07234709 
> vegas(NDIM, NCOMP, integrand,
+                   rel.tol=EPSREL,  abs.tol=EPSABS,
+              flags=list(verbose=0, pseudo.random=1, mersenne.seed=10,
+                vegas.gridno=-4))
integral: 0.6640216 (+-0.00065)
number of evaluations:  13500; probability:  0.07234709 
> 
> proc.time()
   user  system elapsed 
  0.612   0.023   0.629 
