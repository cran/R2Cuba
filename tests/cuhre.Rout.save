
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

> #dyn.load("../src/R2Cuba.so");source("../R/cuhre.R"); source("../R/commoncuba.R")
> library("R2Cuba")
> # DEMO EXAMPLES
> # +++++++++++++++++++++++++++++++
> integrand <- function(arg) {
+   x <- arg[1]
+   y <- arg[2]
+   z <- arg[3]
+   ff <- sin(x)*cos(y)*exp(z);
+ return(ff)
+ } # End integrand
> 
> 
> 
> NDIM <- 3
> NCOMP <- 1
> EPSREL <- 1e-3
> EPSABS <- 1e-12
> VERBOSE <-  2
> LAST <- 1
> 
> 
>  cuhre(NDIM, NCOMP,  integrand,  flags= list(verbose=VERBOSE))
Cuhre input parameters:
  ndim 3
  ncomp 1
  rel.tol 0.001
  abs.tol 4.93038e-32
  pseudo.random  0
  final 0
  verbose 2
  min.eval 0
  max.eval 50000
  key 0
Iteration 1:  127 integrand evaluations so far
[1] 0.66467 +- 7.2682e-10  	chisq 0 (0 df)
Iteration 2:  381 integrand evaluations so far
[1] 0.66467 +- 3.33018e-11  	chisq 0 (1 df)
integral: 0.6646697 (+-3.3e-11)
nregions: 2; number of evaluations:  381; probability:  0 
>  cuhre(NDIM, NCOMP,  integrand,  flags= list(pseudo.random=1))
Iteration 1:  127 integrand evaluations so far
[1] 0.66467 +- 7.2682e-10  	chisq 0 (0 df)
Iteration 2:  381 integrand evaluations so far
[1] 0.66467 +- 3.33018e-11  	chisq 0 (1 df)
integral: 0.6646697 (+-3.3e-11)
nregions: 2; number of evaluations:  381; probability:  0 
> 
> # Essai en déplacant les bornes
> integrand2 <- function(arg) {
+   x <- arg[1]
+   y <- arg[2]
+   z <- arg[3]
+   ff <- sin(x-3)*cos(y-2)*exp(z-1);
+ return(ff)
+ } # End integrand2
>  cuhre(NDIM, NCOMP,  integrand2, lower=c(3,2,1), upper=c(4,3,2), flags=list(verbose=3))
Cuhre input parameters:
  ndim 3
  ncomp 1
  rel.tol 0.001
  abs.tol 4.93038e-32
  pseudo.random  0
  final 0
  verbose 3
  min.eval 0
  max.eval 50000
  key 0
Region (3.000000) - (4.000000)
       (2.000000) - (3.000000)
       (1.000000) - (2.000000)
[1] 0.66467 +- 7.2682e-10
Iteration 1:  127 integrand evaluations so far
[1] 0.66467 +- 7.2682e-10  	chisq 0 (0 df)
Region (3.000000) - (4.000000)
       (2.000000) - (3.000000)
       (1.000000) - (1.500000)
[1] 0.25094 +- 1.25788e-11
Region (3.000000) - (4.000000)
       (2.000000) - (3.000000)
       (1.500000) - (2.000000)
[1] 0.41373 +- 2.0739e-11
Iteration 2:  381 integrand evaluations so far
[1] 0.66467 +- 3.33019e-11  	chisq 0 (1 df)
integral: 0.6646697 (+-3.3e-11)
nregions: 2; number of evaluations:  381; probability:  0 
> 
> 
> proc.time()
   user  system elapsed 
  0.156   0.023   0.173 
