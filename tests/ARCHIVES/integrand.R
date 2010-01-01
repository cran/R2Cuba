# DEMO EXAMPLES
# +++++++++++++++++++++++++++++++
integrand <- function(arg) {
ndim <- arg[1] # Not used  here
ncomp <- arg[ndim+1]# Not used here

  x <- arg[2]
  y <- arg[3]
  z <- arg[4]
  ff <- sin(x)*cos(y)*exp(z);
return(ff)
} # fin integrand
