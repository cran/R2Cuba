\name{R2Cuba-package}
\alias{R2Cuba-package}
\alias{R2Cuba}
\docType{package}
\title{Multidimensional Numerical Integration}
\description{
  R2Cuba implements four general-purpose multidimensional integration
  algorithms: Vegas, Suave, Divonne and Cuhre. It is a wrapper
  around the Cuba-1.6 library by Thomas Hahn available from the URL
  \url{http://www.feynarts.de/cuba/}.
}
\details{
\tabular{ll}{
Package: \tab R2Cuba\cr
Version: \tab 1.0-6\cr
Date: \tab 2010-09-02\cr
Depends: \tab methods\cr
Suggests:\tab deldir\cr
License: \tab GPL (>= 3)\cr}}
\author{
  The Cuba library has been written by Thomas Hahn
  (\url{http://wwwth.mppmu.mpg.de/members/hahn}).
   Interface to R was written by Annie Bouvier and Kiên Kiêu
   (MIA Lab, INRA, Jouy-en-Josas, France
   \url{http://www.jouy.inra.fr/mia_eng/}).

Maintainer: Annie Bouvier <Annie.Bouvier@jouy.inra.fr>}
\references{
  The Cuba library is described at
  \url{http://www.feynarts.de/cuba/}.
User documentation is available in \cite{T. Hahn (2005) CUBA-a library for multidimensional numerical
   integration. \emph{Computer Physics Communications}, \bold{168},
   78-95.} (\url{http://arxiv.org/pdf/hep-ph/0404043}).}
\seealso{The R-package \dQuote{cubature}}
\keyword{ package }
\examples{
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
}
