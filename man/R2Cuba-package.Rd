\name{R2Cuba-package}
\alias{R2Cuba-package}
\alias{R2Cuba}
\docType{package}
\title{
\packageTitle{R2Cuba}
}
\description{
\packageDescription{R2Cuba}
}
\author{
\packageAuthor{R2Cuba}

Maintainer: \packageMaintainer{R2Cuba}
}
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
