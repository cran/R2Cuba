\name{suave}
\alias{suave}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{  Multidimensional numerical integration with SUbregion-Adaptive Vegas algorithm }
\description{
  Suave uses \code{\link{vegas}}-like importance sampling combined
with a globally adaptive subdivision strategy: Until the requested accuracy is reached,
the region with the largest error at the time is bisected in the dimension in which the
fluctuations of the integrand are reduced most. The number of new samples in each half
is prorated for the fluctuation in that half.}
\usage{
suave(ndim, ncomp, integrand, ...,
      lower=rep(0,ndim), upper=rep(1,ndim),
      rel.tol= 0.001, abs.tol = 0,
      flags=list(verbose=1, final=1, pseudo.random=0, smooth=0, mersenne.seed=NULL),
      min.eval=0,  max.eval=50000,
      nnew=1000, flatness= 50)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{ndim}{ same as  \code{\link{cuhre}}}
\item{ ncomp}{ same as  \code{\link{cuhre}}}
\item{ integrand}{ same as  \code{\link{cuhre}};
  But, here, the input argument \code{phw}
  contains the weight of the point being sampled.
                                       This
				       extra value can safely be
				       ignored.
}
\item{ ...}{ same as  \code{\link{cuhre}}}
\item{lower}{ same as    \code{\link{cuhre}}}
\item{upper}{ same as    \code{\link{cuhre}}}
\item{rel.tol}{ same as    \code{\link{cuhre}}}
\item{abs.tol}{ same as    \code{\link{cuhre}}}
\item{flags}{ same as    \code{\link{cuhre}}.
    \code{flags} may have an additional component: 
 \code{smooth}. When \code{ = 0}, apply additional smoothing to the importance
   function, this moderately improves convergence for many integrands.
   When \code{= 1} , use the importance function without smoothing, this should be chosen
if the integrand has sharp edges.
}
\item{min.eval}{ same as    \code{\link{cuhre}}}
\item{max.eval}{ same as    \code{\link{cuhre}}}
  \item{nnew}{ the number of new integrand evaluations in each subdivision.}
  \item{flatness}{This parameter determines how
    prominently \dQuote{outliers}, i.e. individual samples with a large
    fluctuation,
    figure in the total fluctuation,
which in turn determines how a region is split up. As suggested by its name, \code{flatness}
should be chosen large for \dQuote{flat} integrands and small for \dQuote{volatile}
integrands with
high peaks. Note that since \code{flatness} appears in the exponent, one should not use
too large values (say, no more than a few hundred) lest terms be truncated internally
to prevent overflow. More details about this parameter can be found
Hahn's paper from 2005 and in Cuba documentation.}
}
\details{
    See details in the documentation.dériv
}
\value{
  Idem as \code{\link{cuhre}}
}
\references{
  T. Hahn (2005) CUBA-a library for multidimensional numerical
  integration. \emph{Computer Physics Communications}, \bold{168}, 78-95.
  }
\seealso{  \code{\link{cuhre}}, \code{\link{divonne}}, \code{\link{vegas}}}
\examples{
integrand <- function(arg, weight) {
  x <- arg[1]
  y <- arg[2]
  z <- arg[3]
  ff <- sin(x)*cos(y)*exp(z);
return(ff)
} # end integrand
suave(3, 1, integrand, rel.tol=1e-3,  abs.tol=1e-12,
             flags=list(verbose=2, final=0))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ math }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
