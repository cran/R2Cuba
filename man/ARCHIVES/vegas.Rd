\name{vegas}
\alias{vegas}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Multidimensional numerical integration with a Monte Carlo algorithm}
\description{
  Implement a  Monte Carlo  algorithm for multidimensional numerical
  integration.     This algorithm  uses importance sampling as a variance-reduction
technique. Vegas iteratively builds up a piecewise constant weight function, represented
on a rectangular grid. Each iteration consists of a sampling step followed by a refinement
of the grid. }
\usage{
vegas(ndim, ncomp, integrand, ...,
      lower=rep(0,ndim), upper=rep(1,ndim),
      rel.tol= .Machine$double.eps^0.25,  abs.tol = rel.tol,
      flags=list(verbose=1, final=1, pseudo.random=0, smooth=0, mersenne.seed=NULL),
      min.eval=0,  max.eval=50000,
      nstart=1000,  nincrease=500, nbatch=1000, gridno=0, state.file=NULL)
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
  But \code{flags} may have an additional component:
 \code{smooth}. When \code{smooth = 0}, apply additional smoothing to the importance
   function, this moderately improves convergence for many integrands.
   When \code{smooth = 1} , use the importance function without smoothing, this should be chosen
if the integrand has sharp edges.

\emph{Note}:  Value 3 of \code{flags$verbose} has the same effect  as
value 2 (Vegas does not partition the integration region).
}
\item{min.eval}{ same as    \code{\link{cuhre}}}
\item{max.eval}{ same as    \code{\link{cuhre}}}
  \item{nstart}{ the number of integrand evaluations per iteration to
    start with.}
  \item{nincrease}{ the increase in the number of integrand evaluations
    per iteration.}
 \item{nbatch}{    Vegas samples points not all at once, but in batches of a predetermined size, to avoid
excessive memory consumption. \code{nbatch} is the number of
points sampled in each batch. Tuning this number should usually not be necessary as performance
is affected significantly
only as far as the batch of samples fits into the CPU cache.}

\item{gridno}{ an integer.  Vegas
may accelerate convergence to keep the grid accumulated during one integration for
     the next one, if the integrands are reasonably similar to each other. Vegas maintains
     an internal table with space for ten grids for this purpose. 
If \code{gridno} is a number between 1 and 10, the grid is not discarded at the end of
the integration, but stored in the respective slot of the table for a future invocation.
The grid is only re-used if the dimension of the subsequent integration is the same
as the one it originates from.
In repeated invocations it may become necessary to flush a slot in memory. In this
case the negative of the grid number should be set. Vegas will then start with a new
grid and also restore the grid number to its positive value, such that at the end of
the integration the grid is again stored in the indicated slot.}
\item{state.file}{  the name of an external file. Vegas    can store its entire internal state (i.e. all the information to resume an
  interrupted integration)   in an external file.


 The state file is updated after every iteration. If, on a subsequent
 invocation, Vegas finds a file of the specified name, it loads the internal state and continues from the
point it left off. Needless to say, using an existing state file with a different integrand
generally leads to wrong results. Once the integration finishes successfully, i.e. the
prescribed accuracy is attained, the state file is removed.
This feature is useful mainly to define \sQuote{check-points} in long-running integrations
from which the calculation can be restarted.}
}
\details{
    See details in the documentation.
}
\value{
  Idem as \code{\link{cuhre}}, except from \code{nregions} (not present)
}
\seealso{   \code{\link{cuhre}}, \code{\link{suave}}, \code{\link{divonne}}}
\examples{
integrand <- function(arg, weight) {
  x <- arg[1]
  y <- arg[2]
  z <- arg[3]
  ff <- sin(x)*cos(y)*exp(z);
return(ff)
} # end integrand
vegas(3, 1, integrand, rel.tol=1e-3,  abs.tol=1e-12, flags=list(verbose=2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ math }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
