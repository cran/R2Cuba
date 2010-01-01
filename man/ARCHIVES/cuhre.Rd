\name{cuhre}
\alias{cuhre}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Multidimensional numerical integration with a deterministic algorithm }
\description{
  Implement a  deterministic algorithm for multidimensional numerical
  integration. Its algorithm uses one of several cubature rules of
  polynomial degree in a globally adaptive subdivision scheme.
  The subdivision algorithm is similar to
  \code{\link{suave}}'s.
}
\usage{
cuhre(ndim, ncomp, integrand, ...,
     lower=rep(0,ndim), upper=rep(1,ndim),
     rel.tol= .Machine$double.eps^0.25,  abs.tol = rel.tol,
     flags=list(verbose=1, final=1, pseudo.random=0, mersenne.seed=NULL),
      min.eval=0,  max.eval=50000,  key=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ndim}{ the number of dimensions of the integral}
  \item{ncomp}{ the number of components  of the integrand}
  \item{integrand}{ the R function which computes the integrand.
    It is expected to be declared as
    
    \code{integrand <- function(x, ...)} or

    \code{integrand <- function(x, phw, ...)}
  
    where \code{x} is an input vector of length \code{ndim},
    and \code{phw} an ignored argument for compatibility with
    the other 'R2Cuba' functions.

    \code{...}  denotes optional  additional arguments which
    correspond to those
    passed to   the main function in \dQuote{...}.
    
    The value returned by this R function should be a vector of length
    \code{ncomp}.}
  \item{...}{ optional  additional parameters to be passed to
    \code{integrand}, if any }
  \item{lower}{ the lower bounds of the integration region.
    Vector of length \code{ndim}}
  \item{upper}{ the  upper bounds of the integration region.
  Vector of length \code{ndim}}
  \item{rel.tol}{ the requested relative accuracy}
  \item{abs.tol}{ the requested absolute accuracy}
  \item{flags}{flags governing the integration. A list
    with components:
    
 - \code{verbose}:  \code{verbose} encodes the verbosity level, from. 0 to 3.
  Level 0 does not print any output, level 1 prints \dQuote{reasonable} information on the
  progress of the integration, level 2 also echoes the input parameters, and level
  3 further prints the subregion results.

- \code{final}: when \code{ 0}, all sets of samples collected on a subregion during the
  various iterations or phases contribute to the final result.
  When \code{ 1}, only the last (largest) set of samples is used in the final
  result.

- \code{ pseudo.random}:
when \code{  0}, Sobol quasi-random numbers are used for sampling.
When \code{ 1}, Mersenne Twister pseudo-random numbers are used for
sampling.

- \code{ mersenne.seed}:
the seed for the Mersenne Twister algorithm, when \code{
  pseudo.random=1} and when it would be explicitly set.
}
  \item{min.eval}{ the minimum number of integrand evaluations required.}
  \item{max.eval}{ the (approximate) maximum number of integrand
    evaluations allowed.}
  \item{key}{ chooses the basic integration rule:
    
\code{key = 7, 9, 11, 13} selects the cubature rule of degree key. Note that the degree-11
rule is available only in 3 dimensions, the degree-13 rule only in 2
dimensions.

For other values, the default rule is taken, which is the degree-13 rule in 2 dimensions,
the degree-11 rule in 3 dimensions, and the degree-9 rule otherwise.}
}
\details{
  See details in the documentation.
}
\value{
  A list of the S3-class \code{cuba} with components:
  \item{nregions }{the actual number of subregions needed.}
  \item{neval }{the actual number of integrand evaluations needed.}
  \item{ifail}{ an error flag:

    \code{ifail=0} , the desired accuracy was reached,

    \code{ifail= -1}, dimension out of range,

    \code{ifail = 1}, the accuracy goal was not met within the allowed
    maximum number of integrand evaluations.
    

  }
  \item{value}{vector of length \code{ncomp};
               the integral of \code{integrand} over the 
	       hypercube.}
	     \item{abs.error}{vector of length \code{ncomp};
	       the presumed absolute error of \code{value}.}
	     \item{prob}{vector of length \code{ncomp};
              the \eqn{$\chi^2$}{Chi2}-probability (not the \eqn{$\chi^2$}{Chi2}-value itself!)
that \code{abs.error} is not a reliable estimate of the true integration
error.}
\item{message}{\dQuote{OK} or a character string giving the error message.}
\item{call}{The matched call.}
}
\seealso{  \code{\link{vegas}}, \code{\link{suave}}, \code{\link{divonne}} }
\examples{
integrand <- function(arg) {
  x <- arg[1]
  y <- arg[2]
  z <- arg[3]
  ff <- sin(x)*cos(y)*exp(z);
return(ff)
} # End integrand

NDIM <- 3
NCOMP <- 1
cuhre(NDIM, NCOMP,  integrand,
             rel.tol= 1e-3, abs.tol= 1e-12,
             flags= list(verbose=2, final=0))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ math }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
