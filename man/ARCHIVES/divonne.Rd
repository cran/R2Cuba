\name{divonne}
\alias{divonne}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{  Multidimensional numerical integration by stratified sampling for variance reduction}
\description{
  Divonne works by stratified sampling, where the partioning of the
integration
region is aided by methods from numerical optimization.}
\usage{
divonne(ndim, ncomp, integrand, ...,
        lower=rep(0,ndim), upper=rep(1,ndim), 
        rel.tol= .Machine$double.eps^0.25,   abs.tol = rel.tol,
        flags=list(verbose=1, final=1, pseudo.random=0, mersenne.seed=NULL),
        min.eval=0,  max.eval=50000,
        key1=47,  key2=1, key3=1, 
        max.pass=5, border=0, max.chisq=10,
        min.deviation=0.25,
        xgiven=NULL,  nextra=0, peakfinder=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{ndim}{ same as  \code{\link{cuhre}}}
\item{ ncomp}{ same as  \code{\link{cuhre}}}
\item{ integrand}{ same as  \code{\link{cuhre}}.
  But, here, the input argument \code{phw} indicates the integration
  phase:
  
  \code{0}: sampling of the points in \code{xgiven},
  
  \code{1}: partitioning phase,
  
  \code{2}: final integration phase,
  
   \code{3}:  refinement phase.
   
This information might be useful if the integrand takes long to compute and a sufficiently
accurate approximation of the integrand is available. The actual value of the integral is only
of minor importance in the partitioning phase, which is instead much more dependent on
the peak structure of the integrand to find an appropriate tessellation. An approximation
which reproduces the peak structure while leaving out the fine details might hence be a
perfectly viable and much faster substitute when \code{phw < 2}.

In all other instances, \code{phw} can be ignored.
}
\item{ ...}{ same as  \code{\link{cuhre}}}
\item{lower}{ same as    \code{\link{cuhre}}}
\item{upper}{ same as    \code{\link{cuhre}}}
\item{rel.tol}{ same as    \code{\link{cuhre}}}
\item{abs.tol}{ same as    \code{\link{cuhre}}}
\item{flags}{ same as    \code{\link{cuhre}}}
\item{min.eval}{ same as    \code{\link{cuhre}}}
\item{max.eval}{ same as    \code{\link{cuhre}}}
  \item{key1}{   integer that                    determines sampling in
    the partitioning phase*:
    
\code{key1 = 7, 9, 11, 13} selects the cubature rule of degree \code{key1}. Note that the degree-11
rule is available only in 3 dimensions, the degree-13 rule only in 2 dimensions.
For other values of \code{key1}, a quasi-random sample of
\eqn{n=|key1|}{\code{n=|key1|}} points is used, where
the sign of \code{key1} determines the type of sample,

\code{key1 = 0}, use the default rule.

\code{key1 > 0}, use a Korobov quasi-random sample,

 \code{key1 < 0}, use a \dQuote{standard} sample (a Mersenne Twister pseudo-random sample
     if \code{flags$pseudo.random=1}, otherwise a Sobol quasi-random sample).}
   \item{key2}{ integer that                     determines sampling in
     the final integration phase:
     same as \code{key1}, but here
     \eqn{$n=|key2|$}{\code{n = |key2|}}
 determines the number of
 points, \eqn{n  > 39}{\code{n  > 39}}, sample \eqn{n} points,
 \eqn{n < 40}{\code{n < 40}}, sample \eqn{n}{\code{n}}
 \code{nneed} points, where \code{nneed} is the number of points needed to
      reach the prescribed accuracy, as estimated by Divonne from the results of the
      partitioning phase.}
    \item{key3}{ integer that sets the strategy for the refinement
      phase:
      
      \code{key3 = 0}, do not treat the subregion any further.
      
      \code{key3 = 1}, split the subregion up once more.
      
Otherwise, the subregion is sampled a third time with \code{key3} specifying the sampling
parameters exactly as \code{key2} above.}
  \item{max.pass}{   integer that                        controls the thoroughness of the partitioning phase: The
partitioning phase terminates when the estimated total number of integrand evaluations (partitioning plus final integration) does not decrease for \code{max.pass} successive
iterations.

A decrease in points generally indicates that Divonne discovered new structures of
the integrand and was able to find a more effective partitioning.
\code{max.pass} can be
understood as the number of \dQuote{safety} iterations that are performed before the partition is accepted as final and counting consequently restarts at zero whenever new
structures are found.}
  \item{border}{ the relative width of the border of the integration region.
Points falling into the border region will not be sampled directly, but will be extrapolated from two samples from the interior. Use a non-zero \code{border} if the integrand
subroutine cannot produce values directly on the integration
boundary. The  width of the border
is identical in all the dimensions.
For example, set \code{border=0.1} for a border of width equal
to 10\% of the width of the integration region.}
  \item{max.chisq}{ the maximum  \eqn{$\chi^2$}{Chi2} value  a single subregion  is
    allowed to have in the final integration phase. Regions which fail this \eqn{$\chi^2$}{Chi2} test and whose
    sample averages differ by more than \code{min.deviation}
    move on to the refinement phase.}
  \item{min.deviation}{ a bound, given as the fraction of the requested
    error of the entire integral, which determines whether it is
    worthwhile further examining a region that failed the
    \eqn{$\chi^2$}{Chi2} test.
    Only if the two sampling averages
obtained for the region differ by more than this bound is the region further treated.}
\item{xgiven}{ a matrix ( \code{ndim}, \code{ngiven}).
  A list of \code{ngiven} points where the
  integrand might have peaks.

    Divonne will consider these points when partitioning the
integration region. The idea here is to help the integrator find the extrema of the integrand in the presence of very narrow peaks. Even if only the approximate location
of such peaks is known, this can considerably speed up convergence.}
\item{nextra}{the maximum number of extra points the peak-finder
subroutine will return. If \code{nextra} is zero, \code{peakfinder} is
not called and an arbitrary object may be passed in its place, e.g. just
0.}
\item{peakfinder}{ the peak-finder subroutine. This R function is called
whenever a region is up for subdivision and is supposed to point out possible peaks
lying in the region, thus acting as the dynamic counterpart of the static list of points
supplied in \code{xgiven}. It is expected to be declared as

\code{peakfinder <- function(bounds)}

where \code{bounds} is a
matrix of dimension (\code{ndim, 2}) which contains
the upper and lower  bounds of the subregion. The names of the columns
are \code{c("lower", "upper")}.
    
The returned value should be a matrix (\code{ndim, nx})
where  \code{nx} is the actual number of
  points (should be less or equal to
\code{nextra}).} 
}
\details{
   Divonne uses stratified sampling for variance reduction, that is, it partitions the integration
region such that all subregions have an approximately equal value of a quantity called the
spread. 

See details in the documentation.
}
\value{
  Idem as \code{\link{cuhre}}.
  Here \code{ifail} may be \code{>1} when
  the accuracy goal was not met within the allowed maximum number of
integrand evaluations.  Divonne
can estimate the number of points by which
\code{maxeval} needs to be increased to
reach the desired accuracy and returns this value.
}
\seealso{  \code{\link{cuhre}}, \code{\link{suave}}, \code{\link{vegas}}}
\examples{
NDIM <- 3
NCOMP <- 1
integrand <- function(arg, phase) {
  x <- arg[1]
  y <- arg[2]
  z <- arg[3]
  ff <- sin(x)*cos(y)*exp(z);
return(ff)
}
divonne(NDIM, NCOMP, integrand, rel.tol=1e-3,  abs.tol=1e-12,
        flags=list(verbose=2),  key1= 47)

# Example with a peak-finder function
NMAX <- 4

peakf <- function(bounds) {
  print(bounds) # matrix (ndim,2)
  x <- matrix(0, ncol=NMAX, nrow=NDIM)
   pas <- 1/(NMAX-1)
   # 1ier point
   x[,1] <- rep(0, NDIM)
   # Les autres points
   for (i in 2:NMAX) {
      x[,i] <- x[,(i-1)] + pas
    }
  return(x)
} #end peakf

divonne(NDIM, NCOMP, integrand,
               flags=list(verbose=0) ,
                peakfinder=peakf, nextra=NMAX)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ math }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
