divonne <- function(ndim, ncomp,
                   integrand, ...,
                 lower=rep(0,ndim), upper=rep(1,ndim), 
                  rel.tol= 1e-3,
                  abs.tol = 0,
                  flags=list(verbose=1,
                    final=1, pseudo.random=0,
                    mersenne.seed=NULL),
                    min.eval=0,  max.eval=50000,
                    key1=47,
                    key2=1, key3=1, 
                    max.pass=5, border=0, max.chisq=10,
                    min.deviation=0.25,
                    xgiven=NULL,  nextra=0,
                    peakfinder=NULL)
                    
  {
     # Verification
 if (!verif(ndim, ncomp, lower, upper, rel.tol, abs.tol,
            flags, min.eval,  max.eval))
      stop("Error in input: see the warnings")
     # Decode the flags
    lesflags <- decodflags(flags)
if (is.null(flags$mersenne.seed))
  flags$mersenne.seed <- NA
    
 #  xgiven dimensions:
 if (!is.null(xgiven)) {
   if (!is.matrix(xgiven))
     stop("xgiven should be a matrix")
   ngiven <- ncol(xgiven)
   ldxgiven <- nrow(xgiven)
   if (ldxgiven != ndim)
     stop("Matrix xgiven should have ndim rows")
   # Rescaler xgiven dans l'hypercube unité
   for (i in 1:ndim) {
     xgiven[i,] <- xgiven[i,]/(upper[i]-lower[i])
   }
 }
 else {
    ngiven=0; ldxgiven=ndim
  }
      if (nextra <0)
       stop("nextra should be positive or zero")
 
 if ((nextra >0) && is.null(peakfinder))
   stop("nextra positive but not peakfinder subroutine")
 if ((nextra==0) && !is.null(peakfinder))
   warning("peakfinder ignored: nextra is zero")
# 6/6/2012
 if (is.null(peakfinder))
   peakfinder=0
 
  
       # Allocate outputs
    nregions <- 0
    neval <- 0
    fail <- 0
    integral <- rep(0, ncomp)
    error <- rep(0, ncomp)
    prob <- rep(0, ncomp)
  libargs <- c("ndim", "ncomp", 
                 "integrand","lower", "upper",
                  "rel.tol", "abs.tol", "flags",
                    "min.eval","max.eval", "key1",
                    "key2", "key3", 
                    "max.pass", "border", "max.chisq",
                    "min.deviation", "xgiven",  "nextra",
                    "peakfinder")

 ffintegrand <- crff(match.call(), integrand, "divonne", libargs, ...)
  prdbounds <- prod(upper-lower)
  ret <-  .C("Rdivonne", as.integer(ndim),
             as.integer(ncomp),
             ffintegrand, new.env(),
             as.double(lower), as.double(upper),as.double(prdbounds),
             as.double(rel.tol), as.double(abs.tol),
             as.integer(flags$mersenne.seed),
             as.integer(lesflags),
             as.integer(min.eval),  as.integer(max.eval),
             as.integer(key1), as.integer(key2), as.integer(key3),
             as.integer(max.pass), as.double(border),
             as.double(max.chisq), as.double(min.deviation),
             as.integer(ngiven), as.integer(ldxgiven),
             as.double(xgiven), as.integer(nextra),
             peakfinder,
             nregions=as.integer(nregions),
             neval=as.integer(neval), fail=as.integer(fail),
             integral=as.double(integral),
             error=as.double(error), prob=as.double(prob),
             NAOK=TRUE)
#Add to finish the last print:
cat("\n")

# To homogeneize with the R function "integrate", add
    # message and call into the output,
    # ifail rather than fail , abs.error rather than  error,
    # value rather than integral    
    
    if (ret$fail ==0)
      mess ="OK" # OK required to be printed by print.cuba
    else
      mess="Dimension out of range"

    retour <- list(method="divonne",
                   neval=ret$neval, nregions=ret$nregions,
                   ifail=ret$fail, value=ret$integral,
                abs.error=ret$error,
                   neval=ret$neval,
                   prob=ret$prob, message=mess,
                call=match.call())
    attr(retour, "class") = c("cuba")
    return(retour)
  } # End divonne

                
