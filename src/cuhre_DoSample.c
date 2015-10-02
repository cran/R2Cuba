#include "common_stddecl.h"
#include "cuhre_util.h"

/*********************************************************************/
/* The function RIntegrand calls the R user function */
/*********************************************************************/
 void RIntegrand(ccount *ndim, ctreal xx[],
		       ccount  *ncomp,
		       ctreal *lower, ctreal *upper, ctreal prdbounds, 
		       real ff[], SEXP rho, SEXP globf, Glob *globdim)
{


/* rho, globf: La fonction R à intégrer et son environnement d'execution */		

  SEXP args, callsxp,resultsxp;
  int i;
 /*f:  the R function and its environment, rho are global */
   // The input arguments are the x
  PROTECT(args=allocVector(REALSXP, ( *ndim )));
  PROTECT(resultsxp=allocVector(REALSXP, ( *ncomp )));
  /* Fill in the input arguments with rescaling between   0-1,
     according to the bounds */
  for (i =0; i<*ndim; i++) {
    REAL(args)[i] = xx[i] * (globdim->upper_[i] - globdim->lower_[i]) + globdim->lower_[i];
  }

  /* Call the R function */

  PROTECT(callsxp=lang2( globf, args));
  // PROTECT(args=eval(callsxp,rho));
 PROTECT(resultsxp=eval(callsxp,rho));

 UNPROTECT(4);
 if  (length(resultsxp) != *ncomp) {
  error("Function integrand does not return a vector of length ncomp\n Length of returned vector= %d. ncomp=%d\n",
	length(resultsxp), *ncomp);
 }

 for (i =0; i<*ncomp;  i++) {
    ff[i] = REAL(resultsxp)[i] * prdbounds;
 }
} // End RIntegrand


/*********************************************************************/
 void cuhreDoSample(count n, ctreal *x,  real *f,
		    SEXP rho, SEXP globf, Glob *globdim)
{

/* rho, globf: La fonction R à intégrer et son environnement d'execution */
		
  globdim->neval_ += n;

  while( n-- ) {
    RIntegrand(&(globdim->ndim_), x, &(globdim->ncomp_),  globdim->lower_, globdim->upper_, globdim->prdbounds_, f,
	       rho, globf, globdim);
    x += globdim->ndim_;
    f += globdim->ncomp_;
  }
}
