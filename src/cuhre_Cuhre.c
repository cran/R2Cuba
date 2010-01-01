/*
	Cuhre.c
		Adaptive integration using cubature rules
		by Thomas Hahn
		last modified 2 Mar 06 th
*/


#include "cuhre_util.c"
/* Compilation note for R interface: add inclR.h */
#include "inclR.h"

/*  Compilation note for R interface:
 modif #define Print(s) puts(s); fflush(stdout) */
#define Print(s) Rprintf(s)


static Integrand integrand_;

/*********************************************************************/
/* La fonction R à intégrer et son environnement d'execution */
/*********************************************************************/
SEXP rho, f;		
/* Compilation note for R interface:  global, to be used by Sample*/
real *lower_, *upper_,  prdbounds_;

/*********************************************************************/
/* The function RIntegrand calls the R user function */
/*********************************************************************/
static void RIntegrand(ccount *ndim, creal xx[],
		       ccount  *ncomp,
		       creal *lower, creal *upper, creal prdbounds, 
		       real ff[])
{
  SEXP args, callsxp,resultsxp;
  int i;
 /*f:  the R function and its environment, rho are global */
   // The input arguments are the x
  PROTECT(args=allocVector(REALSXP, ( *ndim )));
  PROTECT(resultsxp=allocVector(REALSXP, ( *ncomp )));
  /* Fill in the input arguments with rescaling between   0-1,
     according to the bounds */
  for (i =0; i<*ndim; i++) 
    REAL(args)[i] = xx[i] * (upper_[i] - lower_[i]) + lower_[i];

  /* Call the R function */
  PROTECT(callsxp=lang2( f, args));
  // PROTECT(args=eval(callsxp,rho));
 PROTECT(resultsxp=eval(callsxp,rho));
 UNPROTECT(4);
if  (length(resultsxp) != *ncomp)
   error("Function integrand does not return a vector of length ncomp");

 for (i =0; i<*ncomp;  i++) {
    ff[i] = REAL(resultsxp)[i] * prdbounds;
 }
} // End RIntegrand

/*********************************************************************/

static inline void DoSample(count n, creal *x,  real *f)
{
  neval_ += n;
  while( n-- ) {
    integrand_(&ndim_, x, &ncomp_,  lower_, upper_, prdbounds_, f);
    x += ndim_;
    f += ncomp_;
  }
}

/*********************************************************************/

#include "cuhre_common.h"

Extern void EXPORT(Cuhre)(ccount ndim, ccount ncomp,
  Integrand integrand,
  creal epsrel, creal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  ccount key,
  count *pnregions, number *pneval, int *pfail,
  real *integral, real *error, real *prob)
{
  ndim_ = ndim;
  ncomp_ = ncomp;

  if( BadComponent(ncomp) || BadDimension(ndim) ) *pfail = -1;
  else {
    neval_ = 0;
    integrand_ = integrand;

    *pfail = Integrate( epsrel, Max(epsabs, NOTZERO),
      flags, mineval, maxeval, key,
      integral, error, prob);

    *pnregions = nregions_;
    *pneval = neval_;
  }
}

/*********************************************************************/

void (cuhre)(ccount *pndim, ccount *pncomp,
  Integrand integrand,
  creal *pepsrel, creal *pepsabs,
  cint *pflags, cnumber *pmineval, cnumber *pmaxeval,
  ccount *pkey,
  count *pnregions, number *pneval, int *pfail,
  real *integral, real *error, real *prob)
{

  EXPORT(Cuhre)(*pndim, *pncomp, integrand,
   *pepsrel, *pepsabs,
    *pflags, *pmineval, *pmaxeval,
    *pkey,
    pnregions, pneval, pfail,
    integral, error, prob);
}

/*********************************************************************/
/* Compilation note for R interface:  
add this interface subroutine between R and cuhre */
/*********************************************************************/


void Rcuhre(int *pndim, int *pncomp,
	     void ** integrand, void *env,
   real *lower, real *upper, real *prdbounds,
   double *pepsrel, double *pepsabs,
  int *pmersenneseed, int *pflags, int *pmineval, int *pmaxeval,
	    int *key, 
  int *pnregions, int *pneval, int *pfail,
  double *integral, double *error, double *prob)
{
extern unsigned int mersenneseed;

  /* store the R function and its environment  in a  global*/
  rho= (SEXP)(env);
  f= (SEXP)(integrand);
  lower_ = lower;
  upper_ = upper;
  prdbounds_ = *prdbounds;
  /* mersenneseed is ignored in cuhre */

  if (NA_INTEGER != *pmersenneseed) 
    SUFFIX(mersenneseed)= *pmersenneseed;

 /* call cuhre */

  cuhre((ccount *)pndim, (ccount *)pncomp,
	(Integrand)RIntegrand,
 (creal *)pepsrel, (creal *)pepsabs,
	(cint *)pflags, (cnumber *)pmineval, (cnumber *)pmaxeval,
	(ccount *)key,
	(count *) pnregions,
	(number *)pneval, (int *)pfail,
	(real *)integral, (real *)error, (real *)prob);

} // End Rcuhre
