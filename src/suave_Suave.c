/*
	Suave.c
		Subregion-adaptive Vegas Monte-Carlo integration
		by Thomas Hahn
		last modified 30 Aug 07 th
*/


#include "suave_util.h"
/* Compilation note for R interface: add inclR.h  */
#include "inclR.h"

/* Compilation note for R interface: modif #define Print(s) puts(s); fflush(stdout) */
#define Print(s) Rprintf(s)

static Integrand integrand_;
/*********************************************************************/
/* Compilation note for R interface: 
   The integration R function and its execution environnement */
/*********************************************************************/
SEXP rho, f;		
/*Compilation note for R interface: global, to be used by Sample */
real *lower_, *upper_,  prdbounds_;

/*********************************************************************/
/*The function RIntegrand calls the R user function */
/*********************************************************************/
static void RIntegrand(ccount *ndim, creal xx[],
                      ccount *ncomp, 
       creal *lower, creal *upper, creal prdbounds, real ff[],
		       creal *weight)
{
   SEXP args, argw, s, t, resultsxp;
  int i;

 /*f:  the R function and its environment, rho are global */
  // The input arguments are x +1 weight 
  PROTECT(args=allocVector(REALSXP, ( *ndim )));
PROTECT(argw=allocVector(REALSXP, (1 )));
  PROTECT(resultsxp=allocVector(REALSXP, ( *ncomp )));
  /* Fill in the input arguments with rescaling between   0-1,
     according to the bounds */
  for (i =0; i<*ndim; i++) 
    REAL(args)[i] = xx[i] * (upper[i] - lower[i]) + lower[i];
  REAL(argw)[ 0]=*weight; 

  /* Appel de la fonction R */
 PROTECT(t = s = allocList(3));
         SET_TYPEOF(s, LANGSXP);
         SETCAR(t, f); t = CDR(t);
         SETCAR(t,  args); t = CDR(t);
         SETCAR(t, argw);

 PROTECT(resultsxp=eval(s,rho));
 UNPROTECT(5);
if  (length(resultsxp) != *ncomp)
   error("Function integrand does not return a vector of length ncomp");


 for (i =0; i<*ncomp;  i++) {
   ff[i] = REAL(resultsxp)[i] * prdbounds;
 }
} // End RIntegrand

/*********************************************************************/

static inline void DoSample(number n, creal *w, creal *x, real *f)
{
  neval_ += n;
  while( n-- ) {
    integrand_(&ndim_, x, &ncomp_, lower_, upper_, prdbounds_, f, w++);
    x += ndim_;
    f += ncomp_;
  }
}

/*********************************************************************/

#include "suave_common.h"

Extern void EXPORT(Suave)(ccount ndim, ccount ncomp,
  Integrand integrand,
  creal epsrel, creal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  cnumber nnew, creal flatness,
  count *pnregions, number *pneval, int *pfail,
  real *integral, real *error, real *prob)
{
  ndim_ = ndim;
  ncomp_ = ncomp;

  if( BadComponent(ncomp) || BadDimension(ndim, flags) ) *pfail = -1;
  else {
    neval_ = 0;
    integrand_ = integrand;

    *pfail = Integrate( epsrel, Max(epsabs, NOTZERO),
      flags, mineval, maxeval, nnew, flatness,
      integral, error, prob);
    *pnregions = nregions_;
    *pneval = neval_;
  }
}

/*********************************************************************/

void (suave)(ccount *pndim, ccount *pncomp,
  Integrand integrand,
  creal *pepsrel, creal *pepsabs,
  cint *pflags, cnumber *pmineval, cnumber *pmaxeval,
  cnumber *pnnew, creal *pflatness,
  count *pnregions, number *pneval, int *pfail,
  real *integral, real *error, real *prob)
{
  EXPORT(Suave)(*pndim, *pncomp,
    integrand,
    *pepsrel, *pepsabs,
    *pflags, *pmineval, *pmaxeval,
    *pnnew, *pflatness,
    pnregions, pneval, pfail,
    integral, error, prob);
}

/*********************************************************************/
/* Add Rsuave: interface between R and suave */
/*********************************************************************/


void Rsuave(int *pndim, int *pncomp,
	     void ** integrand, void *env,
  real *lower, real *upper, real *prdbounds,
  double *pepsrel, double *pepsabs,
   int *pmersenneseed, int *pflags, int *pmineval, int *pmaxeval,
  int *pnnew, double *pflatness, 
  int *pnregions, int *pneval, int *pfail,
  double *integral, double *error, double *prob)
{
  /* store the R function and its environment  dans un global*/
  rho= (SEXP)(env);
  f= (SEXP)(integrand);
  lower_ = lower;
  upper_ = upper;
  prdbounds_ = *prdbounds;

   if (NA_INTEGER != *pmersenneseed) 
     SUFFIX(mersenneseed)= *pmersenneseed;

 /* call suave */

  suave((ccount *)pndim, (ccount *)pncomp,
	(Integrand)RIntegrand,
  (creal *)pepsrel, (creal *)pepsabs,
	(cint *)pflags, (cnumber *)pmineval, (cnumber *)pmaxeval,
  (cnumber *)pnnew, (creal *)pflatness, 
	(count *) pnregions,
	(number *)pneval, (int *)pfail,
	(real *)integral, (real *)error, (real *)prob);

} // End Rsuave
