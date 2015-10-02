/*
	Suave.c
		Subregion-adaptive Vegas Monte-Carlo integration
		by Thomas Hahn
		last modified 30 Aug 07 th
*/


#include "common_stddecl.h"
#include "struct_Random.h"
#include "suave_util.h"
#include "inclR.h"
extern bool suaveBadDimension(cint ndim, cint flags);
extern bool suaveBadComponent(cint ncomp);
extern  int suaveIntegrate(  ctreal epsrel, ctreal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  cnumber nnew, ctreal flatness,
		      real *integral, real *erreur, real *prob);

/* Compilation note for R interface: modif #define Print(s) puts(s); fflush(stdout) */
#define Print(s) Rprintf(s)

/*********************************************************************/
/* Compilation note for R interface: 
   The integration R function and its execution environnement */
/*********************************************************************/
SEXP rho, globf;		

/*********************************************************************/
/*The function RIntegrand calls the R user function */
/*********************************************************************/
static void RIntegrand(ccount *ndim, ctreal xx[],
                      ccount *ncomp, 
       ctreal *lower, ctreal *upper, ctreal prdbounds, real ff[],
		       ctreal *weight)
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
         SETCAR(t, globf); t = CDR(t);
         SETCAR(t,  args); t = CDR(t);
         SETCAR(t, argw);

 PROTECT(resultsxp=eval(s,rho));
 UNPROTECT(5);
if  (length(resultsxp) != *ncomp)
  error("Function integrand does not return a vector of length ncomp\n Length of returned vector= %d. ncomp=%d\n",
	length(resultsxp), *ncomp);


 for (i =0; i<*ncomp;  i++) {
   ff[i] = REAL(resultsxp)[i] * prdbounds;
 }
} // End RIntegrand


/*********************************************************************/


 void EXPORT(Suave)(ccount ndim, ccount ncomp,
  Integrand integrand,
  ctreal epsrel, ctreal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  cnumber nnew, ctreal flatness,
  count *pnregions, number *pneval, int *pfail,
  real *integral, real *erreur, real *prob)
{
  ndim_ = ndim;
  ncomp_ = ncomp;

  if( suaveBadComponent(ncomp) || suaveBadDimension(ndim, flags) ) *pfail = -1;
  else {
    neval_ = 0;
    integrand_ = integrand;

    *pfail = suaveIntegrate( epsrel, Max(epsabs, NOTZERO),
      flags, mineval, maxeval, nnew, flatness,
      integral, erreur, prob);
    *pnregions = nregions_;
    *pneval = neval_;
  }
}

/*********************************************************************/

void (suave)(ccount *pndim, ccount *pncomp,
  Integrand integrand,
  ctreal *pepsrel, ctreal *pepsabs,
  cint *pflags, cnumber *pmineval, cnumber *pmaxeval,
  cnumber *pnnew, ctreal *pflatness,
  count *pnregions, number *pneval, int *pfail,
  real *integral, real *erreur, real *prob)
{
  EXPORT(Suave)(*pndim, *pncomp,
    integrand,
    *pepsrel, *pepsabs,
    *pflags, *pmineval, *pmaxeval,
    *pnnew, *pflatness,
    pnregions, pneval, pfail,
    integral, erreur, prob);
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
  double *integral, double *erreur, double *prob)
{
  /* store the R function and its environment  dans un global*/
  rho= (SEXP)(env);
  globf= (SEXP)(integrand);
  lower_ = lower;
  upper_ = upper;
  prdbounds_ = *prdbounds;

   if (NA_INTEGER != *pmersenneseed) 
     SUFFIX(mersenneseed)= *pmersenneseed;

 /* call suave */
  suave((ccount *)pndim, (ccount *)pncomp,
	(Integrand)RIntegrand,
  (ctreal *)pepsrel, (ctreal *)pepsabs,
	(cint *)pflags, (cnumber *)pmineval, (cnumber *)pmaxeval,
  (cnumber *)pnnew, (ctreal *)pflatness, 
	(count *) pnregions,
	(number *)pneval, (int *)pfail,
	(real *)integral, (real *)erreur, (real *)prob);

} // End Rsuave
