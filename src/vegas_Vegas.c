/*
	Vegas.c
		Vegas Monte-Carlo integration
		by Thomas Hahn
		last modified 30 Aug 07 th
*/


#include "struct_Random.h"
#include "vegas_util.h"
//Compilation note for R interface: add inclR.h 
#include "inclR.h"
extern bool vegasBadDimension(cint ndim, cint flags);
extern bool vegasBadComponent(cint ncomp);

extern  int vegasIntegrate(  ctreal *lower, ctreal *upper, ctreal prdbounds,
		       ctreal epsrel, ctreal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  cnumber nstart, cnumber nincrease,
		      real *integral, real *erreur, real *prob);

/* Compilation note for R interface: remove
 int EXPORT(vegasnbatch) = 1000;
 int EXPORT(vegasgridno) = 0;*/
 int EXPORT(vegasnbatch);
 int EXPORT(vegasgridno);
 char EXPORT(vegasstate)[MAXSTATESIZE] = "";

/* Compilation note for R interface: modif #define Print(s) puts(s); fflush(stdout) */
#define Print(s) Rprintf(s)

/*********************************************************************/
/* Compilation note for R interface: 
   The integration R function and its execution environnement */
/*********************************************************************/
SEXP rho, globf;		

/*********************************************************************/
/*  The function RIntegrand calls the R user function */
/*********************************************************************/
static void RIntegrand(ccount *ndim, ctreal  xx[],
		       ccount *ncomp,
		       ctreal *lower, ctreal *upper, ctreal prdbounds,
		        real ff[], ctreal *weight)
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

  /* Call the R function */
 PROTECT(t = s = allocList(3));
         SET_TYPEOF(s, LANGSXP);
         SETCAR(t, globf); t = CDR(t);
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


 void EXPORT(Vegas)(ccount ndim, ccount ncomp,
  Integrand integrand,
  ctreal *lower, ctreal *upper, ctreal prdbounds,
  ctreal epsrel, ctreal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  cnumber nstart, cnumber nincrease,
  number *pneval, int *pfail,
  real *integral, real *erreur, real *prob)
{
  ndim_ = ndim;
  ncomp_ = ncomp;
R_CheckUserInterrupt(); // to allow user interruptions

  if( vegasBadComponent(ncomp) || vegasBadDimension(ndim, flags) ) *pfail = -1;
  else {
    neval_ = 0;
    integrand_ = integrand;

    *pfail = vegasIntegrate(lower, upper, prdbounds,
		       epsrel, epsabs,
      flags, mineval, maxeval, nstart, nincrease,
      integral, erreur, prob);

    *pneval = neval_;
  }
}

/*********************************************************************/

 void (vegas)(ccount *pndim, ccount *pncomp,
  Integrand integrand,
  ctreal *lower, ctreal *upper, ctreal *prdbounds,
  ctreal *pepsrel, ctreal *pepsabs,
  cint *pflags, cnumber *pmineval, cnumber *pmaxeval,
  cnumber *pnstart, cnumber *pnincrease, 
  number *pneval, int *pfail,
  real *integral, real *erreur, real *prob)
{
  /* make sure the filename is null-terminated */
  if( *EXPORT(vegasstate) ) {
    char *p;
    EXPORT(vegasstate)[sizeof(EXPORT(vegasstate)) - 1] = 0;
    if( (p = strchr(EXPORT(vegasstate), ' ')) ) *p = 0;
  }

  EXPORT(Vegas)(*pndim, *pncomp,
		integrand,
		lower, upper, *prdbounds,
		*pepsrel, *pepsabs,
    *pflags, *pmineval, *pmaxeval,
    *pnstart, *pnincrease,
    pneval, pfail,
    integral, erreur, prob);
}

/*********************************************************************/
/* Rvegas :  interface between R and vegas */
/*********************************************************************/


void Rvegas(int *pndim, int *pncomp,
	     void ** integrand, void *env,
	    real *lower, real *upper, real *prdbounds,
	    double *pepsrel, double *pepsabs,
	    int *pmersenneseed, int *pvegasnbatch,
	    int *pvegasgridno, 
	    int *pflags, int *pmineval, int *pmaxeval,
	    int *pnstart, int *pnincrease, 
	    char **state,
	    int *pneval, int *pfail,
	    double *integral, double *erreur, double *prob)
{


  /* store the R function and its environment in a global*/
  rho= (SEXP)(env);
  globf= (SEXP)(integrand);

  if (strlen(*state) >0) {
    strncpy(vegasstate, *state, 128);
   }


  if (NA_INTEGER != *pmersenneseed) 
    SUFFIX(mersenneseed) = *pmersenneseed;
  if (NA_INTEGER != *pvegasnbatch) 
    EXPORT(vegasnbatch)= *pvegasnbatch;
  else
    EXPORT(vegasnbatch)=1000;

  if (NA_INTEGER != *pvegasgridno) 
    EXPORT(vegasgridno)= *pvegasgridno;
  else
    EXPORT(vegasgridno)=0;

  /* call vegas */
  vegas((ccount *)pndim, (ccount *)pncomp,
	(Integrand)RIntegrand,
(ctreal *)lower, (ctreal *)upper, (ctreal *)prdbounds,
  (ctreal *)pepsrel, (ctreal *)pepsabs,
	(cint *)pflags, (cnumber *)pmineval, (cnumber *)pmaxeval,
  (cnumber *)pnstart, (cnumber *)pnincrease, 
	(number *)pneval, (int *)pfail,
	(real *)integral, (real *)erreur, (real *)prob);

} // End Rvegas
