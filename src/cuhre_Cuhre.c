/*
	Cuhre.c
		Adaptive integration using cubature rules
		by Thomas Hahn
		last modified 2 Mar 06 th
*/


#include "cuhre_util.h"
#include "inclR.h"

extern bool cuhreBadDimension(ccount ndim);
extern bool cuhreBadComponent(cint ncomp);
extern int cuhreIntegrate(  ctreal epsrel, ctreal epsabs,
  cint flags, number mineval, cnumber maxeval, ccount key,
			    real *integral, real *erreur, real *prob,
			    SEXP rho, SEXP globf, Glob *globdim);


/*********************************************************************/

 void EXPORT(Cuhre)(ccount ndim, ccount ncomp,
  ctreal epsrel, ctreal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  ccount key,
  count *pnregions, number *pneval, int *pfail,
  real *integral, real *erreur, real *prob,
		    SEXP rho, SEXP globf, Glob *globdim)
{

/* rho, globf: La fonction R à intégrer et son environnement d'execution */		
  globdim->ndim_ = ndim;
  globdim->ncomp_ = ncomp;

  if( cuhreBadComponent(ncomp) || cuhreBadDimension(ndim) ) *pfail = -1;
  else {
    globdim->neval_ = 0;

    *pfail = cuhreIntegrate( epsrel, Max(epsabs, NOTZERO),
      flags, mineval, maxeval, key,
			     integral, erreur, prob, 
			     rho, globf, globdim);
    *pnregions = globdim->nregions_;
    *pneval =globdim->neval_;
  }
}

/*********************************************************************/

void (cuhre)(ccount *pndim, ccount *pncomp,
  ctreal *pepsrel, ctreal *pepsabs,
  cint *pflags, cnumber *pmineval, cnumber *pmaxeval,
  ccount *pkey,
  count *pnregions, number *pneval, int *pfail,
	     real *integral, real *erreur, real *prob,
	     SEXP rho, SEXP globf, Glob *globdim)
{

/* rho, globf: La fonction R à intégrer et son environnement d'execution */		

  EXPORT(Cuhre)(*pndim, *pncomp,
   *pepsrel, *pepsabs,
    *pflags, *pmineval, *pmaxeval,
    *pkey,
    pnregions, pneval, pfail,
		integral, erreur, prob, rho, globf, globdim);
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
  double *integral, double *erreur, double *prob,
  SEXP rho, SEXP globf)
{

/* rho, globf: La fonction R à intégrer et son environnement d'execution */		

  /* store the R function and its environment */
  rho= (SEXP)(env);
  globf= (SEXP)(integrand);
  Glob globdim;

  globdim.lower_ = lower;
  globdim.upper_ = upper;
  globdim.prdbounds_ = *prdbounds;
  /* mersenneseed is ignored in cuhre

  if (NA_INTEGER != *pmersenneseed) 
    SUFFIX(mersenneseed)= *pmersenneseed;
 */

 /* call cuhre */
  cuhre((ccount *)pndim, (ccount *)pncomp,
 (ctreal *)pepsrel, (ctreal *)pepsabs,
	(cint *)pflags, (cnumber *)pmineval, (cnumber *)pmaxeval,
	(ccount *)key,
	(count *) pnregions,
	(number *)pneval, (int *)pfail,
	(real *)integral, (real *)erreur, (real *)prob, 
	rho, globf, &globdim);

} // End Rcuhre
