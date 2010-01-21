/*
	Divonne.c
		Multidimensional integration by partitioning
		originally by J.H. Friedman and M.H. Wright
		(CERNLIB subroutine D151)
		this version by Thomas Hahn
		last modified 2 Mar 06 th
*/

#include "common_stddecl.h"
#include "divonne_util.h"
#include "struct_Random.h"

#include "inclR.h"
extern bool divonneBadDimension(ccount ndim, cint flags, ccount key);
  extern bool divonneBadComponent(cint ncomp);
extern void divonneDoSample(number n, ccount ldx, ctreal *x,real *f);
extern int divonneIntegrate(ctreal epsrel, ctreal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  int key1, int key2, int key3, ccount maxpass, 
  ctreal maxchisq, ctreal mindeviation,
			    real *integral, real *erreur, real *prob);

/*********************************************************************/
/* Compilation note for R interface: 
   The integration R function and its execution environnement */
/*********************************************************************/
SEXP rho, globf, peakf;		


/* verif_: to verify the dimensions of the structures returned by
   the user functions: only used at the first call */
bool verif_; 

/*********************************************************************/
/*  The function RIntegrand calls the R user function */
/*********************************************************************/
static void RIntegrand(ccount *ndim, ctreal xx[],
		       ccount  *ncomp, 
		       ctreal *lower, ctreal *upper, ctreal prdbounds,
		       real ff[],
		       cint *phase)
{
  SEXP args, argphase, s, t, resultsxp;
  int i;
 /*f:  the R function and its environment, rho are global */
  // The input arguments are x +1 phase 
   PROTECT(args=allocVector(REALSXP, ( *ndim )));
PROTECT(argphase=allocVector(REALSXP, (1 )));
// The output  arguments are ncomp values
  PROTECT(resultsxp=allocVector(REALSXP, ( *ncomp )));
 /* Fill in the input arguments with rescaling between   0-1,
     according to the bounds */
  for (i =0; i<*ndim; i++) 
    REAL(args)[i] = xx[i] * (upper_[i] - lower_[i]) + lower_[i];
  REAL(argphase)[ 0] = *phase;

  /* Call the R function */
PROTECT(t = s = allocList(3));
         SET_TYPEOF(s, LANGSXP);
         SETCAR(t, globf); t = CDR(t);
         SETCAR(t,  args); t = CDR(t);
         SETCAR(t, argphase); 

 PROTECT(resultsxp=eval(s,rho));
 UNPROTECT(5);


 if (verif_== true) {
 if  (length(resultsxp) != *ncomp)
   error("Function integrand does not return a vector of length ncomp");
 verif_= false; // do not verify the next time
 }


 for (i =0; i<*ncomp;  i++) {
    ff[i] = REAL(resultsxp)[i] * prdbounds;
 }

} // End RIntegrand

/*********************************************************************/
/* The function Rpeakf calls the R user function */
/*********************************************************************/
static void Rpeakf(const int *ndim, 
		   const double b[],
		    int *n, double x[])
{
  SEXP args, callsxp,resultsxp, dim;
  int i,j, p,l,nl, kl, nx;

 /*peakf:  the R function and its environment, rho are global */
  PROTECT(args=allocVector(REALSXP, ( 2*(*ndim) )));
  /*  n is the maximum number of points in output,
i.e the argument nextra of the function divonne */
  PROTECT(resultsxp=allocVector(REALSXP, (1+ ((*n) * (*ndim)))));
  /* Rescaler les bornes inf et sup qui sont dans l'hypercube unité
     selon l'échelle de l'utilisateur */
  /* Les binf et bsup: les mettre dans une matrice R (ndim,2) */

 kl=0;
  j=0;
  for (i =0; i< *ndim; i++) {
    l=kl;
    kl++;
  for (p=0; p<2; p++) {
    REAL(args)[l] = b[j++] * (upper_[i] - lower_[i]) + lower_[i];
    l=l+ *ndim;
  }
  } // fin i

  /* Put the bounds into a R matrix */
PROTECT(dim = allocVector(INTSXP, 2));
       INTEGER(dim)[0] = *ndim; INTEGER(dim)[1] = 2;
       setAttrib(args, R_DimSymbol, dim);
     
       /* Affect labels to the columns */
       SEXP dimnames, elmt;
       // Create the list dimnames
PROTECT(dimnames = allocVector(VECSXP, 2));
// Create a vector of length 2 to store the columns labels
 PROTECT( elmt = allocVector(STRSXP, 2));
 // Affect the labels
 SET_STRING_ELT( elmt,0, mkChar("lower"));
 SET_STRING_ELT( elmt,1, mkChar("upper"));
 // Affect the vector to the 2nd component of dimnames
 SET_VECTOR_ELT(dimnames, 1, elmt);
 // Affect  dimnames to the argument
 setAttrib(args, R_DimNamesSymbol, dimnames);

 /* Call the R function */
  PROTECT(callsxp=lang2( peakf, args));
 PROTECT(resultsxp=eval(callsxp,rho));
 UNPROTECT(7);

 nx=INTEGER(getAttrib(resultsxp, R_DimSymbol))[1];


 /* Verify the dimensions of the result */
  if (verif_== true) {
 nl=INTEGER(getAttrib(resultsxp, R_DimSymbol))[0];

 if (nl != *ndim)
   error("peakfinder does not return a matrix with ndim rows.");


 // In output, n is the effective number of pointd
 // Verify it is <= nextra
 if (nx > (*n))
   error("peakfinder returns more than nextra points");
  } // fin verif

 *n =nx;
 j=0;
   /* Rescale the points into the unit hypercube */
 /* In R x is a matrix (ndim, npoints),
    so, the values are stored point by point */
   for (p=0; p< *n; p++) {
 for (i =0; i< *ndim; i++) {
  switch(TYPEOF(resultsxp)) {
         case REALSXP:
     x[j] = (REAL(resultsxp)[j] - lower_[i]) / (upper_[i] - lower_[i]);
     break;
  case INTSXP:
     x[j] = (INTEGER(resultsxp)[j] - lower_[i]) / (upper_[i] - lower_[i]);
     break;
  default:
    error("peakfinder does not return a real or integer matrix");
  } // fin switch

     j++;
   } // fin i
 } // fin p

} // fin Rpeakf



/*********************************************************************/


 void EXPORT(Divonne)(ccount ndim, ccount ncomp,
  Integrand integrand,
  ctreal epsrel, ctreal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  cint key1, cint key2, cint key3, ccount maxpass,
  ctreal border, ctreal maxchisq, ctreal mindeviation,
  cnumber ngiven, ccount ldxgiven, real *xgiven,
  cnumber nextra, PeakFinder peakfinder,
  int *pnregions, number *pneval, int *pfail,
  real *integral, real *erreur, real *prob)
{

  ndim_ = ndim;
  ncomp_ = ncomp;

  if( divonneBadComponent(ncomp) ||
      divonneBadDimension(ndim, flags, key1) ||
      divonneBadDimension(ndim, flags, key2) ||
      ((key3 & -2) && divonneBadDimension(ndim, flags, key3)) ) *pfail = -1;
  else {
    neval_ = neval_opt_ = neval_cut_ = 0;
    integrand_ = integrand;
    peakfinder_ = peakfinder;
    border_.lower = border;
    border_.upper = 1 - border_.lower;
    ngiven_ = ngiven;
    xgiven_ = NULL;
    ldxgiven_ = IMax(ldxgiven, ndim_);
    nextra_ = nextra;

    if( ngiven + nextra ) {
      cnumber nxgiven = ngiven*ldxgiven;
      cnumber nxextra = nextra*ldxgiven;
      cnumber nfgiven = ngiven*ncomp;
      cnumber nfextra = nextra*ncomp;

      Alloc(xgiven_, nxgiven + nxextra + nfgiven + nfextra);
      xextra_ = xgiven_ + nxgiven;
      fgiven_ = xextra_ + nxextra;
      fextra_ = fgiven_ + nfgiven;

      if( nxgiven ) {
        phase_ = 0;
        Copy(xgiven_, xgiven, nxgiven);
        divonneDoSample(ngiven_, ldxgiven_, xgiven_, fgiven_);
      }
    }

    *pfail = divonneIntegrate( epsrel, Max(epsabs, NOTZERO),
      flags, mineval, maxeval, key1, key2, key3, maxpass,
      maxchisq, mindeviation,
      integral, erreur, prob);

    *pnregions = nregions_;
    *pneval = neval_;

    if( xgiven_ ) free(xgiven_);
  }
}

/*********************************************************************/
void (divonne)(ccount *pndim, ccount *pncomp,
  Integrand integrand,
 ctreal *pepsrel, ctreal *pepsabs,
  cint *pflags, cnumber *pmineval, cnumber *pmaxeval,
  cint *pkey1, cint *pkey2, cint *pkey3, ccount *pmaxpass,
  ctreal *pborder, ctreal *pmaxchisq, ctreal *pmindeviation,
  cnumber *pngiven, ccount *pldxgiven, real *xgiven,
  cnumber *pnextra, PeakFinder peakfinder,
  int *pnregions, number *pneval, int *pfail,
  real *integral, real *erreur, real *prob)
{

  EXPORT(Divonne)(*pndim, *pncomp,
    integrand,
    *pepsrel, *pepsabs,
    *pflags, *pmineval, *pmaxeval,
    *pkey1, *pkey2, *pkey3, *pmaxpass,
    *pborder, *pmaxchisq, *pmindeviation,
    *pngiven, *pldxgiven, xgiven,
    *pnextra, peakfinder,
    pnregions, pneval, pfail,
    integral, erreur, prob);
}

/*********************************************************************/
/*  Compilation note for R interface: add this subroutine
  interface between R  and divonne */
/*********************************************************************/


void Rdivonne(int *pndim, int *pncomp,
	     void ** integrand, void *env,
    real *lower, real *upper, real *prdbounds,
double *pepsrel, double *pepsabs,
   int *pmersenneseed, int *pflags, int *pmineval, int *pmaxeval,
	      int *key1, int *key2, int *key3,
	      int *maxpass, double *border,
double *pmaxchisq, double *pmindeviation,
  int *pngiven, int *pldxgiven, double *xgiven,
	      int *pnextra, void ** peakfinder,
  int *pnregions, int *pneval, int *pfail,
  double *integral, double *erreur, double *prob)
{


  /* store the R function and its environment  into a global*/
  rho= (SEXP)(env);
  globf= (SEXP)(integrand);
  peakf = (SEXP)peakfinder;
  lower_ = lower;
  upper_ = upper;
  prdbounds_ = *prdbounds;

   if (NA_INTEGER != *pmersenneseed) 
    SUFFIX(mersenneseed)= *pmersenneseed;

    verif_= true;

  /* call divonne */

  divonne((ccount *)pndim, (ccount *)pncomp,
	(Integrand)RIntegrand,
(ctreal *)pepsrel, (ctreal *)pepsabs,
	(cint *)pflags, (cnumber *)pmineval, (cnumber *)pmaxeval,
 (cint *)key1, (cint *)key2, (cint *)key3, 
(ccount *)maxpass,
  (ctreal *)border, (ctreal *)pmaxchisq, (ctreal *)pmindeviation,
	  (cnumber *)pngiven, (ccount *)pldxgiven, (double *)xgiven,
	  (cnumber *)pnextra, (PeakFinder) Rpeakf,
	(count *) pnregions,
	(number *)pneval, (int *)pfail,
	(real *)integral, (real *)erreur, (real *)prob);

} // End Rdivonne
