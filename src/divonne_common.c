/*
	common.c
		includes most of the modules
		this file is part of Divonne
		last modified 5 May 09 th
*/

#include "common_stddecl.h"
#include "divonne_util.h"
#include "divonne_KorobovCoeff.h"

#define IsRule(k, d) (k == 9 || k == 7 || (k == 11 && d == 3) || (k == 13 && d == 2))



 bool divonneBadDimension(ccount ndim, cint flags, ccount key)
{
#if NDIM > 0
  if( ndim > NDIM ) return true;
#endif
  if( IsSobol(key) ) return
    ndim < SOBOL_MINDIM || (!PSEUDORNG && ndim > SOBOL_MAXDIM);
  if( IsRule(key, ndim) ) return ndim < 1;
  return ndim < KOROBOV_MINDIM || ndim > KOROBOV_MAXDIM;
}


 bool divonneBadComponent(cint ncomp)
{
#if NCOMP > 0
  if( ncomp > NCOMP ) return true;
#endif
  return ncomp < 1;
}

