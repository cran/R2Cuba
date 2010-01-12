#ifndef __divonne_common_h__
#define __divonne_common_h__

// Compilation note for R interface:  move common.c into divonne_common.h

/*
	common.c
		includes most of the modules
		this file is part of Divonne
		last modified 5 May 09 th
*/


static bool Explore(count iregion, cSamples *samples, cint depth, cint flags);

static void Split(count iregion, int depth);


#include "common_Random.h"
#include "common_ChiSquare.h"
#include "divonne_Rule.h"
#include "divonne_Sample.h"
#include "divonne_FindMinimum.h"
#include "divonne_Explore.h"
#include "divonne_Split.h"
#include "divonne_Integrate.h"


static inline bool BadDimension(ccount ndim, cint flags, ccount key)
{
#if NDIM > 0
  if( ndim > NDIM ) return true;
#endif
  if( IsSobol(key) ) return
    ndim < SOBOL_MINDIM || (!PSEUDORNG && ndim > SOBOL_MAXDIM);
  if( IsRule(key, ndim) ) return ndim < 1;
  return ndim < KOROBOV_MINDIM || ndim > KOROBOV_MAXDIM;
}


static inline bool BadComponent(cint ncomp)
{
#if NCOMP > 0
  if( ncomp > NCOMP ) return true;
#endif
  return ncomp < 1;
}

#endif
