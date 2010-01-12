#ifndef __vegas_common_h__
#define __vegas_common_h__
//Compilation note for R interface: move into a .h
/*
	common.c
		include most of the modules
		this file is part of Vegas
		last modified 14 Feb 05 th
*/


#include "common_Random.h"
#include "common_ChiSquare.h"
#include "vegas_Grid.h"
#include "vegas_Integrate.h"


static inline bool BadDimension(cint ndim, cint flags)
{
#if NDIM > 0
  if( ndim > NDIM ) return true;
#endif
  return ndim < SOBOL_MINDIM || (!PSEUDORNG && ndim > SOBOL_MAXDIM);
}


static inline bool BadComponent(cint ncomp)
{
#if NCOMP > 0
  if( ncomp > NCOMP ) return true;
#endif
  return ncomp < 1;
}
#endif
