#ifndef __cuhre_common_h__
#define __cuhre_common_h__

// Compilation note for R interface: move common.c into cuhre_common.h
/*
	common.c
		includes most of the modules
		this file is part of Cuhre
		last modified 14 Feb 05 th
*/


#include "common_ChiSquare.h"
#include "cuhre_Rule.h"
#include "cuhre_Integrate.h"


static inline bool BadDimension(ccount ndim)
{
#if NDIM > 0
  if( ndim > NDIM ) return true;
#endif
  return ndim < 2;
}


static inline bool BadComponent(cint ncomp)
{
#if NCOMP > 0
  if( ncomp > NCOMP ) return true;
#endif
  return ncomp < 1;
}

#endif
