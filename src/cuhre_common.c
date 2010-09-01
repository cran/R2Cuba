/*
	common.c
		includes most of the modules
		this file is part of Cuhre
		last modified 14 Feb 05 th
*/

#include "common_stddecl.h"
#include "cuhre_util.h"
#include "inclR.h"


 bool cuhreBadDimension(ccount ndim)
{
#if NDIM > 0
  if( ndim > NDIM ) return true;
#endif
return ndim < 2;
}


 bool cuhreBadComponent(cint ncomp)
{
#if NCOMP > 0
  if( ncomp > NCOMP ) return true;
#endif
  return ncomp < 1;
}


