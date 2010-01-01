//Compilation note for R interface: move into a .h

/*
	common.c
		includes most of the modules
		this file is part of Suave
		last modified 14 Feb 05 th
*/


#include "common_Random.h"
#include "common_ChiSquare.c"
#include "suave_Grid.h"
#include "suave_Sample.h"
#include "suave_Fluct.h"
#include "suave_Integrate.h"


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

