#ifndef __cuhre_util_h__
#define __cuhre_util_h__
/*
	util.c
		Utility functions
		this file is part of Cuhre
		last modified 9 Jan 05 th
*/


#include "cuhre_decl.h"
/* globals */
 count ndim_, ncomp_, nregions_;
 number neval_;
 real *lower_, *upper_,  prdbounds_;
 Integrand integrand_;
enum { nrules = 5 };

#define CUHRETYPEDEFSET \
  typedef struct { \
    count n; \
    real weight[5], scale[5], norm[5]; \
    real gen[MAXNDIM]; \
  } Set

#ifdef DEBUG
#include "common_debug.h"
#endif
#endif
