/*
	decl.h
		Type declarations
		this file is part of Suave
		last modified 30 Aug 07 th
*/

//Compilation note for R interface: add ifndef
#ifndef SUAVE_DECLH
#define SUAVE_DECLH

#include "common_stddecl.h"

#define MINSAMPLES 10

#define NBINS 64

typedef unsigned char bin_t;
/* Note: bin_t must be wide enough to hold the numbers 0..NBINS */

typedef const bin_t cbin_t;

typedef real Grid[NBINS];

typedef const Grid cGrid;

typedef struct {
  real avg, err, sigsq, chisq;
} Result;

typedef const Result cResult;


typedef struct {
  real lower, upper, mid;
  Grid grid;
} Bounds;

typedef const Bounds cBounds;


#define TYPEDEFREGION \
  typedef struct region { \
    struct region *next; \
    count div, df; \
    number n; \
    Result result[NCOMP]; \
    Bounds bounds[NDIM]; \
    real fluct[NCOMP][NDIM][2]; \
    real w[]; \
  } Region


typedef void (*Integrand)(ccount *, creal *, ccount *,   creal *lower, creal *upper, creal prdbounds,  real *,creal *);

#endif
