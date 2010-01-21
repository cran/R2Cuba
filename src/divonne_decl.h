
/*
	decl.h
		Type declarations
		this file is part of Divonne
		last modified 25 May 09 th
*/
//Compilation note for R interface: add ifndef
#ifndef __divonne_decl_h__
#define __divonne_decl_h__


#include "common_stddecl.h"

#define EXTRAPOLATE_EPS (.25*border_.lower)
/*#define EXTRAPOLATE_EPS 0x1p-26*/


typedef struct {
  real lower, upper;
} Bounds;

typedef const Bounds cBounds;


typedef struct {
  real avg, spreadsq;
  real spread, secondspread;
  real nneed, maxerrsq, mindevsq;
  int iregion;
} Totals;


typedef struct {
  void *first, *last;
  real errcoeff[3];
  count n;
} Rule;

typedef const Rule cRule;


typedef struct samples {
  real weight;
  real *x, *f, *avg, *err;
  void (*sampler)(const struct samples *, cBounds *, ctreal);
  cRule *rule;
  count coeff;
  number n, neff;
} Samples;

typedef const Samples cSamples;


#define DIVONNETYPEDEFREGION \
  typedef struct { \
    real avg, err, spread, chisq; \
    real fmin, fmax; \
    real xmin[MAXNDIM], xmax[MAXNDIM]; \
  } Result; \
  typedef const Result cResult; \
  typedef struct region { \
    count cutcomp, depth, xmajor; \
    real fmajor, fminor, vol; \
    Bounds bounds[MAXNDIM]; \
    Result result[MAXNCOMP]; \
  } Region

#define CHUNKSIZE 4096


typedef void (*Integrand)(ccount *, ctreal *, ccount *,  ctreal *lower, ctreal *upper, ctreal prdbounds, real *, cint *);

typedef void (*PeakFinder)(ccount *, cBounds *, number *, real *);

#endif
