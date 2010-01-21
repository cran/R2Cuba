#ifndef __divonne_util_h__
#define  __divonne_util_h__
/*
	util.c
		Utility functions
		this file is part of Divonne
		last modified 9 Apr 09 th
*/


#include "divonne_decl.h"
/* globals */
  real *lower_, *upper_,  prdbounds_;
  count ndim_, ncomp_, nregions_;

  Integrand integrand_;
  PeakFinder peakfinder_;

  int selectedcomp_;
 number  neval_, neval_opt_, neval_cut_;
 int sign_, phase_;

 Bounds border_;

 Samples samples_[3];
 Rule rule7_, rule9_, rule11_, rule13_;
 real *xgiven_, *fgiven_, *xextra_, *fextra_;
 count ldxgiven_;
 number ngiven_, nextra_;

 Totals *totals_;

 void *voidregion_;
#define region_ ((Region *)voidregion_)
 count size_;


#define IsSobol(k) NegQ(k)
#define FIRST -INT_MAX

#define MARKMASK 0xfffffff
#define Marked(x) ((x) & ~MARKMASK)
#define Unmark(x) ((x) & MARKMASK)
#define MarkLast(x) (x | Marked(INT_MAX))
#define MEM(samples) (samples)->x

#define DIVONNETYPEDEFSET \
  typedef struct { \
    count n; \
    real weight[5], scale[5], norm[5]; \
    real gen[MAXNDIM]; \
  } Set

#ifdef DEBUG
#include "common_debug.h"
#endif
#endif
