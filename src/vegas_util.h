//Compilation note for R interface: move into a .h
//Compilation note for R interface: add ifndef
#ifndef __vegas_util_h__
#define __vegas_util_h__
/*
	util.c
		Utility functions
		this file is part of Vegas
		last modified 2 Mar 06 th
*/


#include "vegas_decl.h"

static count ndim_, ncomp_;
static number neval_;
static Grid *gridptr_[MAXGRIDS];
static count griddim_[MAXGRIDS];
/* Compilation note for R interface: remove
 int EXPORT(vegasnbatch) = 1000;
 int EXPORT(vegasgridno) = 0;*/
int EXPORT(vegasnbatch);
int EXPORT(vegasgridno);
char EXPORT(vegasstate)[MAXSTATESIZE] = "";


#define SamplesAlloc(p, n) \
  MemAlloc(p, (n)*((ndim_ + ncomp_ + 1)*sizeof(real) + ndim_*sizeof(bin_t)))

#ifdef DEBUG
#include "common_debug.h" 
#endif
#endif
