/*
	Explore.c
		sample region, determine min and max, split if necessary
		this file is part of Divonne
		last modified 25 May 09 th
*/

#include "divonne_decl.h"
#include "divonne_util.h"
extern void Split(count iregion, int depth);
extern real FindMinimum(cBounds *b, real *xmin, real fminiii);
extern count SampleExtra( cBounds *b);
/*********************************************************************/

typedef struct {
  real fmin, fmax;
  real *xmin, *xmax;
} Extrema;

/*********************************************************************/

 bool Explore(count iregion, cSamples *samples, cint depth, cint flags)
{
#define SPLICE (flags & 1)
#define HAVESAMPLES (flags & 2)

  DIVONNETYPEDEFREGION;

  count n, dim, comp, maxcomp;
  Extrema extrema[NCOMP];
  Result *r;
  real *x, *f;
  real halfvol, maxerr;
  Region *region;
  Bounds *bounds;
  Result *result;

  /* needed as of gcc 3.3 to make gcc correctly address region #@$&! */
  /*  sizeof(*region); */


  if( SPLICE ) {
    if( nregions_ == size_ ) {
      size_ += CHUNKSIZE;
      ReAlloc(voidregion_, size_*sizeof(Region));
    }
    VecCopy(region_[nregions_].bounds, region_[iregion].bounds);
    iregion = nregions_++;
  }
  region = &region_[iregion];
  bounds = region->bounds;
  result = region->result;

  for( comp = 0; comp < ncomp_; ++comp ) {
    Extrema *e = &extrema[comp];
    e->fmin = INFTY;
    e->fmax = -INFTY;
    e->xmin = e->xmax = NULL;
  }

  if( !HAVESAMPLES ) {
    real vol = 1;
    for( dim = 0; dim < ndim_; ++dim ) {
      cBounds *b = &bounds[dim];
      vol *= b->upper - b->lower;
    }
    region->vol = vol;

    for( comp = 0; comp < ncomp_; ++comp ) {
      Result *ri = &result[comp];
      ri->fmin = INFTY;
      ri->fmax = -INFTY;
    }

    x = xgiven_;
    f = fgiven_;
    n = ngiven_;
    if( nextra_ ) n += SampleExtra(bounds);


    for( ; n; --n ) {
      for( dim = 0; dim < ndim_; ++dim ) {
        cBounds *b = &bounds[dim];
        if( x[dim] < b->lower || x[dim] > b->upper ) goto skip;
      }

      for( comp = 0; comp < ncomp_; ++comp ) {
        Extrema *e = &extrema[comp];
        ctreal y = f[comp];
        if( y < e->fmin ) e->fmin = y, e->xmin = x;
        if( y > e->fmax ) e->fmax = y, e->xmax = x;
      }
skip:
      x += ldxgiven_;
      f += ncomp_;
    }

    samples->sampler(samples, bounds, vol);
  }

  x = samples->x;
  f = samples->f;
  for( n = samples->n; n; --n ) {
    for( comp = 0; comp < ncomp_; ++comp ) {
      Extrema *e = &extrema[comp];
      ctreal y = *f++;
      if( y < e->fmin ) e->fmin = y, e->xmin = x;
      if( y > e->fmax ) e->fmax = y, e->xmax = x;
    }
    x += ndim_;
  }
  neval_opt_ -= neval_;

  halfvol = .5*region->vol;
  maxerr = -INFTY;
  maxcomp = -1;

  for( comp = 0; comp < ncomp_; ++comp ) {
    Extrema *e = &extrema[comp];
    Result *ri = &result[comp];
    real xtmp[NDIM], ftmp, err;

    if( e->xmin ) {	/* not all NaNs */
      selectedcomp_ = comp;

      sign_ = 1;
      VecCopy(xtmp, e->xmin);
      ftmp = FindMinimum(bounds, xtmp, e->fmin);
      if( ftmp < ri->fmin ) {
        ri->fmin = ftmp;
        VecCopy(ri->xmin, xtmp);
      }

      sign_ = -1;
      VecCopy(xtmp, e->xmax);
      ftmp = -FindMinimum(bounds, xtmp, -e->fmax);
      if( ftmp > ri->fmax ) {
        ri->fmax = ftmp;
        VecCopy(ri->xmax, xtmp);
      }
    }

    ri->avg = samples->avg[comp];
    ri->err = samples->err[comp];
    ri->spread = halfvol*(ri->fmax - ri->fmin);

    err = ri->spread/Max(fabs(ri->avg), NOTZERO);
    if( err > maxerr ) {
      maxerr = err;
      maxcomp = comp;
    }
  }

  neval_opt_ += neval_;

  if( maxcomp == -1 ) {		/* all NaNs */
    region->depth = 0;
    return false;
  }

  region->cutcomp = maxcomp;
  r = &region->result[maxcomp];
  if( halfvol*(r->fmin + r->fmax) > r->avg ) {
    region->fminor = r->fmin;
    region->fmajor = r->fmax;
    region->xmajor = r->xmax - (real *)region->result;
  }
  else {
    region->fminor = r->fmax;
    region->fmajor = r->fmin;
    region->xmajor = r->xmin - (real *)region->result;
  }

  region->depth = IDim(depth);

  if( !HAVESAMPLES ) {
    if( samples->weight*r->spread < r->err ||
        r->spread < totals_[maxcomp].secondspread ) region->depth = 0;
    if( region->depth == 0 )
      for( comp = 0; comp < ncomp_; ++comp )
        totals_[comp].secondspread =
          Max(totals_[comp].secondspread, result[comp].spread);
  }

  if( region->depth ) Split(iregion, region->depth);
  return true;
}


