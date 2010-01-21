
#include "common_stddecl.h"
#include "divonne_util.h"
/*********************************************************************/
 void divonneDoSample(number n, ccount ldx, ctreal *x,real *f)
{
  neval_ += n;

  while( n-- ) {
    integrand_(&ndim_, x, &ncomp_,   lower_, upper_, prdbounds_,
            f, &phase_);
    x += ldx;
    f += ncomp_;
  }
}
/*********************************************************************/
 count SampleExtra( cBounds *b)
{
  number n = nextra_;
  peakfinder_(&ndim_, b, &n, xextra_);
  divonneDoSample(n, ldxgiven_, xextra_, fextra_);
  return n;
}

