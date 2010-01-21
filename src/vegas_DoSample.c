#include "common_stddecl.h"
#include "vegas_util.h"

/*********************************************************************/


 void vegasDoSample(number n, ctreal *w, ctreal *x,
 ctreal *lower, ctreal *upper, ctreal prdbounds, real *f)
{
  neval_ += n;
  while( n-- ) {
    integrand_(&ndim_, x, &ncomp_,  lower, upper, prdbounds, f, w++);
    x += ndim_;
    f += ncomp_;
  }


}

