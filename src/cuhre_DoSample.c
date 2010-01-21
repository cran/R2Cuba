#include "common_stddecl.h"
#include "cuhre_util.h"
/*********************************************************************/
 void cuhreDoSample(count n, ctreal *x,  real *f)
{
  neval_ += n;
  while( n-- ) {
    integrand_(&ndim_, x, &ncomp_,  lower_, upper_, prdbounds_, f);
    x += ndim_;
    f += ncomp_;
  }
}
