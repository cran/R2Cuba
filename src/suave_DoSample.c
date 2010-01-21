#include "common_stddecl.h"
#include "suave_util.h"
/*********************************************************************/

 void suaveDoSample(number n, ctreal *w, ctreal *x, real *f)
{
  neval_ += n;
  while( n-- ) {
    integrand_(&ndim_, x, &ncomp_, lower_, upper_, prdbounds_, f, w++);
    x += ndim_;
    f += ncomp_;
  }
}
