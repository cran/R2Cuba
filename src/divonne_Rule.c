/*
	Rule.c
		integration with cubature rules
		code lifted with minor modifications from DCUHRE
		by J. Berntsen, T. Espelid, and A. Genz
		this file is part of Divonne
		last modified 9 Feb 05 th
*/

#include "divonne_decl.h"
#include "divonne_util.h"

extern void divonneDoSample(number n, ccount ldx, ctreal *x,real *f);

enum { nrules = 5 };


/*********************************************************************/

 void RuleIni(Rule *rule)
{
  rule->first = NULL;
}

/*********************************************************************/

 void divonneRuleFree(Rule *rule)
{
  if( rule->first ) free(rule->first);
}


/*********************************************************************/

 real *divonneExpandFS(cBounds *b, real *g, real *x)
{
  count dim, ndim = ndim_;

next:
  /* Compute centrally symmetric sum for permutation of G */

  for( dim = 0; dim < ndim; ++dim )
    *x++ = (.5 + g[dim])*b[dim].lower + (.5 - g[dim])*b[dim].upper;

  for( dim = 0; dim < ndim; ) {
    g[dim] = -g[dim];
    if( g[dim++] < 0 ) goto next;
  }

  /* Find next distinct permutation of G and loop back for next sum.
     Permutations are generated in reverse lexicographic order. */

  for( dim = 1; dim < ndim; ++dim ) {
    ctreal gd = g[dim];
    count big = dim - 1;
    if( g[big] > gd ) {
      count i, j = dim, lastbig = big;
      for( i = 0; i < --j; ++i ) {
        ctreal tmp = g[i];
        g[i] = g[j];
        g[j] = tmp;
        if( tmp <= gd ) --big;
        if( g[i] > gd ) lastbig = i;
      }
      if( g[big] <= gd ) big = lastbig;
      g[dim] = g[big];
      g[big] = gd;
      goto next;
    }
  }

  /* Restore original order to generators */

  for( dim = 0; dim < --ndim; ++dim ) {
    ctreal tmp = g[dim];
    g[dim] = g[ndim];
    g[ndim] = tmp;
  }

  return x;
}



/*********************************************************************/

 void SampleRule(cSamples *samples, cBounds *b, ctreal vol)
{
  DIVONNETYPEDEFSET;

  real *x = samples->x, *f = samples->f;
  Set *first = (Set *)samples->rule->first;
  Set *last = (Set *)samples->rule->last;
  Set *s;
  ctreal *errcoeff = samples->rule->errcoeff;
  count comp, rul, n;

  for( s = first; s <= last; ++s )
    if( s->n ) x = divonneExpandFS(b, s->gen, x);

  divonneDoSample(samples->n, ndim_, samples->x, f);

  for( comp = 0; comp < ncomp_; ++comp ) {
    real sum[nrules];
    ctreal *f1 = f++;

    Zap(sum);
    for( s = first; s <= last; ++s )
      for( n = s->n; n; --n ) {
        ctreal fun = *f1;
        f1 += ncomp_;
        for( rul = 0; rul < nrules; ++rul )
          sum[rul] += fun*s->weight[rul];
      }

    /* Search for the null rule, in the linear space spanned by two
       successive null rules in our sequence, which gives the greatest
       error estimate among all normalized (1-norm) null rules in this
       space. */

    for( rul = 1; rul < nrules - 1; ++rul ) {
      real maxerr = 0;
      for( s = first; s <= last; ++s )
        maxerr = Max(maxerr,
          fabs(sum[rul + 1] + s->scale[rul]*sum[rul])*s->norm[rul]);
      sum[rul] = maxerr;
    }

    samples->avg[comp] = vol*sum[0];
    samples->err[comp] = vol*(
      (errcoeff[0]*sum[1] <= sum[2] && errcoeff[0]*sum[2] <= sum[3]) ?
        errcoeff[1]*sum[1] :
        errcoeff[2]*Max(Max(sum[1], sum[2]), sum[3]));
  }
}
