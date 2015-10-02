/*
	Integrate.c
		partition the integration region until each region
		has approximately equal spread = 1/2 vol (max - min),
		then do a main integration over all regions
		this file is part of Divonne
		last modified 8 May 09 th
*/
#include "divonne_util.h"
#include "divonne_KorobovCoeff.h"
#include "common_ChiSquare.h"



extern bool Explore(count iregion, cSamples *samples, cint depth, cint flags);
extern void divonneRuleFree(Rule *rule);
extern void IniRandom(cnumber n, cint flags, count ndim);
extern void RuleIni(Rule *rule);
extern void SamplesAlloc(Samples *samples);
extern void SamplesFree(cSamples *samples);
extern void SamplesIni(Samples *samples);
extern count SamplesLookup(Samples *samples, cint key,
			   cnumber nwant, cnumber nmax, number nmin);
extern void SampleRule(cSamples *samples, cBounds *b, ctreal vol);

extern void Split(count iregion, int depth);




extern void decodflags(cint flags, 
		int *smooth,
		int *pseudorandom,
		int *final,
		int *verbose);

#define EPS 0x1p-52 // pow(2, -52)
#define INIDEPTH 3
#define DEPTH 5
#define POSTDEPTH 15
// To rescale into the user unit
#define RESCALE(a, d) (a * (upper_[d] - lower_[d]) + lower_[d])

/*********************************************************************/

 int divonneIntegrate(ctreal epsrel, ctreal epsabs,
  cint flags, cnumber mineval, cnumber maxeval,
  int key1, int key2, int key3, ccount maxpass, 
  ctreal maxchisq, ctreal mindeviation,
  real *integral, real *erreur, real *prob)
{
  DIVONNETYPEDEFREGION;

  Totals totals[NCOMP];
  real nneed, weight;
  count dim, comp, iter, pass = 0, err, iregion;
  number nwant, nmin = INT_MAX;
  int fail = -1;

  if( VERBOSE > 1 ) {
    char s[512];
    //decode  flags: smooth ignored here
    int smooth, pseudorandom, final, verbose;
    decodflags( flags, 
		&smooth,
		&pseudorandom,
		&final,
		&verbose);
    sprintf(s, "Divonne input parameters:\n"
      "  ndim " COUNT "\n  ncomp " COUNT "\n"
      "  rel.tol " REEL "\n  abs.tol " REEL "\n"
      "  pseudo.random  %d\n  final %d\n  verbose %d\n  min.eval " NUMBER "\n  max.eval " NUMBER "\n"
      "  key1 %d\n  key2 %d\n  key3 %d\n  max.pass " COUNT "\n"
      "  border " REEL "\n  max.chisq " REEL "\n  min.deviation " REEL "\n"
      "  ngiven " NUMBER "\n  nextra " NUMBER "\n",
      ndim_, ncomp_,
      epsrel, epsabs,
      pseudorandom, final, verbose, mineval, maxeval,
      key1, key2, key3, maxpass,
      border_.lower, maxchisq, mindeviation,
      ngiven_, nextra_);
    Print(s);
  }

  size_ = CHUNKSIZE;
 MemAlloc(voidregion_, size_*sizeof(Region));
  for( dim = 0; dim < ndim_; ++dim ) {
    Bounds *b = &region_->bounds[dim];
    b->lower = 0;
    b->upper = 1;
  }
  nregions_ = 0;

  RuleIni(&rule7_);
  RuleIni(&rule9_);
  RuleIni(&rule11_);
  RuleIni(&rule13_);
  SamplesIni(&samples_[0]);
  SamplesIni(&samples_[1]);
  SamplesIni(&samples_[2]);

#ifdef MLVERSION
  if( setjmp(abort_) ) goto abort;
#endif

  /* Step 1: partition the integration region */

  if( VERBOSE ) Print("Partitioning phase:\n");

  if( IsSobol(key1) || IsSobol(key2) || IsSobol(key3) )
    IniRandom(2*maxeval, flags, ndim_);

  SamplesLookup(&samples_[0], key1,
    (number)47, (number)INT_MAX, (number)0);
  SamplesAlloc(&samples_[0]);

  totals_ = totals;
  Zap(totals);
  phase_ = 1;

  Explore(0, &samples_[0], INIDEPTH, 1);

  for( iter = 1; ; ++iter ) {
    Totals *maxtot;
    count valid;

    for( comp = 0; comp < ncomp_; ++comp ) {
      Totals *tot = &totals[comp];
      tot->avg = tot->spreadsq = 0;
      tot->spread = tot->secondspread = -INFTY;
    }

    for( iregion = 0; iregion < nregions_; ++iregion ) {
      Region *region = &region_[iregion];
      for( comp = 0; comp < ncomp_; ++comp ) {
        cResult *r = &region->result[comp];
        Totals *tot = &totals[comp];
        tot->avg += r->avg;
        tot->spreadsq += Sq(r->spread);
        if( r->spread > tot->spread ) {
          tot->secondspread = tot->spread;
          tot->spread = r->spread;
          tot->iregion = iregion;
        }
        else if( r->spread > tot->secondspread )
          tot->secondspread = r->spread;
      }
    }

    maxtot = totals;
    valid = 0;
    for( comp = 0; comp < ncomp_; ++comp ) {
      Totals *tot = &totals[comp];
      integral[comp] = tot->avg;
      valid += tot->avg == tot->avg;
      if( tot->spreadsq > maxtot->spreadsq ) maxtot = tot;
      tot->spread = sqrt(tot->spreadsq);
      erreur[comp] = tot->spread*samples_[0].weight;
    }


    if( VERBOSE ) {
      char s[128 + 64*NCOMP], *p = s;

      p += sprintf(p, 
        "Iteration " COUNT " (pass " COUNT "):  " COUNT " regions\n"
        NUMBER7 " integrand evaluations so far,\n"
        NUMBER7 " in optimizing regions,\n"
        NUMBER7 " in finding cuts",
        iter, pass, nregions_, neval_, neval_opt_, neval_cut_);

      for( comp = 0; comp < ncomp_; ++comp )
        p += sprintf(p, "\n[" COUNT "] "
          REEL " +- " REEL "\n",
          comp + 1, integral[comp], erreur[comp]);

      Print(s);
    }

    if( valid == 0 ) goto abort;	/* all NaNs */

    if( neval_ > maxeval ) break;

    nneed = maxtot->spread/MaxErr(maxtot->avg);
    if( nneed < MAXPRIME ) {
      cnumber n = neval_ + nregions_*(number)ceil(nneed);
      if( n < nmin ) {
        nmin = n;
        pass = 0;
      }
      else if( ++pass > maxpass && n >= mineval ) break;
    }

    Split(maxtot->iregion, DEPTH);
  }

  /* Step 2: do a "full" integration on each region */

/* nneed = samples_[0].neff + 1; */
  nneed = 2*samples_[0].neff;
  for( comp = 0; comp < ncomp_; ++comp ) {
    Totals *tot = &totals[comp];
    ctreal maxerr = MaxErr(tot->avg);
    tot->nneed = tot->spread/maxerr;
    nneed = Max(nneed, tot->nneed);
    tot->maxerrsq = Sq(maxerr);
    tot->mindevsq = tot->maxerrsq*Sq(mindeviation);
  }
  nwant = (number)Min(ceil(nneed), MARKMASK/40.);

  err = SamplesLookup(&samples_[1], key2, nwant,
    (maxeval - neval_)/nregions_ + 1, samples_[0].n + 1);

  /* the number of points needed to reach the desired accuracy */
  fail = Unmark(err)*nregions_;

  if( Marked(err) ) {
    if( VERBOSE ) Print("\nNot enough samples left for main integration.");
    for( comp = 0; comp < ncomp_; ++comp )
      prob[comp] = -999;
    weight = samples_[0].weight;
  }
  else {
    bool can_adjust = (key3 == 1 && samples_[1].sampler != SampleRule &&
      (key2 < 0 || samples_[1].neff < MAXPRIME));
    count df, nlimit;

    SamplesAlloc(&samples_[1]);

    if( VERBOSE ) {
      char s[128];
      sprintf(s, "\nMain integration on " COUNT
        " regions with " NUMBER " samples per region.",
        nregions_, samples_[1].neff);
      Print(s);
    }

    ResClear(integral);
    ResClear(erreur);
    ResClear(prob);

    nlimit = maxeval - nregions_*samples_[1].n;
    df = 0;

    for( iregion = 0; iregion < nregions_; ++iregion ) {
      Region *region = &region_[iregion];
      char s[64*NDIM + 256*NCOMP], *p = s;
      int todo;

refine:
      phase_ = 2;
      samples_[1].sampler(&samples_[1], region->bounds, region->vol);

      if( can_adjust )
        for( comp = 0; comp < ncomp_; ++comp )
          totals[comp].spreadsq -= Sq(region->result[comp].spread);

      nlimit += samples_[1].n;
      todo = 0;

      for( comp = 0; comp < ncomp_; ++comp ) {
        cResult *r = &region->result[comp];
        Totals *tot = &totals[comp];

        samples_[0].avg[comp] = r->avg;
        samples_[0].err[comp] = r->err;

        if( neval_ < nlimit ) {
          ctreal avg2 = samples_[1].avg[comp];
          ctreal err2 = samples_[1].err[comp];
          ctreal diffsq = Sq(avg2 - r->avg);

#define Var(s) Sq((s.err[comp] == 0) ? r->spread*s.weight : s.err[comp])

          if( err2*tot->nneed > r->spread ||
              diffsq > Max(maxchisq*(Var(samples_[0]) + Var(samples_[1])),
                           EPS*Sq(avg2)) ) {
            if( key3 && diffsq > tot->mindevsq ) {
              if( key3 == 1 ) {
                ccount xregion = nregions_;

                if( VERBOSE > 2 ) Print("\nSplit");

                phase_ = 1;
                Explore(iregion, &samples_[1], POSTDEPTH, 2);

                if( can_adjust ) {
                  number nnew;
                  count ireg, xreg;

                  for( ireg = iregion, xreg = xregion;
                       ireg < nregions_; ireg = xreg++ ) {
                    cResult *result = region_[ireg].result;
                    count c;
                    for( c = 0; c < ncomp_; ++c )
                      totals[c].spreadsq += Sq(result[c].spread);
                  }

                  nnew = (tot->spreadsq/Sq(MARKMASK) > tot->maxerrsq) ?
                    MARKMASK :
                    (number)ceil(sqrt(tot->spreadsq/tot->maxerrsq));
                  if( nnew > nwant + nwant/64 ) {
                    ccount erri = SamplesLookup(&samples_[1], key2, nnew,
                      (maxeval - neval_)/nregions_ + 1, samples_[1].n);
                    fail += Unmark(erri)*nregions_;
                    nwant = nnew;
                    SamplesFree(&samples_[1]);
                    SamplesAlloc(&samples_[1]);

                    if( key2 > 0 && samples_[1].neff >= MAXPRIME )
                      can_adjust = false;

                    if( VERBOSE > 2 ) {
                      char si[128];
                      sprintf(si, "Sampling remaining " COUNT
                        " regions with " NUMBER " points per region.",
                        nregions_, samples_[1].neff);
                      Print(si);
                    }
                  }
                }

                goto refine;
              }
              todo |= 3;
            }
            todo |= 1;
          }
        }
      }

      if( can_adjust ) {
        for( comp = 0; comp < ncomp_; ++comp )
          totals[comp].maxerrsq -=
            Sq(region->result[comp].spread*samples_[1].weight);
      }

      switch( todo ) {
      case 1:	/* get spread right */
        Explore(iregion, &samples_[1], 0, 2);
        break;

      case 3:	/* sample region again with more points */
        if( MEM(&samples_[2]) == NULL ) {
          SamplesLookup(&samples_[2], key3,
            nwant, (number)INT_MAX, (number)0);
          SamplesAlloc(&samples_[2]);
        }
        phase_ = 3;
        samples_[2].sampler(&samples_[2], region->bounds, region->vol);
        Explore(iregion, &samples_[2], 0, 2);
        ++region->depth;	/* misused for df here */
        ++df;
      }

      ++region->depth;	/* misused for df here */

      if( VERBOSE > 2 ) {
        for( dim = 0; dim < ndim_; ++dim ) {
          cBounds *b = &region->bounds[dim];
          p += sprintf(p,
            (dim == 0) ? "\nRegion (" REALF ") - (" REALF ")" :
                         "\n       (" REALF ") - (" REALF ")",
            RESCALE(b->lower,dim),  RESCALE(b->upper,dim));
        }
      }

      for( comp = 0; comp < ncomp_; ++comp ) {
        Result *r = &region->result[comp];

        ctreal x1 = samples_[0].avg[comp];
        ctreal s1 = Var(samples_[0]);
        ctreal x2 = samples_[1].avg[comp];
        ctreal s2 = Var(samples_[1]);
        ctreal r2 = (s1 == 0) ? Sq(samples_[1].neff*samples_[0].weight) : s2/s1;

        real norm = 1 + r2;
        real avg = x2 + r2*x1;
        real sigsq = s2;
        real chisq = Sq(x2 - x1);
        real chiden = s1 + s2;

        if( todo == 3 ) {
          ctreal x3 = samples_[2].avg[comp];
          ctreal s3 = Var(samples_[2]);
          ctreal r3 = (s2 == 0) ? Sq(samples_[2].neff*samples_[1].weight) : s3/s2;

          norm = 1 + r3*norm;
          avg = x3 + r3*avg;
          sigsq = s3;
          chisq = s1*Sq(x3 - x2) + s2*Sq(x3 - x1) + s3*chisq;
          chiden = s1*s2 + s3*chiden;
        }

        avg = LAST ? r->avg : (sigsq *= norm = 1/norm, avg*norm);
        if( chisq > EPS ) chisq /= Max(chiden, NOTZERO);

#define Out(s) s.avg[comp], r->spread*s.weight, s.err[comp]

        if( VERBOSE > 2 ) {
          p += sprintf(p, "\n[" COUNT "] "
            REEL " +- " REEL "(" REEL ")\n    "
            REEL " +- " REEL "(" REEL ")",
            comp + 1, Out(samples_[0]), Out(samples_[1]));
          if( todo == 3 ) p += sprintf(p, "\n    "
            REEL " +- " REEL "(" REEL ")\n",
            Out(samples_[2]));
          p += sprintf(p, "  \tchisq " REEL, chisq);
        }

        integral[comp] += avg;
        erreur[comp] += sigsq;
        prob[comp] += chisq;

        r->avg = avg;
        r->spread = sqrt(sigsq);
        r->chisq = chisq;
      }

      if( VERBOSE > 2 ) Print(s);
    }

    for( comp = 0; comp < ncomp_; ++comp )
      erreur[comp] = sqrt(erreur[comp]);

    df += nregions_;

    if( VERBOSE > 2 ) {
      char s[16 + 128*NCOMP], *p = s;

      p += sprintf(p, "\nTotals:");

      for( comp = 0; comp < ncomp_; ++comp )
        p += sprintf(p, "\n[" COUNT "] "
          REEL " +- " REEL "  \tchisq " REEL " (" COUNT " df)\n",
          comp + 1, integral[comp], erreur[comp], prob[comp], df);

      Print(s);
    }

    for( comp = 0; comp < ncomp_; ++comp )
      prob[comp] = ChiSquare(prob[comp], df);

    weight = 1;
  }

#ifdef MLVERSION
  if( REGIONS ) {
    MLPutFunction(stdlink, "List", 2);
    MLPutFunction(stdlink, "List", nregions_);
    for( iregion = 0; iregion < nregions_; ++iregion ) {
      Region *region = &region_[iregion];
      cBounds *b = region->bounds;
      real lower[NDIM], upper[NDIM];

      for( dim = 0; dim < ndim_; ++dim ) {
        lower[dim] = b[dim].lower;
        upper[dim] = b[dim].upper;
      }

      MLPutFunction(stdlink, "Cuba`Divonne`region", 4);

      MLPutRealList(stdlink, lower, ndim_);
      MLPutRealList(stdlink, upper, ndim_);

      MLPutFunction(stdlink, "List", ncomp_);
      for( comp = 0; comp < ncomp_; ++comp ) {
        cResult *r = &region->result[comp];
        real res[] = {r->avg, r->spread*weight, r->chisq};
        MLPutRealList(stdlink, res, Elements(res));
      }

      MLPutInteger(stdlink, region->depth);  /* misused for df */
    }
  }
#endif

abort:

  SamplesFree(&samples_[2]);
  SamplesFree(&samples_[1]);
  SamplesFree(&samples_[0]);
  divonneRuleFree(&rule13_);
  divonneRuleFree(&rule11_);
  divonneRuleFree(&rule9_);
  divonneRuleFree(&rule7_);

  free(region_);

  return fail;
}

