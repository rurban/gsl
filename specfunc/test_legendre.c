/* specfunc/test_legendre.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2004 Gerard Jungman
 * Copyright (C) 2013, 2021 Patrick Alken
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* Author:  G. Jungman */

#include <config.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_sf.h>
#include "test_sf.h"

static double
test_legendre_dx(const size_t l)
{
  const double dx_max = 0.4;
  double dx;

  if (l < 1000)
    dx = exp((double)l / 1000.0) / exp(2.0);
  else
    dx = dx_max;

  return dx;
}

/*
test_legendre_sum()
  This routine computes the sum:

  Sum_{m=0}^l [P(l,m)(x)]^2

This sum should equate to 1.0 for Schmidt semi-normalized
ALFs for all l.
*/

static double
test_legendre_sum(const int indexl, const size_t lmax,
                  const size_t l, double *p)
{
  double sum = 0.0;
  size_t idx;
  size_t m;

  for (m = 0; m <= l; ++m)
    {
      idx = indexl ? gsl_sf_legendre_array_index(l, m) : gsl_sf_legendre_array_index_m(l, m, lmax);
      sum += p[idx] * p[idx];
    }

  return sum;
}

static void
test_value(const size_t flags, const size_t lmax, const size_t l, const size_t m,
           const double *p, const double expected, const double tol,
           const char *desc, const char *desc2)
{
  size_t idx;
  double value;

  if (l > lmax)
    return;

  idx = (flags & GSL_SF_LEGENDRE_FLG_INDEXL) ? gsl_sf_legendre_array_index(l, m) : gsl_sf_legendre_array_index_m(l, m, lmax);
  value = p[idx];

  gsl_test_rel(value, expected, tol, "%s %s lmax=%zu l=%zu m=%zu", desc, desc2, lmax, l, m);
}

/*
test_legendre_cross()
  This routine cross tests the ALF normalizations against each other
*/

static int
test_legendre_cross(const double tol, const size_t lmax,
                    const size_t flags, const char *desc)
{
  int s = 0;
  const size_t lmax_P = GSL_MIN(lmax, 140);
  const size_t nlm = gsl_sf_legendre_nlm(lmax);
  const size_t plm_size = gsl_sf_legendre_array_n(lmax);
  const double dx = test_legendre_dx(lmax);
  const double b = 1.0 / M_SQRT2 / M_SQRTPI; /* Ylm = b * Nlm */
  const int indexl = (flags & GSL_SF_LEGENDRE_FLG_INDEXL);
  size_t l, m;
  double x;
  double *p = malloc(plm_size * sizeof(double));
  double *dp = malloc(nlm * sizeof(double));
  double *d2p = malloc(nlm * sizeof(double));
  double *p_schmidt = malloc(plm_size * sizeof(double));
  double *dp_schmidt = malloc(nlm * sizeof(double));
  double *d2p_schmidt = malloc(nlm * sizeof(double));
  double *p_spharm = malloc(plm_size * sizeof(double));
  double *dp_spharm = malloc(nlm * sizeof(double));
  double *d2p_spharm = malloc(nlm * sizeof(double));
  double *p_full = malloc(plm_size * sizeof(double));
  double *dp_full = malloc(nlm * sizeof(double));
  double *d2p_full = malloc(nlm * sizeof(double));

  gsl_sf_legendre_precompute(GSL_SF_LEGENDRE_NONE, lmax, flags, p);
  gsl_sf_legendre_precompute(GSL_SF_LEGENDRE_SCHMIDT, lmax, flags, p_schmidt);
  gsl_sf_legendre_precompute(GSL_SF_LEGENDRE_SPHARM, lmax, flags, p_spharm);
  gsl_sf_legendre_precompute(GSL_SF_LEGENDRE_FULL, lmax, flags, p_full);

  for (x = -1.0; x <= 1.0; x += dx)
    {
      gsl_sf_legendre_deriv2_alt_arrayx(GSL_SF_LEGENDRE_NONE, lmax, x, p, dp, d2p);
      gsl_sf_legendre_deriv2_alt_arrayx(GSL_SF_LEGENDRE_SCHMIDT, lmax, x,
                                        p_schmidt, dp_schmidt, d2p_schmidt);
      gsl_sf_legendre_deriv2_alt_arrayx(GSL_SF_LEGENDRE_SPHARM, lmax, x,
                                        p_spharm, dp_spharm, d2p_spharm);
      gsl_sf_legendre_deriv2_alt_arrayx(GSL_SF_LEGENDRE_FULL, lmax, x,
                                        p_full, dp_full, d2p_full);

      for (l = 0; l <= lmax; ++l)
        {
          double a_lm = sqrt(2.0 / (double)l / (l + 1.0));
          size_t l0idx = indexl ? gsl_sf_legendre_array_index(l, 0) : gsl_sf_legendre_array_index_m(l, 0, lmax);
          double cl0 = sqrt((4.0 * M_PI) / (2.0 * l + 1.0));
          double clm = sqrt((8.0 * M_PI) / (2.0 * l + 1.0));

          if (lmax <= lmax_P)
            {
              /* test S(l,0) = P(l,0) */
              gsl_test_rel(p[l0idx], p_schmidt[l0idx], tol,
                           "%s l=%zu, m=0, x=%f", desc, l, x);
              gsl_test_rel(dp[l0idx], dp_schmidt[l0idx], tol,
                           "%s deriv l=%zu, m=0, x=%f", desc, l, x);

              if (l > 1)
                {
                  gsl_test_rel(d2p[l0idx], d2p_schmidt[l0idx], tol,
                               "%s deriv2 l=%zu, m=0, x=%f", desc, l, x);
                }
              else
                {
                  gsl_test_abs(d2p[l0idx], d2p_schmidt[l0idx], tol,
                               "%s deriv2 l=%zu, m=0, x=%f", desc, l, x);
                }
            }

          gsl_test_rel(cl0 * p_spharm[l0idx], p_schmidt[l0idx], tol,
                       "%s spharm l=%zu, m=%zu, x=%f", desc, l, m, x);
          gsl_test_rel(cl0 * dp_spharm[l0idx], dp_schmidt[l0idx], tol,
                       "%s spharm deriv l=%zu, m=%zu, x=%f", desc, l, m, x);
          gsl_test_rel(cl0 * d2p_spharm[l0idx], d2p_schmidt[l0idx], tol,
                       "%s spharm deriv2 l=%zu, m=%zu, x=%f", desc, l, m, x);

          gsl_test_rel(b * p_full[l0idx], p_spharm[l0idx], tol,
                       "%s full l=%zu, m=%zu, x=%f", desc, l, m, x);
          gsl_test_rel(b * dp_full[l0idx], dp_spharm[l0idx], tol,
                       "%s full deriv l=%zu, m=%zu, x=%f", desc, l, m, x);
          gsl_test_rel(b * d2p_full[l0idx], d2p_spharm[l0idx], tol,
                       "%s full deriv2 l=%zu, m=%zu, x=%f", desc, l, m, x);

          /* test S(l,m) = a_{lm} * P(l,m) for m > 0 */
          for (m = 1; m <= l; ++m)
            {
              size_t lmidx = indexl ? gsl_sf_legendre_array_index(l, m) : gsl_sf_legendre_array_index_m(l, m, lmax);

              if (lmax <= lmax_P)
                {
                  gsl_test_rel(a_lm * p[lmidx], p_schmidt[lmidx], tol,
                               "%s schmidt l=%zu, m=%zu, x=%f", desc, l, m, x);
                  gsl_test_abs(a_lm * dp[lmidx], dp_schmidt[lmidx], tol,
                               "%s schmidt deriv l=%zu, m=%zu, x=%f", desc, l, m, x);
                  gsl_test_abs(a_lm * d2p[lmidx], d2p_schmidt[lmidx], tol,
                               "%s schmidt deriv2 l=%zu, m=%zu, x=%f", desc, l, m, x);

                  a_lm /= sqrt((double) (l + m + 1)) *
                          sqrt((double) (l - m));
                }

              gsl_test_rel(clm * p_spharm[lmidx], p_schmidt[lmidx], tol,
                           "%s spharm l=%zu, m=%zu, x=%f", desc, l, m, x);
              gsl_test_rel(clm * dp_spharm[lmidx], dp_schmidt[lmidx], tol,
                           "%s spharm deriv l=%zu, m=%zu, x=%f", desc, l, m, x);
              gsl_test_rel(clm * d2p_spharm[lmidx], d2p_schmidt[lmidx], tol,
                           "%s spharm deriv2 l=%zu, m=%zu, x=%f", desc, l, m, x);

              gsl_test_rel(b * p_full[lmidx], p_spharm[lmidx], tol,
                           "%s full l=%zu, m=%zu, x=%f", desc, l, m, x);
              gsl_test_rel(b * dp_full[lmidx], dp_spharm[lmidx], tol,
                           "%s full deriv l=%zu, m=%zu, x=%f", desc, l, m, x);
              gsl_test_rel(b * d2p_full[lmidx], d2p_spharm[lmidx], tol,
                           "%s full deriv2 l=%zu, m=%zu, x=%f", desc, l, m, x);
            }
        }
    }

  free(p);
  free(dp);
  free(d2p);
  free(p_schmidt);
  free(dp_schmidt);
  free(d2p_schmidt);
  free(p_spharm);
  free(dp_spharm);
  free(d2p_spharm);
  free(p_full);
  free(dp_full);
  free(d2p_full);

  return s;
}

int
test_legendre_eps(const double tol, const gsl_sf_legendre_t norm, const size_t lmax,
                  const size_t flags, const char * desc)
{
  int s = 0;
  const size_t nlm = gsl_sf_legendre_nlm(lmax);
  const size_t plm_size = gsl_sf_legendre_array_n(lmax);
  const double dx = test_legendre_dx(lmax);
  const int indexl = (flags & GSL_SF_LEGENDRE_FLG_INDEXL);
  double *Plm = malloc(plm_size * sizeof(double));
  double *dPlm = malloc(nlm * sizeof(double));
  double *d2Plm = malloc(nlm * sizeof(double));
  double *Plm_alt = malloc(plm_size * sizeof(double));
  double *dPlm_alt = malloc(nlm * sizeof(double));
  double *d2Plm_alt = malloc(nlm * sizeof(double));
  size_t l, i;
  char buf[32];
  double x;

  const double * P_norm, * dP_norm, * d2P_norm;

  const size_t nx = 4;
  const double x_vals[] = { -0.38, 0.86, -1.0, 1.0 };

  const double P_schmidt[] = {
     1.0000000000000000,  1.0000000000000000,  1.0000000000000000,  1.0000000000000000, /* P(0,0) */
    -0.3800000000000000,  0.8600000000000000, -1.0000000000000000,  1.0000000000000000, /* P(1,0) */
     0.924986486387774,   0.510294032886923,   0.0000000000000000,  0.0000000000000000, /* P(1,1) */
    -0.2834000000000000,  0.6094000000000000,  1.0000000000000000,  1.0000000000000000, /* P(2,0) */
    -0.6088069644805322,  0.7601154649130617,  0.0000000000000000,  0.0000000000000000, /* P(2,1) */
     0.7409713354779655,  0.2255130151454678,  0.0000000000000000,  0.0000000000000000, /* P(2,2) */
     0.432820000000000,   0.3001400000000000, -1.0000000000000000,  1.0000000000000000, /* P(3,0) */
    -0.1574692712880831,  0.843098019568306,   0.0000000000000000,  0.0000000000000000, /* P(3,1) */
    -0.6296076646928626,  0.433665691241537,   0.0000000000000000,  0.0000000000000000, /* P(3,2) */
     0.6256712113882177,  0.105051311462542,   0.0000000000000000,  0.0000000000000000  /* P(3,3) */
  };
  const double dP_schmidt[] = {
     0.0000000000000000,  0.0000000000000000,  0.0000000000000000,  0.0000000000000000, /* d/dt P(0,0) */
    -0.9249864863877740, -0.5102940328869230,  0.0000000000000000,  0.0000000000000000, /* d/dt P(1,0) */
    -0.3800000000000000,  0.8600000000000000, -1.0000000000000000,  1.0000000000000000, /* d/dt P(1,1) */
     1.0544845944820630, -1.3165586048482610,  0.0000000000000000,  0.0000000000000000, /* d/dt P(2,0) */
    -1.2318345343429852,  0.8299987469870067,  1.732050807568877,   1.732050807568877,  /* d/dt P(2,1) */
    -0.6088069644805323,  0.7601154649130618,  0.0000000000000000,  0.0000000000000000, /* d/dt P(2,2) */
     0.3857193648237012, -2.0651599510933770,  0.0000000000000000,  0.0000000000000000, /* d/dt P(3,0) */
     2.0556852768359271,  0.0495041877016476, -2.449489742783178,   2.449489742783178,  /* d/dt P(3,1) */
    -1.015268386703732,   1.204393961359827,   0.0000000000000000,  0.0000000000000000, /* d/dt P(3,2) */
    -0.7711087583214189,  0.5311298312465608,  0.0000000000000000,  0.0000000000000000  /* d/dt P(3,3) */
  };
  const double d2P_schmidt[] = {
     0.0000000000000000,  0.0000000000000000,  0.0000000000000000,  0.0000000000000000, /* d2/dt2 P(0,0) */
     0.3800000000000000, -0.8600000000000000,  1.0000000000000000, -1.0000000000000000, /* d2/dt2 P(1,0) */
    -0.924986486387774,  -0.510294032886923,   0.0000000000000000,  0.0000000000000000, /* d2/dt2 P(1,1) */
     2.1336000000000000, -1.4376000000000000, -3.0000000000000000, -3.0000000000000000, /* d2/dt2 P(2,0) */
     2.435227857922129,  -3.040461859652247,   0.0000000000000000,  0.0000000000000000, /* d2/dt2 P(2,1) */
    -1.231834534342985,   0.829998746987006,   1.732050807568877,   1.732050807568877,  /* d2/dt2 P(2,2) */
    -5.0353800000000000, -0.1212600000000000,  6.0000000000000000, -6.0000000000000000, /* d2/dt2 P(3,0) */
     2.550095896902701,  -6.962902176434765,   0.0000000000000000,  0.0000000000000000, /* d2/dt2 P(3,1) */
     4.194735310676946,  -0.572225543435454,  -3.872983346207417,   3.872983346207417,  /* d2/dt2 P(3,2) */
    -1.243444749701409,   1.475075327310439,   0.0000000000000000,  0.0000000000000000  /* d2/dt2 P(3,3) */
  };

  const double P_spharm[] = {
     0.2820947917738781,  0.2820947917738781,  0.2820947917738781,  0.2820947917738781, /* P(0,0) */
    -0.1856689545231096,  0.4201981602365111, -0.4886025119029199,  0.4886025119029199, /* P(1,0) */
     0.3195774193870231,  0.1763036028725652,  0.0000000000000000,  0.0000000000000000, /* P(1,1) */
    -0.1787639391851283,  0.3843992397297712,  0.6307831305050398,  0.6307831305050398, /* P(2,0) */
    -0.2715467968528703,  0.3390350830030173,  0.0000000000000000,  0.0000000000000000, /* P(2,1) */
     0.3304962072510409,  0.1005858022068386,  0.0000000000000000,  0.0000000000000000, /* P(2,2) */
     0.3230363605433075,  0.2240102889271943, -0.7463526651802306,  0.7463526651802306, /* P(3,0) */
    -0.0831045702267447,  0.4449458488130187,  0.0000000000000000,  0.0000000000000000, /* P(3,1) */
    -0.3322760939818003,  0.228867515534375,   0.0000000000000000,  0.0000000000000000, /* P(3,2) */
     0.3301986266929495,  0.05544093790133093, 0.0000000000000000,  0.0000000000000000  /* P(3,3) */
  };
  const double dP_spharm[] = {
     0.0000000000000000,  0.0000000000000000,  0.0000000000000000,  0.0000000000000000, /* d/dt P(0,0) */
    -0.4519507207253225, -0.2493309462776218,  0.0000000000000000,  0.0000000000000000, /* d/dt P(1,0) */
    -0.1312877767991075,  0.2971249685453485, -0.3454941494713355,  0.3454941494713355, /* d/dt P(1,1) */
     0.6651510935767333, -0.830462958259534,   0.0000000000000000,  0.0000000000000000, /* d/dt P(2,0) */
    -0.5494364249577846,  0.3702051952190248,  0.7725484040463791,  0.7725484040463791, /* d/dt P(2,1) */
    -0.2715467968528703,  0.3390350830030173,  0.0000000000000000,  0.0000000000000000, /* d/dt P(2,2) */
     0.2878826759477952, -1.541337633522017,   0.0000000000000000,  0.0000000000000000, /* d/dt P(3,0) */
     1.0848900236564757,  0.02612588608378763, -1.292720736456603,  1.292720736456603,  /* d/dt P(3,1) */
    -0.5358089375257002,  0.6356201544832385,  0.0000000000000000,  0.0000000000000000, /* d/dt P(3,2) */
    -0.4069534419902396,  0.2803043158788606,  0.0000000000000000,  0.0000000000000000  /* d/dt P(3,3) */
  };
  const double d2P_spharm[] = {
     0.0000000000000000,  0.0000000000000000,  0.0000000000000000,  0.0000000000000000, /* d2/dt2 P(0,0) */
     0.1856689545231096, -0.4201981602365112,  0.4886025119029199, -0.4886025119029199, /* d2/dt2 P(1,0) */
    -0.3195774193870231, -0.1763036028725652,  0.0000000000000000,  0.0000000000000000, /* d2/dt2 P(1,1) */
     1.345838887245553,  -0.906813828414046,  -1.89234939151512,   -1.89234939151512,   /* d2/dt2 P(2,0) */
     1.086187187411481,  -1.356140332012069,   0.0000000000000000,  0.0000000000000000, /* d2/dt2 P(2,1) */
    -0.5494364249577846,  0.3702051952190248,  0.7725484040463791,  0.7725484040463791, /* d2/dt2 P(2,2) */
    -3.758169283195231,  -0.0905027241797544,  4.478115991081385,  -4.478115991081385,  /* d2/dt2 P(3,0) */
     1.345815737988507,  -3.674678800315671,   0.0000000000000000,  0.0000000000000000, /* d2/dt2 P(3,1) */
     2.213775883746908,  -0.3019926203441294, -2.043970952866565,   2.043970952866565,  /* d2/dt2 P(3,2) */
    -0.6562292482803778,  0.7784725243564758,  0.0000000000000000,  0.0000000000000000  /* d2/dt2 P(3,3) */
  };

  const double P_full[] = {
     0.7071067811865475,  0.7071067811865475,  0.7071067811865475,  0.7071067811865475, /* P(0,0) */
    -0.4654030511288038,  1.053280589396766,  -1.224744871391589,   1.224744871391589,  /* P(1,0) */
     0.801061795369121,   0.441927595879687,   0.0000000000000000,  0.0000000000000000, /* P(1,1) */
    -0.4480947444458591,  0.963546003053305,   1.581138830084189,   1.581138830084189,  /* P(2,0) */
    -0.680666878876885,   0.849834925147231,   0.0000000000000000,  0.0000000000000000, /* P(2,1) */
     0.828431137753766,   0.2521312158381029,  0.0000000000000000,  0.0000000000000000, /* P(2,2) */
     0.809732075071749,   0.5615105240331647, -1.87082869338697,    1.87082869338697,   /* P(3,0) */
    -0.2083122654814157,  1.115313845314403,   0.0000000000000000,  0.0000000000000000, /* P(3,1) */
    -0.832892652158728,   0.5736857855830142,  0.0000000000000000,  0.0000000000000000, /* P(3,2) */
     0.827685213912874,   0.1389698225155376,  0.0000000000000000,  0.0000000000000000  /* P(3,3) */
  };
  const double dP_full[] = {
     0.0000000000000000,  0.0000000000000000,  0.0000000000000000,  0.0000000000000000, /* d/dt P(0,0) */
    -1.132872455309952,  -0.6249799996799897,  0.0000000000000000,  0.0000000000000000, /* d/dt P(1,0) */
    -0.3290896534380868,  0.7447818472546172, -0.866025403784439,   0.866025403784439,  /* d/dt P(1,1) */
     1.66728653806117,   -2.081661932207053,   0.0000000000000000,  0.0000000000000000, /* d/dt P(2,0) */
    -1.377232877911357,   0.927966809751297,   1.936491673103708,   1.936491673103708,  /* d/dt P(2,1) */
    -0.6806668788768851,  0.849834925147231,   0.0000000000000000,  0.0000000000000000, /* d/dt P(2,2) */
     0.7216148553071772, -3.863560492939122,   0.0000000000000000,  0.0000000000000000, /* d/dt P(3,0) */
     2.719416008162415,   0.06548788475741096, -3.24037034920393,   3.24037034920393,   /* d/dt P(3,1) */
    -1.343073832601915,   1.59326345115301,    0.0000000000000000,  0.0000000000000000, /* d/dt P(3,2) */
    -1.020081004151141,   0.7026187236830514,  0.0000000000000000,  0.0000000000000000  /* d/dt P(3,3) */
  };
  const double d2P_full[] = {
     0.0000000000000000,  0.0000000000000000,  0.0000000000000000,  0.0000000000000000, /* d2/dt2 P(0,0) */
     0.4654030511288039, -1.053280589396766,   1.224744871391589,  -1.224744871391589,  /* d2/dt2 P(1,0) */
    -0.801061795369121,  -0.441927595879687,   0.0000000000000000,  0.0000000000000000, /* d2/dt2 P(1,1) */
     3.373517807867627,  -2.273045182129031,  -4.743416490252569,  -4.743416490252569,  /* d2/dt2 P(2,0) */
     2.722667515507541,  -3.399339700588925,   0.0000000000000000,  0.0000000000000000, /* d2/dt2 P(2,1) */
    -1.377232877911357,   0.927966809751297,   1.936491673103708,   1.936491673103708,  /* d2/dt2 P(2,2) */
    -9.42033338610689,   -0.2268566873601031,  11.22497216032182,  -11.22497216032182,  /* d2/dt2 P(3,0) */
     3.373459781285376,  -9.21105378105839,    0.000000000000000,   0.000000000000000,  /* d2/dt2 P(3,1) */
     5.549113223896231,  -0.7569832408844999, -5.123475382979799,   5.123475382979799,  /* d2/dt2 P(3,2) */
    -1.644922788379442,   1.951341240575312,   0.0000000000000000,  0.0000000000000000  /* d2/dt2 P(3,3) */
  };

  const double P_none[] = {
     1.0000000000000000,  1.0000000000000000,  1.0000000000000000,  1.0000000000000000, /* P(0,0) */
    -0.3800000000000000,  0.8600000000000000, -1.0000000000000000,  1.0000000000000000, /* P(1,0) */
     0.924986486387774,   0.510294032886923,   0.0000000000000000,  0.0000000000000000, /* P(1,1) */
    -0.2834000000000000,  0.6094000000000000,  1.0000000000000000,  1.0000000000000000, /* P(2,0) */
    -1.054484594482063,   1.316558604848261,   0.0000000000000000,  0.0000000000000000, /* P(2,1) */
     2.5668000000000000,  0.7812000000000000,  0.0000000000000000,  0.0000000000000000, /* P(2,2) */
     0.432820000000000,   0.3001400000000000, -1.0000000000000000,  1.0000000000000000, /* P(3,0) */
    -0.3857193648237012,  2.065159951093377,   0.0000000000000000,  0.0000000000000000, /* P(3,1) */
    -4.8769200000000000,  3.359160000000000,   0.0000000000000000,  0.0000000000000000, /* P(3,2) */
     11.87127656630069,   1.993208492456322,   0.0000000000000000,  0.0000000000000000  /* P(3,3) */
  };
  const double dP_none[] = {
     0.0000000000000000,  0.0000000000000000,  0.0000000000000000,  0.0000000000000000, /* d/dt P(0,0) */
    -0.9249864863877740, -0.5102940328869230,  0.0000000000000000,  0.0000000000000000, /* d/dt P(1,0) */
    -0.3800000000000000,  0.8600000000000000, -1.0000000000000000,  1.0000000000000000, /* d/dt P(1,1) */
     1.0544845944820630, -1.3165586048482610,  0.0000000000000000,  0.0000000000000000, /* d/dt P(2,0) */
    -2.1336000000000000,  1.4376000000000000,  3.0000000000000000,  3.0000000000000000, /* d/dt P(2,1) */
    -2.108969188964126,   2.633117209696523,   0.0000000000000000,  0.0000000000000000, /* d/dt P(2,2) */
     0.3857193648237012, -2.0651599510933770,  0.0000000000000000,  0.0000000000000000, /* d/dt P(3,0) */
     5.0353800000000000,  0.1212600000000000, -6.0000000000000000,  6.0000000000000000, /* d/dt P(3,1) */
    -7.864235107268853,   9.32919550923872,    0.0000000000000000,  0.0000000000000000, /* d/dt P(3,2) */
    -14.630760000000000,  10.077480000000000,  0.0000000000000000,  0.0000000000000000  /* d/dt P(3,3) */
  };
  const double d2P_none[] = {
     0.0000000000000000,  0.0000000000000000,  0.0000000000000000,  0.0000000000000000, /* d2/dt2 P(0,0) */
     0.3800000000000000, -0.8600000000000000,  1.0000000000000000, -1.0000000000000000, /* d2/dt2 P(1,0) */
    -0.924986486387774,  -0.510294032886923,   0.0000000000000000,  0.0000000000000000, /* d2/dt2 P(1,1) */
     2.1336000000000000, -1.4376000000000000, -3.0000000000000000, -3.0000000000000000, /* d2/dt2 P(2,0) */
     4.217938377928252,  -5.266234419393045,   0.0000000000000000,  0.0000000000000000, /* d2/dt2 P(2,1) */
    -4.267200000000000,   2.875200000000000,   6.0000000000000000,  6.0000000000000000, /* d2/dt2 P(2,2) */
    -5.0353800000000000, -0.1212600000000000,  6.0000000000000000, -6.0000000000000000, /* d2/dt2 P(3,0) */
     6.246433742576635,  -17.05555746117962,   0.0000000000000000,  0.0000000000000000, /* d2/dt2 P(3,1) */
     32.49228000000000,  -4.432440000000000,  -30.000000000000000,  30.000000000000000, /* d2/dt2 P(3,2) */
    -23.59270532180656,   27.98758652771618,   0.0000000000000000,  0.0000000000000000  /* d2/dt2 P(3,3) */
  };

  /* test specific values */
  gsl_sf_legendre_precompute(norm, lmax, flags, Plm);

  if (norm == GSL_SF_LEGENDRE_SCHMIDT)
    {
      P_norm = P_schmidt;
      dP_norm = dP_schmidt;
      d2P_norm = d2P_schmidt;
    }
  else if (norm == GSL_SF_LEGENDRE_SPHARM)
    {
      P_norm = P_spharm;
      dP_norm = dP_spharm;
      d2P_norm = d2P_spharm;
    }
  else if (norm == GSL_SF_LEGENDRE_FULL)
    {
      P_norm = P_full;
      dP_norm = dP_full;
      d2P_norm = d2P_full;
    }
  else if (norm == GSL_SF_LEGENDRE_NONE)
    {
      P_norm = P_none;
      dP_norm = dP_none;
      d2P_norm = d2P_none;
    }
  else
    {
      fprintf(stderr, "test_legendre_eps: invalid norm\n");
      return -1;
    }

  for (i = 0; i < nx; ++i)
    {
      gsl_sf_legendre_deriv2_alt_arrayx(norm, lmax, x_vals[i], Plm, dPlm, d2Plm);

      for (l = 0; l <= 3; ++l)
        {
          size_t m;

          for (m = 0; m <= l; ++m)
            {
              size_t lmidx = gsl_sf_legendre_array_index(l, m);
              size_t idx = nx * lmidx + i;
              double csfac = 1.0;

              if ((flags & GSL_SF_LEGENDRE_FLG_CSPHASE) && (m & 1))
                csfac = -1.0;

              sprintf(buf, "x=%g", x_vals[i]);
              test_value(flags, lmax, l, m, Plm, csfac * P_norm[idx], tol, desc, buf);

              sprintf(buf, "deriv x=%g", x_vals[i]);
              test_value(flags, lmax, l, m, dPlm, csfac * dP_norm[idx], tol, desc, buf);

              sprintf(buf, "deriv2 x=%g", x_vals[i]);
              test_value(flags, lmax, l, m, d2Plm, csfac * d2P_norm[idx], tol, desc, buf);
            }
        }
    }

  /* test array routines */

  gsl_sf_legendre_precompute(norm, lmax, flags, Plm);
  gsl_sf_legendre_precompute(norm, lmax, flags, Plm_alt);

  for (x = -1.0 + dx; x < 1.0; x += dx)
    {
      double u = sqrt((1.0 + x) * (1.0 - x)); /* sin(theta) */

      s += gsl_sf_legendre_deriv2_arrayx(norm, lmax, x, Plm, dPlm, d2Plm);
      s += gsl_sf_legendre_deriv2_alt_arrayx(norm, lmax, x, Plm_alt, dPlm_alt, d2Plm_alt);

      for (l = 0; l <= lmax; ++l)
        {
          double sum = test_legendre_sum(indexl, lmax, l, Plm);
          double rhs;
          size_t m;

          /* test arrays vs alternate arrays */
          for (m = 0; m <= l; ++m)
            {
              size_t idx = indexl ? gsl_sf_legendre_array_index(l, m) : gsl_sf_legendre_array_index_m(l, m, lmax);
              gsl_test_rel(Plm[idx], Plm_alt[idx], tol, "%s Plm lmax=%zu l=%zu m=%zu x=%g", desc, lmax, l, m, x);
              gsl_test_rel(-u * dPlm[idx], dPlm_alt[idx], tol, "%s dPlm lmax=%zu l=%zu m=%zu x=%g", desc, lmax, l, m, x);

              if (fabs(d2Plm_alt[idx]) > 1.0e-10)
                {
                  gsl_test_rel(u*u * d2Plm[idx] - x * dPlm[idx], d2Plm_alt[idx], tol,
                               "%s d2Plm lmax=%zu l=%zu m=%zu x=%g", desc, lmax, l, m, x);
                }
            }

          if (norm == GSL_SF_LEGENDRE_SCHMIDT)
            rhs = 1.0;
          else if (norm == GSL_SF_LEGENDRE_SPHARM)
            {
              if (l == 0)
                rhs = (2.0 * l + 1.0) / (4.0 * M_PI);
              else
                {
                  double Pl = gsl_sf_legendre_Pl(l, x);
                  rhs = (2.0 * l + 1.0) / (8.0 * M_PI) * (1.0 + Pl*Pl);
                }
            }
          else if (norm == GSL_SF_LEGENDRE_FULL)
            {
              if (l == 0)
                rhs = 0.5 * (2.0 * l + 1.0);
              else
                {
                  double Pl = gsl_sf_legendre_Pl(l, x);
                  rhs = 0.25 * (2.0 * l + 1.0) * (1.0 + Pl*Pl);
                }
            }

          if (norm != GSL_SF_LEGENDRE_NONE)
            gsl_test_rel(sum, rhs, tol, "%s lmax=%zu l=%zu, x=%f, sum=%.12e", desc, lmax, l, x, sum);
        }
    }

  free(Plm);
  free(dPlm);
  free(d2Plm);
  free(Plm_alt);
  free(dPlm_alt);
  free(d2Plm_alt);

  return s;
}

static int
test_legendre_all(const size_t lmax)
{
  int s = 0;
  const double tol = 1.0e-10;
  size_t flags;

  flags = GSL_SF_LEGENDRE_FLG_INDEXL;
  s += test_legendre_eps(tol, GSL_SF_LEGENDRE_SCHMIDT, lmax, flags, "schmidt L nocsphase");
  s += test_legendre_eps(tol, GSL_SF_LEGENDRE_SPHARM, lmax, flags, "spharm L nocsphase");
  s += test_legendre_eps(tol, GSL_SF_LEGENDRE_FULL, lmax, flags, "full L nocsphase");
  if (lmax <= 100)
    {
      s += test_legendre_eps(tol, GSL_SF_LEGENDRE_NONE, lmax, flags, "unnorm L nocsphase");
      s += test_legendre_cross(tol, lmax, flags, "legendre cross L nocsphase");
    }

  flags = GSL_SF_LEGENDRE_FLG_INDEXL | GSL_SF_LEGENDRE_FLG_CSPHASE;
  s += test_legendre_eps(tol, GSL_SF_LEGENDRE_SCHMIDT, lmax, flags, "schmidt L csphase");
  s += test_legendre_eps(tol, GSL_SF_LEGENDRE_SPHARM, lmax, flags, "spharm L csphase");
  s += test_legendre_eps(tol, GSL_SF_LEGENDRE_FULL, lmax, flags, "full L csphase");
  if (lmax <= 100)
    {
      s += test_legendre_eps(tol, GSL_SF_LEGENDRE_NONE, lmax, flags, "unnorm L csphase");
      s += test_legendre_cross(tol, lmax, flags, "legendre cross L csphase");
    }

  flags = 0;
  s += test_legendre_eps(tol, GSL_SF_LEGENDRE_SCHMIDT, lmax, flags, "schmidt M nocsphase");
  s += test_legendre_eps(tol, GSL_SF_LEGENDRE_SPHARM, lmax, flags, "spharm M nocsphase");
  s += test_legendre_eps(tol, GSL_SF_LEGENDRE_FULL, lmax, flags, "full M nocsphase");
  if (lmax <= 100)
    {
      s += test_legendre_eps(tol, GSL_SF_LEGENDRE_NONE, lmax, flags, "unnorm M nocsphase");
      s += test_legendre_cross(tol, lmax, flags, "legendre cross M nocsphase");
    }

  flags = GSL_SF_LEGENDRE_FLG_CSPHASE;
  s += test_legendre_eps(tol, GSL_SF_LEGENDRE_SCHMIDT, lmax, flags, "schmidt M csphase");
  s += test_legendre_eps(tol, GSL_SF_LEGENDRE_SPHARM, lmax, flags, "spharm M csphase");
  s += test_legendre_eps(tol, GSL_SF_LEGENDRE_FULL, lmax, flags, "full M csphase");
  if (lmax <= 100)
    {
      s += test_legendre_eps(tol, GSL_SF_LEGENDRE_NONE, lmax, flags, "unnorm M csphase");
      s += test_legendre_cross(tol, lmax, flags, "legendre cross M csphase");
    }

  return s;
}

int
test_legendre(void)
{
  gsl_sf_result r;
  double L[256], DL[256];
  int s = 0;
  int sa;

  TEST_SF(s,  gsl_sf_legendre_P1_e, (-0.5, &r), -0.5, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_legendre_P1_e, ( 0.5, &r), 0.5, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_legendre_P2_e, (0.0, &r), -0.5  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_legendre_P2_e, (0.5, &r), -0.125, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_legendre_P2_e, (1.0, &r), 1.0  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_legendre_P2_e, (100.0, &r), 14999.5  , TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s,  gsl_sf_legendre_P3_e, ( -0.5, &r), 0.4375, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_legendre_P3_e, (  0.5, &r), -0.4375, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_legendre_P3_e, (  1.0, &r), 1.0        , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s,  gsl_sf_legendre_P3_e, (100.0, &r), 2.49985e+06, TEST_TOL0, GSL_SUCCESS);


  TEST_SF(s, gsl_sf_legendre_Pl_e, (1, -0.5, &r), -0.5, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Pl_e, (1,  1.0e-8, &r), 1.0e-08, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Pl_e, (1,  0.5, &r), 0.5, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Pl_e, (1,  1.0, &r), 1.0, TEST_TOL0, GSL_SUCCESS);
 
  TEST_SF(s, gsl_sf_legendre_Pl_e, (10, -0.5, &r), -0.1882286071777345, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Pl_e, (10,  1.0e-8, &r), -0.24609374999999864648, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Pl_e, (10,  0.5, &r), -0.18822860717773437500, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Pl_e, (10,  1.0, &r), 1.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_Pl_e, (99, -0.5, &r), 0.08300778172138770477, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Pl_e, (99,  1.0e-8, &r), -7.958923738716563193e-08, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Pl_e, (99,  0.5, &r), -0.08300778172138770477, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Pl_e, (99,  0.999, &r), -0.3317727359254778874, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Pl_e, (99,  1.0, &r), 1.0, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_Pl_e, (1000, -0.5, &r),   -0.019168251091650277878, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Pl_e, (1000,  1.0e-8, &r), 0.0252250181770982897470252620,  TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Pl_e, (1000,  0.5, &r),   -0.019168251091650277878, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Pl_e, (1000,  1.0, &r),    1.0,                     TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_Pl_e, (4000, -0.5, &r), -0.009585404456573080972, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Pl_e, (4000,  0.5, &r), -0.009585404456573080972, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Pl_e, (4000,  1.0, &r), 1.0, TEST_TOL0, GSL_SUCCESS);

  sa = 0;
  gsl_sf_legendre_Pl_array(100, 0.5, L);
  TEST_SF_VAL(sa, L[0],   +0.0,   1.0, TEST_TOL1);
  TEST_SF_VAL(sa, L[10],  +0.0,  -0.18822860717773437500, TEST_TOL1);
  TEST_SF_VAL(sa, L[100], +0.0,  -0.06051802596186118687, TEST_TOL1);
  gsl_test(sa, "gsl_sf_legendre_Pl_array(100, 0.5)");
  s += sa;

  sa = 0;
  gsl_sf_legendre_Pl_deriv_array(100, 0.5, L, DL);
  TEST_SF_VAL(sa, DL[0],   +0.0,   0.0, TEST_TOL1);
  TEST_SF_VAL(sa, DL[1],   +0.0,   1.0, TEST_TOL1);
  TEST_SF_VAL(sa, DL[10],  +0.0,  -2.3171234130859375000, TEST_TOL1);
  TEST_SF_VAL(sa, DL[100], +0.0,  -7.0331691653942815112, TEST_TOL1);
  gsl_test(sa, "gsl_sf_legendre_Pl_deriv_array(100, 0.5)");
  s += sa;
  sa = 0;

  gsl_sf_legendre_Pl_deriv_array(10, 1.0, L, DL);
  TEST_SF_VAL(sa, DL[0],   +0.0,   0.0, TEST_TOL1);
  TEST_SF_VAL(sa, DL[1],   +0.0,   1.0, TEST_TOL1);
  TEST_SF_VAL(sa, DL[10],  +0.0,  55.0, TEST_TOL1);
  gsl_test(sa, "gsl_sf_legendre_Pl_deriv_array(10, 1.0)");
  s += sa;

  gsl_sf_legendre_Pl_deriv_array(10, 1.0 - 1.0e-11, L, DL);
  TEST_SF_VAL(sa, DL[0],   +0.0,   0.0, TEST_TOL1);
  TEST_SF_VAL(sa, DL[1],   +0.0,   1.0, TEST_TOL1);
  TEST_SF_VAL(sa, DL[10],  +0.0,  54.999999985150000001, TEST_TOL1);
  gsl_test(sa, "gsl_sf_legendre_Pl_deriv_array(10, 1.0 - 1.0e-11)");
  s += sa;

  gsl_sf_legendre_Pl_deriv_array(10, -1.0, L, DL);
  TEST_SF_VAL(sa, DL[0],   +0.0,   0.0, TEST_TOL1);
  TEST_SF_VAL(sa, DL[1],   +0.0,   1.0, TEST_TOL1);
  TEST_SF_VAL(sa, DL[10],  +0.0, -55.0, TEST_TOL1);
  gsl_test(sa, "gsl_sf_legendre_Pl_deriv_array(10, -1.0)");
  s += sa;

  gsl_sf_legendre_Pl_deriv_array(10, -1.0 + 1.0e-11, L, DL);
  TEST_SF_VAL(sa, DL[0],   +0.0,   0.0, TEST_TOL1);
  TEST_SF_VAL(sa, DL[1],   +0.0,   1.0, TEST_TOL1);
  TEST_SF_VAL(sa, DL[10],  +0.0, -54.999999985150000001, TEST_TOL1);
  gsl_test(sa, "gsl_sf_legendre_Pl_deriv_array(10, -1.0 + 1.0e-11)");
  s += sa;

  TEST_SF(s, gsl_sf_legendre_Plm_e, (10, 0, -0.5, &r), -0.18822860717773437500, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Plm_e, (10, 0, 1.0e-08, &r), -0.24609374999999864648, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Plm_e, (10, 0, 0.5, &r), -0.18822860717773437500, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_Plm_e, (10, 1, -0.5, &r), -2.0066877394361256516, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Plm_e, (10, 1, 1.0e-08, &r), -2.7070312499999951725e-07, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Plm_e, (10, 1, 0.5, &r), 2.0066877394361256516, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_Plm_e, (10, 5, -0.5, &r),    -30086.169706116174977,    TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Plm_e, (10, 5, 1.0e-08, &r), -0.0025337812499999964949, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Plm_e, (10, 5, 0.5, &r),      30086.169706116174977,    TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Plm_e, (10, 5, 0.999, &r),   -0.5036411489013270406,    TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_Plm_e, (100, 5, -0.5, &r), -6.617107444248382171e+08, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Plm_e, (100, 5, 1.0e-08, &r), 817.8987598063712851, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Plm_e, (100, 5, 0.5, &r), 6.617107444248382171e+08, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Plm_e, (100, 5, 0.999, &r), -1.9831610803806212189e+09, TEST_TOL2, GSL_SUCCESS);

#ifndef GSL_DISABLE_DEPRECATED

  sa = 0;
  gsl_sf_legendre_Plm_deriv_array(100, 2, -1.0 + 1.0/1125899906842624.0, L, DL);
  TEST_SF_VAL(sa, L[0],   +0.0,  5.3290705182007490275e-15, TEST_TOL1);
  TEST_SF_VAL(sa, L[1],   +0.0, -2.6645352591003721471e-14, TEST_TOL1);
  TEST_SF_VAL(sa, L[98],  +0.0,  2.2646284847349109694e-08, TEST_TOL2);
  gsl_test(sa, "gsl_sf_legendre_Plm_deriv_array(100, 2, -1.0 + 2^(-50)");
  s += sa;

  sa = 0;
  gsl_sf_legendre_Plm_deriv_array(100, 2, 1.0 - 1.0/1125899906842624.0, L, DL);
  TEST_SF_VAL(sa, L[0],   +0.0,  5.3290705182007490275e-15, TEST_TOL1);
  TEST_SF_VAL(sa, L[1],   +0.0,  2.6645352591003721471e-14, TEST_TOL1);
  TEST_SF_VAL(sa, L[10],  +0.0,  5.3343995887188313290e-12, TEST_TOL1);
  TEST_SF_VAL(sa, L[98],  +0.0,  2.2646284847349109694e-08, TEST_TOL2);
  gsl_test(sa, "gsl_sf_legendre_Plm_deriv_array(100, 2, 1.0 - 2^(-50)");
  s += sa;

  sa = 0;
  gsl_sf_legendre_Plm_array(100, 5, 0.5, L);
  TEST_SF_VAL(sa, L[0],  +0.0, -460.3466286991656682, TEST_TOL1);
  TEST_SF_VAL(sa, L[10], +0.0,  38852.51334152290535, TEST_TOL1 );
  TEST_SF_VAL(sa, L[95], +0.0,  6.617107444248382171e+08, TEST_TOL1);
  gsl_test(sa, "gsl_sf_legendre_Plm_array(100, 5, 0.5)");
  s += sa;

  sa = 0;
  gsl_sf_legendre_Plm_array(100, 5, 0.999, L);
  TEST_SF_VAL(sa, L[0],  +0.0,  -0.00016883550990916552255, TEST_TOL2);
  TEST_SF_VAL(sa, L[10], +0.0,  -30.651334850159821525, TEST_TOL2 );
  TEST_SF_VAL(sa, L[95], +0.0,  -1.9831610803806212189e+09, TEST_TOL2);
  gsl_test(sa, "gsl_sf_legendre_Plm_array(100, 5, 0.999)");
  s += sa;

  sa = 0;
  gsl_sf_legendre_Plm_array(100, 5, -0.999, L);
  TEST_SF_VAL(sa, L[0],  +0.0,  -0.00016883550990916552255, TEST_TOL2);
  TEST_SF_VAL(sa, L[10], +0.0,  -30.651334850159821525, TEST_TOL2 );
  TEST_SF_VAL(sa, L[95], +0.0,   1.9831610803806212189e+09, TEST_TOL2);
  gsl_test(sa, "gsl_sf_legendre_Plm_array(100, 5, -0.999)");
  s += sa;

  sa = 0;
  gsl_sf_legendre_Plm_deriv_array(100, 2, 0.999, L, DL);
  TEST_SF_VAL(sa, L[0],   +0.0,  0.00599700000000000000, TEST_TOL1);
  TEST_SF_VAL(sa, L[1],   +0.0,  0.02995501500000000000, TEST_TOL1 );
  TEST_SF_VAL(sa, DL[0],  +0.0,  -5.9940000000000000000, TEST_TOL1);
  TEST_SF_VAL(sa, DL[1],  +0.0,  -29.910045000000000000, TEST_TOL1 );
  TEST_SF_VAL(sa, DL[2],  +0.0,  -89.490629790000000000, TEST_TOL1 );
  TEST_SF_VAL(sa, DL[10], +0.0,  -5703.9461633355291972, TEST_TOL1 );
  TEST_SF_VAL(sa, DL[95], +0.0,  6.4518473603456858414E+06, TEST_TOL3);
  gsl_test(sa, "gsl_sf_legendre_Plm_deriv_array(100, 2, 0.999)");
  s += sa;

  sa = 0;
  gsl_sf_legendre_Plm_deriv_array(100, 2, 1.0 - 1.0e-15, L, DL);
  TEST_SF_VAL(sa, DL[0],  +0.0,  -5.9999999999999940000, TEST_TOL1);
  TEST_SF_VAL(sa, DL[1],  +0.0,  -29.999999999999910000, TEST_TOL1 );
  TEST_SF_VAL(sa, DL[2],  +0.0,  -89.999999999999490000, TEST_TOL1 );
  TEST_SF_VAL(sa, DL[10], +0.0,  -6005.9999999996936940, TEST_TOL1 );
  TEST_SF_VAL(sa, DL[95], +0.0,  -2.2586255999928454270e+07, TEST_TOL3);
  gsl_test(sa, "gsl_sf_legendre_Plm_deriv_array(100, 2, 1.0 - 1.0e-15)");
  s += sa;

  sa = 0;
  gsl_sf_legendre_Plm_deriv_array(100, 2, -1.0 + 1.0e-15, L, DL);
  TEST_SF_VAL(sa, DL[0],  +0.0,   5.9999999999999940000, TEST_TOL1);
  TEST_SF_VAL(sa, DL[1],  +0.0,  -29.999999999999910000, TEST_TOL1 );
  TEST_SF_VAL(sa, DL[95], +0.0,  -2.2586255999928454270e+07, TEST_TOL3);
  gsl_test(sa, "gsl_sf_legendre_Plm_deriv_array(100, 2, -1.0 + 1.0e-15)");
  s += sa;

  sa = 0;
  gsl_sf_legendre_Plm_deriv_array(100, 5, 0.999, L, DL);
  TEST_SF_VAL(sa, DL[0],  +0.0,  0.42187762481054616565, TEST_TOL1);
  TEST_SF_VAL(sa, DL[1],  +0.0,  4.6341560284340909936, TEST_TOL1 );
  TEST_SF_VAL(sa, DL[2],  +0.0,  27.759505566959219127, TEST_TOL1 );
  TEST_SF_VAL(sa, DL[10], +0.0,  76051.795860179545484, TEST_TOL1 );
  TEST_SF_VAL(sa, DL[95], +0.0,  3.0344503083851936814e+12, TEST_TOL3);
  gsl_test(sa, "gsl_sf_legendre_Plm_deriv_array(100, 5, 0.999)");
  s += sa;

#endif /* !GSL_DISABLE_DEPRECATED */

  TEST_SF(s, gsl_sf_legendre_sphPlm_e, (10, 0, -0.5, &r), -0.24332702369300133776, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_sphPlm_e, (10, 0, 0.5, &r), -0.24332702369300133776, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_sphPlm_e, (10, 0, 0.999, &r), 1.2225754122797385990, TEST_TOL1, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_sphPlm_e, (10, 5, -0.5, &r),    -0.3725739049803293972,     TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_sphPlm_e, (10, 5, 1.0e-08, &r), -3.1377233589376792243e-08, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_sphPlm_e, (10, 5, 0.5, &r),      0.3725739049803293972,     TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_sphPlm_e, (10, 5, 0.999, &r),   -6.236870674727370094e-06,  TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_sphPlm_e, (10, 10, -0.5, &r), 0.12876871185785724117, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_sphPlm_e, (10, 10, 0.5, &r), 0.12876871185785724117,  TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_sphPlm_e, (10, 10, 0.999, &r), 1.7320802307583118647e-14, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_sphPlm_e, (200, 1, -0.5, &r),   0.3302975570099492931, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_sphPlm_e, (200, 1, 0.5, &r),   -0.3302975570099492931, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_sphPlm_e, (200, 1, 0.999, &r), -1.4069792055546256912, TEST_TOL2, GSL_SUCCESS);



  /* Test case from alberto@physik.fu-berlin.de */

  TEST_SF(s, gsl_sf_legendre_sphPlm_e, (3, 1, 0.0, &r), 0.323180184114150653007, TEST_TOL2, GSL_SUCCESS);

  /* Other test cases */

  TEST_SF(s, gsl_sf_legendre_sphPlm_e, (200, 1, -0.5, &r), 0.3302975570099492931418227583, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_sphPlm_e, (140,135,1,&r), 0.0, TEST_TOL2, GSL_SUCCESS);

#ifdef EXTENDED
  TEST_SF(s, gsl_sf_legendre_sphPlm_e, (140,135,0.99998689456491752,&r), -6.54265253269093276310395668335e-305, TEST_TOL6, GSL_SUCCESS);
#endif

#ifndef GSL_DISABLE_DEPRECATED

  sa = 0;
  gsl_sf_legendre_sphPlm_array(100, 5, 0.5, L);
  TEST_SF_VAL(sa, L[0],  +0.0, -0.22609703187800460722, TEST_TOL1);
  TEST_SF_VAL(sa, L[10], +0.0,  0.07452710323813558940, TEST_TOL1);
  TEST_SF_VAL(sa, L[95], +0.0,  0.25865355990880161717, TEST_TOL1);
  gsl_test(sa, "gsl_sf_legendre_sphPlm_array(100, 5, 0.5)");
  s += sa;

  sa = 0;
  gsl_sf_legendre_sphPlm_array(100, 2, 1.0 - 1.0/1125899906842624.0, L);
  TEST_SF_VAL(sa, L[0],  +0.0, 6.8616082064776657177e-16, TEST_TOL2);
  TEST_SF_VAL(sa, L[10], +0.0, 4.8543150313086787324e-14, TEST_TOL2);
  TEST_SF_VAL(sa, L[95], +0.0, 8.3138984963650838973e-12, TEST_TOL2);
  gsl_test(sa, "gsl_sf_legendre_sphPlm_array(100, 2, 1.0 - 2^(-50))");
  s += sa;

  sa = 0;
  gsl_sf_legendre_sphPlm_array(100, 2, -1.0 + 1.0/1125899906842624.0, L);
  TEST_SF_VAL(sa, L[0],  +0.0,  6.8616082064776657177e-16, TEST_TOL2);
  TEST_SF_VAL(sa, L[95], +0.0, -8.3138984963650838973e-12, TEST_TOL2);
  gsl_test(sa, "gsl_sf_legendre_sphPlm_array(100, 2, -1.0 + 2^(-50))");
  s += sa;

  sa = 0;
  gsl_sf_legendre_sphPlm_deriv_array(100, 0, 0.5, L, DL);
  TEST_SF_VAL(sa, DL[0],  +0.0,  0.0, TEST_TOL1);
  TEST_SF_VAL(sa, DL[10], +0.0, -2.9953934850252897591, TEST_TOL1);
  TEST_SF_VAL(sa, DL[95], +0.0, -36.411811015111761007, TEST_TOL1);
  gsl_test(sa, "gsl_sf_legendre_sphPlm_deriv_array(100, 0, 0.5)");
  s += sa;

  sa = 0;
  gsl_sf_legendre_sphPlm_deriv_array(100, 1, 0.5, L, DL);
  TEST_SF_VAL(sa, DL[0],  +0.0,  0.19947114020071633897, TEST_TOL1);
  TEST_SF_VAL(sa, DL[1],  +0.0, -0.44603102903819277863, TEST_TOL1);
  TEST_SF_VAL(sa, DL[10], +0.0,  1.3658895325030216565, TEST_TOL1);
  TEST_SF_VAL(sa, DL[99], +0.0, -27.925571865639037118, TEST_TOL1);
  gsl_test(sa, "gsl_sf_legendre_sphPlm_deriv_array(100, 1, 0.5)");
  s += sa;

  sa = 0;
  gsl_sf_legendre_sphPlm_deriv_array(100, 1, 1.0 - 1.0/1125899906842624.0, L, DL);
  TEST_SF_VAL(sa, DL[0],  +0.0,  8.1973898803378530946e+06, TEST_TOL1);
  TEST_SF_VAL(sa, DL[1],  +0.0,  1.8329921010504257405e+07, TEST_TOL1);
  TEST_SF_VAL(sa, DL[10], +0.0,  1.8439572562895384115e+08, TEST_TOL1);
  TEST_SF_VAL(sa, DL[99], +0.0,  4.7682463136232210552e+09, TEST_TOL3);
  gsl_test(sa, "gsl_sf_legendre_sphPlm_deriv_array(100, 1, 1.0 - 2^(-50))");
  s += sa;

  sa = 0;
  gsl_sf_legendre_sphPlm_deriv_array(100, 2, 0.5, L, DL);
  TEST_SF_VAL(sa, DL[0],  +0.0, -0.38627420202318958034, TEST_TOL1);
  TEST_SF_VAL(sa, DL[1],  +0.0,  0.25549636910832059085, TEST_TOL1);
  TEST_SF_VAL(sa, DL[2],  +0.0,  1.5053547230039006279, TEST_TOL1);
  TEST_SF_VAL(sa, DL[10], +0.0,  0.73576559668648243477, TEST_TOL1);
  TEST_SF_VAL(sa, DL[98], +0.0, 28.444589950264378407, TEST_TOL1);
  gsl_test(sa, "gsl_sf_legendre_sphPlm_deriv_array(100, 2, 0.5)");
  s += sa;

  sa = 0;
  gsl_sf_legendre_sphPlm_deriv_array(100, 5, 0.5, L, DL);
  TEST_SF_VAL(sa, DL[0],  +0.0,  0.75365677292668202407, TEST_TOL1);
  TEST_SF_VAL(sa, DL[1],  +0.0,  0.54346962777757450534, TEST_TOL1);
  TEST_SF_VAL(sa, DL[2],  +0.0, -0.98309969029001383773, TEST_TOL1);
  TEST_SF_VAL(sa, DL[3],  +0.0, -2.7728270988954534293, TEST_TOL1);
  TEST_SF_VAL(sa, DL[10], +0.0, -5.7407133315443482193, TEST_TOL1);
  TEST_SF_VAL(sa, DL[95], +0.0, -25.893934624747394561, TEST_TOL1);
  gsl_test(sa, "gsl_sf_legendre_sphPlm_deriv_array(100, 5, 0.5)");
  s += sa;
  sa = 0;

  gsl_sf_legendre_sphPlm_deriv_array(100, 5, 1.0 - 1.0/1125899906842624.0, L, DL);
  TEST_SF_VAL(sa, DL[0],  +0.0,  1.7374288379067753301e-22, TEST_TOL1);
  TEST_SF_VAL(sa, DL[1],  +0.0,  6.2643887625426827113e-22, TEST_TOL1);
  TEST_SF_VAL(sa, DL[2],  +0.0,  1.6482697200734667281e-21, TEST_TOL1);
  TEST_SF_VAL(sa, DL[95], +0.0,  3.9890549466071349506e-15, TEST_TOL2);
  gsl_test(sa, "gsl_sf_legendre_sphPlm_deriv_array(100, 5, 1.0 - 2^(-50))");
  s += sa;

  gsl_sf_legendre_sphPlm_deriv_array(100, 5, -1.0 + 1.0/1125899906842624.0, L, DL);
  TEST_SF_VAL(sa, DL[0],  +0.0, -1.7374288379067753301e-22, TEST_TOL1);
  TEST_SF_VAL(sa, DL[1],  +0.0,  6.2643887625426827113e-22, TEST_TOL1);
  TEST_SF_VAL(sa, DL[2],  +0.0, -1.6482697200734667281e-21, TEST_TOL1);
  TEST_SF_VAL(sa, DL[95], +0.0,  3.9890549466071349506e-15, TEST_TOL3);
  gsl_test(sa, "gsl_sf_legendre_sphPlm_deriv_array(100, 5, -1.0 + 2^(-50))");
  s += sa;

#endif /* !GSL_DISABLE_DEPRECATED */

  TEST_SF(s, gsl_sf_conicalP_half_e, (0.0, -0.5, &r),   0.8573827581049917129, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_half_e, (0.0,  0.5, &r),   0.8573827581049917129, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_half_e, (0.0,  2.0, &r),   0.6062611623284649811, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_half_e, (0.0,  100.0, &r), 0.07979045091636735635, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_conicalP_half_e, (10.0, -0.5, &r),    5.345484922591867188e+08, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_half_e, (10.0,  0.5, &r),    15137.910380385258370, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_half_e, (10.0,  2.0, &r),    0.4992680691891618544, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_half_e, (10.0,  100.0, &r), -0.07272008163718195685, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_conicalP_half_e, (200.0, -1.0e-3, &r),  1.3347639529084185010e+136, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_half_e, (200.0,  1.0e-8, &r),  1.0928098010940058507e+136, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_half_e, (200.0,  0.5, &r),     3.895546021611205442e+90,   TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_half_e, (200.0,  10.0, &r),   -0.04308567180833581268,     TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_half_e, (200.0,  100.0, &r),  -0.04694669186576399194,     TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_half_e, (200.0,  1000.0, &r),  0.023698140704121273277,    TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_half_e, (200.0,  1.0e+8, &r), -0.00006790983312124277891,  TEST_TOL3, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_conicalP_half_e, (1.0e+8,  1.1, &r),   1.1599311133054742944,  TEST_SQRT_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_half_e, (1.0e+8,  100.0, &r), 0.07971967557381557875, TEST_SQRT_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_conicalP_mhalf_e, (0.0, -0.5, &r),  1.7956982494514644808, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_mhalf_e, (0.0,  0.5, &r),  0.8978491247257322404, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_mhalf_e, (0.0,  2.0, &r),  0.7984204253272901551, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_mhalf_e, (0.0,  100.0, &r),  0.4227531369388072584, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_conicalP_mhalf_e, (10.0, -0.5, &r),  5.345484922591867181e+07, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_mhalf_e, (10.0,  0.5, &r),  1513.7910356104985334, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_mhalf_e, (10.0,  2.0, &r),  0.03439243987215615642, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_mhalf_e, (10.0,  100.0, &r),  0.003283756665952609624, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_conicalP_mhalf_e, (200.0, -0.5, &r),  1.7699538115312304280e+179, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_mhalf_e, (200.0,  1.0e-8, &r),  5.464049005470029253e+133, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_mhalf_e, (200.0,  0.5, &r),  1.9477730108056027211e+88, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_mhalf_e, (200.0,  10.0, &r),  0.0012462575917716355362, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_mhalf_e, (200.0,  100.0, &r),  -0.0003225881344802625149, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_mhalf_e, (200.0,  1000.0, &r), -0.00004330652890886567623, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_mhalf_e, (200.0,  1.0e+8, &r),  2.0943091278037078483e-07, TEST_TOL3, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_conicalP_mhalf_e, (1.0e+8,  1.1, &r), 2.092320445620989618e-09, 16.0*TEST_SQRT_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_mhalf_e, (1.0e+8,  100.0, &r),  -3.359967833599016923e-11, 256.0*TEST_SQRT_TOL0, GSL_SUCCESS);


  TEST_SF(s, gsl_sf_conicalP_0_e, (0.0, -0.5, &r),  1.3728805006183501647, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_0_e, (0.0,  0.5, &r),  1.0731820071493643751, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_0_e, (0.0,  2.0, &r),  0.9012862993604472987, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_0_e, (0.0,  100.0, &r),  0.30091748588199264556, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_conicalP_0_e, (10.0, -0.5, &r),  1.6795592815421804669e+08, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_0_e, (10.0,  0.5, &r),  4826.034132009618240,      TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_0_e, (10.0,  2.0, &r),  0.18798468917758716146,    TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_0_e, (10.0,  100.0, &r), -0.008622130749987962529, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_conicalP_0_e, (200.0,  -0.5, &r), 2.502194818646823e+180, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_conicalP_0_e, (1000.0,  100.0, &r),   0.0017908817653497715844, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_0_e, (1000.0,  1000.0, &r), -0.0006566893804926284301, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_0_e, (1000.0,  1.0e+8, &r),  2.3167213561756390068e-06, TEST_TOL4, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_conicalP_1_e, (0.0, -0.5, &r),    0.4939371126656998499,  TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_1_e, (0.0,  0.5, &r),    0.14933621085538265636, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_1_e, (0.0,  2.0, &r),   -0.13666874968871549533, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_1_e, (0.0,  100.0, &r), -0.10544528203156629098, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_conicalP_1_e, (10.0, -0.5, &r),    1.7253802958788312520e+09, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_1_e, (10.0,  0.5, &r),    46781.02294059967988,      TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_1_e, (10.0,  2.0, &r),    0.26613342643657444400,    TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_1_e, (10.0,  100.0, &r), -0.23281959695501029796,    TEST_TOL2, GSL_SUCCESS);


  /* FIXME: Mathematica gets some brain-damaged numbers for
   * these x < 0 points. I have checked what I am doing in detail,
   * and it must be right because you can do it by summing
   * manifestly positive definite quantities.
   */
  TEST_SF(s, gsl_sf_conicalP_1_e, (200.0, -0.999, &r), 2.71635193199341135e+270, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_1_e, (200.0, -0.9, &r),   4.2952493176812905e+234,  TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_1_e, (200.0, -0.5, &r),   5.01159205956053439e+182, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_1_e, (200.0,  0.999, &r), 195733.0396081538,        TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_1_e, (200.0,  10.0, &r), -2.9272610662414349553,    TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_conicalP_1_e, (1000.0, 100.0, &r),  -1.7783258105862399857,    TEST_TOL6, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_1_e, (1000.0, 1000.0, &r),  0.4535161075156427179,    TEST_TOL4, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_1_e, (1000.0, 1.0e+8, &r),  0.0009983414549874888478, TEST_SQRT_TOL0, GSL_SUCCESS);


  TEST_SF(s, gsl_sf_conicalP_sph_reg_e, (2,  1.0, -0.5, &r),  1.6406279287008789526,      TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_sph_reg_e, (10, 1.0, -0.5, &r),  0.000029315266725049129448, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_sph_reg_e, (20, 1.0, -0.5, &r),  7.335769429462034431e-15,   TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_sph_reg_e, (30, 1.0, -0.5, &r),  1.3235612394267378871e-26,  TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_sph_reg_e, (10, 1.0, 0.5, &r),  2.7016087199857873954e-10, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_sph_reg_e, (20, 1.0, 0.5, &r),  1.1782569701435933399e-24, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_sph_reg_e, (30, 1.0, 0.5, &r),  3.636240588303797919e-41,  TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_sph_reg_e, (10, 1.0, 2.0, &r),  2.4934929626284934483e-10, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_sph_reg_e, (20, 1.0, 2.0, &r),  1.1284762488012616191e-24, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_sph_reg_e, (30, 100.0, 100.0, &r),  -1.6757772087159526048e-64, TEST_TOL6, GSL_SUCCESS);


  TEST_SF(s, gsl_sf_conicalP_cyl_reg_e, (2, 1.0, -0.5, &r),   2.2048510472375258708,       TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_cyl_reg_e, (10, 1.0, -0.5, &r),  0.00007335034531618655690,   TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_cyl_reg_e, (20, 1.0, -0.5, &r),  2.5419860619212164696e-14,   TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_cyl_reg_e, (30, 1.0, -0.5, &r),  5.579714972260536827e-26,    TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_cyl_reg_e, (10, 1.0, 0.5, &r),  1.1674078819646475282e-09,    TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_cyl_reg_e, (20, 1.0, 0.5, &r),  7.066408031229072207e-24,     TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_cyl_reg_e, (30, 1.0, 0.5, &r),  2.6541973286862588488e-40,    TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_cyl_reg_e, (10, 1.0, 2.0, &r),  1.0736109751890863051e-09,    TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_cyl_reg_e, (20, 1.0, 2.0, &r),  6.760965304863386741e-24,     TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_conicalP_cyl_reg_e, (30, 100.0, 100.0, &r), -4.268753482520651007e-63, TEST_TOL4, GSL_SUCCESS);


  TEST_SF(s, gsl_sf_legendre_H3d_0_e, (1.0e-06, 1.0e-06, &r), 0.9999999999998333333    , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_0_e, (1.0, 0.0, &r), 1.0                      , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_0_e, (1.0, 1.0, &r), 0.7160229153604338713    , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_0_e, (1.0, 100.0, &r), -3.767437313149604566e-44 , TEST_TOL2, GSL_SUCCESS);  
  TEST_SF(s, gsl_sf_legendre_H3d_0_e, (1.0, 500.0, &r), -6.665351935878582205e-218, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_0_e, (100.0, 1.0, &r), -0.004308757035378200029  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_0_e, (100.0, 10.0, &r), 7.508054627912986427e-07 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_0_e, (1000.0, 1.0, &r), 0.0007036067909088818319 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_0_e, (1.0e+08, 1.0, &r), 7.927485371429105968e-09 , TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_0_e, (1.0e+08, 100.0, &r), -3.627118904186918957e-52 , 32.0*TEST_SQRT_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_H3d_1_e, (1.0e-06, 1.0e-06, &r), 3.333333333334222222e-07, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_1_e, (1.0, 1.0e-10, &r), 4.714045207910316829e-11,     TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_1_e, (1.0, 1.0, &r), 0.3397013994799344639,           TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_1_e, (1.0, 100.0, &r), -7.200624449531811272e-44,     TEST_TOL2, GSL_SUCCESS);  
  TEST_SF(s, gsl_sf_legendre_H3d_1_e, (1.0, 500.0, &r), 4.192260336821728677e-218,     TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_1_e, (100.0, 0.01, &r), 0.30117664944267412324   , TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_1_e, (100.0, 1.0, &r), -0.007393833425336299309  , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_1_e, (100.0, 10.0, &r), -5.031062029821254982e-07 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_1_e, (1000.0, 0.001, &r), 0.30116875865090396421   , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_1_e, (1000.0, 1.0, &r), -0.0004776144516074971885 , TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_1_e, (1.0e+08, 1.0e-08, &r), 0.30116867893975679722   , TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_1_e, (1.0e+08, 1.0, &r), 3.0921097047369081582e-09, TEST_TOL4, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_1_e, (1.0e+08, 100.0, &r), -6.496142701296286936e-52 , 32.0*TEST_SQRT_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 1.0e-06, 1.0e-06, &r),  1.1544011544013627977e-32, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 1.0, 1.0e-10, &r),      2.0224912016958766992e-52, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 1.0, 1.0, &r),          0.011498635037491577728,   TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 1.0, 5.0, &r),          0.0020696945662545205776,  TEST_TOL4, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 1.0, 7.0, &r),     -0.0017555303787488993676,   TEST_TOL4, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 1.0, 10.0, &r),     0.00008999979724504887101,  TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 1.0, 100.0, &r),   -4.185397793298567945e-44,   TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 1.0, 500.0, &r),    1.4235113901091961263e-217, TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 100.0, 0.001, &r),  9.642762597222417946e-10,   TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 100.0, 0.002, &r),  3.0821201254308036109e-08,  TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 100.0, 0.01, &r),   0.00009281069019005840532,  TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 100.0, 1.0, &r),   -0.008043100696178624653,    TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 100.0, 10.0, &r),  -3.927678432813974207e-07,   TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 1000.0, 0.001, &r),  0.00009256365284253254503, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 1000.0, 0.01, &r),  -0.05553733815473079983, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 1.0e+08, 1.0e-08, &r),   0.00009256115861125841299, TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_H3d_e, (5, 1.0e+08, 100.0, &r),    -6.496143209092860765e-52 , 128.0*TEST_SQRT_TOL0, GSL_SUCCESS);

#if FIXME
  sa = 0;
  gsl_sf_legendre_H3d_array(100, 1.0, 3.0, L);
  TEST_SF_VAL(sa, L[0], +0.0, gsl_sf_legendre_H3d(0, 1.0, 3.0), 1.0e-12);
  TEST_SF_VAL(sa, L[1], +0.0, gsl_sf_legendre_H3d(1, 1.0, 3.0), 1.0e-12);
  TEST_SF_VAL(sa, L[10], +0.0, gsl_sf_legendre_H3d(10, 1.0, 3.0), 1.0e-12);
  TEST_SF_VAL(sa, L[100], +0.0, gsl_sf_legendre_H3d(100, 1.0, 3.0), 1.0e-12);
  gsl_test(sa, "  gsl_sf_legendre_H3d_array(100, 1.0, 3.0)");
  s += sa;
#endif

  /* x = -1 + 2^-16 */
  TEST_SF(s, gsl_sf_legendre_Q0_e, (-0.9999847412109375, &r), -5.8917472200477175158028143531855, TEST_TOL4, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Q0_e, (-0.5, &r), -0.5493061443340548457, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Q0_e, (-1e-10, &r), -1.000000000000000000e-10, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Q0_e, (0.0, &r), 0.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Q0_e, (1e-10, &r), 1.000000000000000000e-10, TEST_TOL0, GSL_SUCCESS);
  /* x = 1 - 2^-16 */
  TEST_SF(s, gsl_sf_legendre_Q0_e, (0.9999847412109375, &r), 5.8917472200477175158028143531855, TEST_TOL4, GSL_SUCCESS);
  /* x = 1 + 2^-16 */
  TEST_SF(s, gsl_sf_legendre_Q0_e, ( 1.0000152587890625, &r), 5.8917548494422489138325509750429, TEST_TOL4, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Q0_e, ( 1.5, &r), 0.8047189562170501873, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Q0_e, ( 9.99, &r), 0.1004364599660005447, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Q0_e, ( 10.0, &r), 0.1003353477310755806, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Q0_e, ( 10.01, &r), 0.1002344395571710243, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Q0_e, ( 100, &r), 0.010000333353334762015, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Q0_e, ( 1e10, &r), 1.000000000000000000e-10, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_Q1_e, (-0.9999847412109375, &r), 4.8916573191196772369, TEST_TOL4, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Q1_e, (-0.5, &r), -0.7253469278329725772, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Q1_e, (-0.01, &r), -0.9998999966664666524, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Q1_e, (-1e-10, &r), -0.999999999999999999, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Q1_e, (0.0, &r), -1.0, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Q1_e, (1e-10, &r), -0.999999999999999999, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Q1_e, (0.0001, &r), -0.9999999899999999667, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Q1_e, (0.01, &r), -0.9998999966664666524, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Q1_e, (0.5, &r), -0.7253469278329725772, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Q1_e, (0.9999847412109375, &r), 4.8916573191196772369, TEST_TOL4, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Q1_e, (1.0000152587890625, &r), 4.8918447504867045145, TEST_TOL4, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Q1_e, ( 1.5, &r), 0.20707843432557528095, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Q1_e, ( 9.99, &r), 3.360235060345441639e-3, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Q1_e, ( 10.0, &r), 3.353477310755806357e-3, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Q1_e, ( 10.01, &r), 3.346739967281953346e-3, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Q1_e, ( 100.0, &r), 3.333533347620158821e-5, TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Q1_e, ( 1e10, &r), 3.333333333333333333e-21, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_Ql_e, (10, -0.5, &r), -0.29165813966586752393,    TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Ql_e, (10,  0.5, &r), 0.29165813966586752393,     TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Ql_e, (10,  1.5, &r), 0.000014714232718207477406, TEST_TOL0, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_Ql_e, (100, -0.5, &r), -0.09492507395207282096,   TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Ql_e, (100,  0.5, &r), 0.09492507395207282096,    TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Ql_e, (100,  1.5, &r), 1.1628163435044121988e-43, TEST_TOL2, GSL_SUCCESS);

  TEST_SF(s, gsl_sf_legendre_Ql_e, (1000, -0.5, &r), -0.030105074974005303500, TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Ql_e, (1000,  0.5, &r), 0.030105074974005303500,  TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_legendre_Ql_e, (1000,  1.1, &r), 1.0757258447825356443e-194, TEST_TOL3, GSL_SUCCESS);

  /* test associated legendre functions */
  {
    size_t l;

    for (l = 0; l <= 10; ++l)
      test_legendre_all(l);

    test_legendre_all(140);
    test_legendre_all(1000);
    /*test_legendre_all(2700);*/
  }

  return s;
}
