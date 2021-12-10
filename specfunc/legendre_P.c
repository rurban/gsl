/* specfunc/legendre_P.c
 * 
 * Copyright (C) 2009-2021 Patrick Alken
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

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_legendre.h>

/*
 * The routines in this module compute associated Legendre functions
 * (ALFs) up to order and degree 2700, using the method described
 * in
 *
 * [1] S. A. Holmes and W. E. Featherstone, A unified approach
 *     to the Clenshaw summation and the recursive computation of very
 *     high degree and order normalised associated Legendre functions,
 *     Journal of Geodesy, 76, pg. 279-299, 2002.
 *
 * Computations of the alternative derivatives, d^k/dtheta^k Plm
 * use the methods described in,
 *
 * [2] W. Bosch, On the computation of derivatives of Legendre functions,
 *     Phys. Chem. Earth (A), 25, 9-11, pg. 655-659, 2000.
 *
 * Further information on ALFs can be found in
 *
 * [3] Abramowitz and Stegun, Handbook of Mathematical Functions,
 *     Chapter 8, 1972.
 */

static int legendre_derivk_alt_array(const gsl_sf_legendre_t norm,
                                     const size_t lmax, const double csfac,
                                     const double sqrts[], const double input_array[],
                                     double output_array[]);
static void legendre_sqrts(const size_t lmax, double *array);

#define LEGENDRE
#include "legendre_source.c"
#undef LEGENDRE

#define LEGENDRE_DERIV
#include "legendre_source.c"
#undef LEGENDRE_DERIV

#define LEGENDRE_DERIV_ALT
#include "legendre_source.c"
#undef LEGENDRE_DERIV_ALT

#define LEGENDRE_DERIV2
#include "legendre_source.c"
#undef LEGENDRE_DERIV2

#define LEGENDRE_DERIV2_ALT
#include "legendre_source.c"
#undef LEGENDRE_DERIV2_ALT

/*
gsl_sf_legendre_precompute()
  Precompute recurrence factors for ALFs

Inputs: norm         - ALF normalization
        lmax         - maximum degree
        csphase      - 1 to use Condon-Shortley phase factor; 0 to omit
        output_array - (output) output array, length determined by
                       gsl_sf_legendre_array_n()

Notes:
1) Precomputed factors are stored at the back of output_array, skipping
nlm entries which will be used for the ALF values

2) The precomputed factors are:
alm:   length nlm
blm:   length nlm
cl:    length lmax+1
dl:    length lmax+1
sqrts: length 2*lmax+2

Special storage:

alm[idx(0,0) = 0] = P(0,0)
alm[idx(1,0) = 1] = P(1,0)/x
alm[idx(1,1) = 2] = unused
blm[idx(0,0) = 0] = unused
blm[idx(1,0) = 1] = unused
blm[idx(1,1) = 2] = unused
dl[0]             = csfac
*/

int
gsl_sf_legendre_precompute(const gsl_sf_legendre_t norm, const size_t lmax,
                           const int csphase, double output_array[])
{
  if (csphase != 0 && csphase != 1)
    {
      GSL_ERROR ("csphase must be 0 or 1", GSL_EDOM);
    }
  else
    {
      const double csfac = csphase ? -1.0 : 1.0;
      const size_t nlm = gsl_sf_legendre_nlm(lmax);
      double * alm = &output_array[nlm]; /* save nlm entries for the ALFs */
      double * cl = alm + 2 * nlm;
      double * dl = cl + lmax + 1;
      double * sqrts = dl + lmax + 1;
      size_t l, m, k;

      /* compute square root factors */
      for (l = 0; l <= 2 * lmax + 1; ++l)
        sqrts[l] = sqrt((double) l);

      /* store csfac in dl[0] which is not needed in the recurrence relations */
      dl[0] = csfac;

      if (norm == GSL_SF_LEGENDRE_SCHMIDT)
        {
          cl[0] = 1.0;
          alm[0] = 1.0; /* S(0,0) */

          if (lmax == 0)
            return GSL_SUCCESS;

          cl[1] = M_SQRT3;
          dl[1] = csfac;
          alm[1] = 1.0; /* S(1,0)/x */

          /* m = 0 terms */
          k = 2; /* idx(2,0) */
          for (l = 2; l <= lmax; ++l)
            {
              /* a_l^0 */
              alm[2*k] = 2.0 - 1.0 / (double) l;

              /* b_l^0 */
              alm[2*k + 1] = -(1.0 - 1.0 / (double) l);

              ++k;
            }

          /* m > 0 terms */
          k = lmax + 1; /* idx(1,1) */
          for (m = 1; m <= lmax; ++m)
            {
              /* alm and blm are unused for l=m and l=m+1 */
              k += 2;

              for (l = m + 2; l <= lmax; ++l)
                {
                  /* a_l^m */
                  alm[2*k] = ((2.0*l - 1.0) / sqrts[l + m]) / sqrts[l - m];

                  /* b_l^m */
                  alm[2*k + 1] = -(sqrts[l + m - 1] / sqrts[l + m]) *
                                  (sqrts[l - m - 1] / sqrts[l - m]);

                  ++k;
                }
            }

          for (l = 2; l <= lmax; ++l)
            {
              cl[l] = sqrts[2 * l + 1];
              dl[l] = csfac * sqrt(1.0 - 0.5 / l);
            }
        }
      else if (norm == GSL_SF_LEGENDRE_SPHARM)
        {
          alm[0] = 0.5 / M_SQRTPI;    /* Y(0,0) */
          cl[0] = M_SQRT3;

          if (lmax == 0)
            return GSL_SUCCESS;

          alm[1] = sqrt(0.75 / M_PI); /* Y(1,0)/x */
          cl[1] = sqrt(5.0);
          dl[1] = csfac * (M_SQRT3 / M_SQRT2);

          k = 2; /* idx(2,0) */
          for (m = 0; m <= lmax; ++m)
            {
              if (m > 0)
                {
                  /* alm and blm are unused for l=m and l=m+1 */
                  k += 2;
                }

              for (l = m + 2; l <= lmax; ++l)
                {
                  /* a_l^m */
                  alm[2*k] = (sqrts[2 * l + 1] / sqrts[l + m]) *
                             (sqrts[2 * l - 1] / sqrts[l - m]);

                  /* b_l^m */
                  alm[2*k + 1] = -(sqrts[l + m - 1] / sqrts[l + m]) *
                                  (sqrts[l - m - 1] / sqrts[l - m]) *
                                  (sqrts[2*l + 1] / sqrts[2*l - 3]);

                  ++k;
                }
            }

          for (l = 2; l <= lmax; ++l)
            {
              cl[l] = sqrt(2.0 * l + 3.0);
              dl[l] = csfac * sqrt(1.0 + 0.5 / l);
            }
        }
      else if (norm == GSL_SF_LEGENDRE_FULL)
        {
          cl[0] = M_SQRT3;
          alm[0] = M_SQRT1_2;         /* N(0,0) */

          if (lmax == 0)
            return GSL_SUCCESS;

          cl[1] = sqrt(5.0);
          dl[1] = csfac * (M_SQRT3 / M_SQRT2);
          alm[1] = M_SQRT3 / M_SQRT2; /* N(1,0)/x */

          k = 2; /* idx(2,0) */
          for (m = 0; m <= lmax; ++m)
            {
              if (m > 0)
                {
                  /* alm and blm are unused for l=m and l=m+1 */
                  k += 2;
                }

              for (l = m + 2; l <= lmax; ++l)
                {
                  /* a_l^m */
                  alm[2*k] = (sqrts[2 * l + 1] / sqrts[l + m]) *
                             (sqrts[2 * l - 1] / sqrts[l - m]);

                  /* b_l^m */
                  alm[2*k + 1] = -(sqrts[l + m - 1] / sqrts[l + m]) *
                                  (sqrts[l - m - 1] / sqrts[l - m]) *
                                  (sqrts[2*l + 1] / sqrts[2*l - 3]);

                  ++k;
                }
            }

          for (l = 2; l <= lmax; ++l)
            {
              cl[l] = sqrt(2.0 * l + 3.0);
              dl[l] = csfac * sqrt(1.0 + 0.5 / l);
            }
        }
      else if (norm == GSL_SF_LEGENDRE_NONE)
        {
          cl[0] = 1.0;
          alm[0] = 1.0;         /* P(0,0) */

          if (lmax == 0)
            return GSL_SUCCESS;

          cl[1] = 3.0;
          dl[1] = csfac;
          alm[1] = 1.0;         /* P(1,0)/x */

          k = 2; /* idx(2,0) */
          for (m = 0; m <= lmax; ++m)
            {
              if (m > 0)
                {
                  /* alm and blm are unused for l=m and l=m+1 */
                  k += 2;
                }

              for (l = m + 2; l <= lmax; ++l)
                {
                  /* a_l^m */
                  alm[2*k] = (2.0 * l - 1.0) / ((double) (l - m));

                  /* b_l^m */
                  alm[2*k + 1] = -(l + m - 1.0) / ((double) (l - m));

                  ++k;
                }
            }

          for (l = 2; l <= lmax; ++l)
            {
              cl[l] = 2.0 * l + 1.0;
              dl[l] = csfac * (2.0 * l - 1.0);
            }
        }
      else
        {
          GSL_ERROR ("unknown normalization", GSL_EDOM);
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_sf_legendre_array_n()
  Compute total size of array needed for array functions

size = nlm +     // for ALFs
       nlm +     // for alm factors
       nlm +     // for blm factors
       lmax+1    // for cl factors
       lmax+1    // for dl factors
       2*lmax+2  // for sqrt factors

Inputs: lmax - maximum degree
*/

size_t
gsl_sf_legendre_array_n(const size_t lmax)
{
  const size_t nlm = gsl_sf_legendre_nlm(lmax);
  const size_t size = 3 * nlm + 4 * (lmax + 1);
  return size;
}

/* number of P_{lm} functions for a given lmax */
size_t
gsl_sf_legendre_nlm(const size_t lmax)
{
  return ((lmax + 1) * (lmax + 2) / 2);
}

/*
gsl_sf_legendre_arrayx()
  Compute array of associated Legendre functions at a
given point x

Inputs: norm         - ALF normalization
        lmax         - maximum degree
        x            - input point
        result_array - (input/output) output array of ALFs

Notes:
1) result_array must be of length returned by gsl_sf_legendre_array_n()
and be initialized by gsl_sf_legendre_precompute()

2) result_array is indexed by gsl_sf_legendre_array_index(l,m),
result_array[index(l,m)] = Plm(x)
*/

int
gsl_sf_legendre_arrayx(const gsl_sf_legendre_t norm, const size_t lmax,
                       const double x, double result_array[])
{
  if (x < -1.0 || x > 1.0)
    {
      GSL_ERROR ("x is outside [-1,1]", GSL_EDOM);
    }
  else if (x == -1.0 || x == 1.0) /* endpoints */
    {
      const size_t nlm = gsl_sf_legendre_nlm(lmax);
      const double * alm = &result_array[nlm];
      const double * cl = alm + 2 * nlm;
      const double * dl = cl + lmax + 1;
      const double * sqrts = dl + lmax + 1;
      size_t l, m, k;

      /* l=0 m=0 term */
      result_array[0] = alm[0];

      /* check for quick return */
      if (lmax == 0)
        return GSL_SUCCESS;

      /* set P(l,m) = 0 for m > 0 */
      k = 0;
      for (l = 1; l <= lmax; ++l)
        {
          k += l;
          for (m = 1; m <= l; ++m)
            result_array[k + m] = 0.0;
        }

      /* set m=0 terms according to normalization */
      k = 0; /* idx(0,0) */
      if (norm == GSL_SF_LEGENDRE_SCHMIDT)
        {
          if (x == 1.0)
            {
              for (l = 1; l <= lmax; ++l)
                {
                  k += l;
                  result_array[k] = 1.0;
                }
            }
          else
            {
              for (l = 1; l <= lmax; ++l)
                {
                  k += l;
                  if (l & 1)
                    result_array[k] = -1.0;
                  else
                    result_array[k] =  1.0;
                }
            }
        }
      else if (norm == GSL_SF_LEGENDRE_SPHARM)
        {
          if (x == 1.0)
            {
              for (l = 1; l <= lmax; ++l)
                {
                  k += l;
                  result_array[k] = alm[0] * sqrts[2*l + 1];
                }
            }
          else
            {
              for (l = 1; l <= lmax; ++l)
                {
                  k += l;
                  if (l & 1)
                    result_array[k] = -alm[0] * sqrts[2*l + 1];
                  else
                    result_array[k] =  alm[0] * sqrts[2*l + 1];
                }
            }
        }
      else if (norm == GSL_SF_LEGENDRE_FULL)
        {
          if (x == 1.0)
            {
              for (l = 1; l <= lmax; ++l)
                {
                  k += l;
                  result_array[k] = sqrt(l + 0.5);
                }
            }
          else
            {
              for (l = 1; l <= lmax; ++l)
                {
                  k += l;
                  if (l & 1)
                    result_array[k] = -sqrt(l + 0.5);
                  else
                    result_array[k] =  sqrt(l + 0.5);
                }
            }
        }
      else if (norm == GSL_SF_LEGENDRE_NONE)
        {
          if (x == 1.0)
            {
              for (l = 1; l <= lmax; ++l)
                {
                  k += l;
                  result_array[k] = 1.0;
                }
            }
          else
            {
              for (l = 1; l <= lmax; ++l)
                {
                  k += l;
                  if (l & 1)
                    result_array[k] = -1.0;
                  else
                    result_array[k] =  1.0;
                }
            }
        }
      else
        {
          GSL_ERROR ("unknown normalization", GSL_EDOM);
        }

      return GSL_SUCCESS;
    }
  else
    {
      /* interior point */

      const double eps = 1.0e-280;
      const double u = sqrt((1.0 - x) * (1.0 + x)); /* sin(theta) */
      const size_t Lp2 = lmax + 2;
      const size_t nlm = gsl_sf_legendre_nlm(lmax);
      const double * alm = &result_array[nlm];
      const double * cl = alm + 2 * nlm;
      const double * dl = cl + lmax + 1;
      size_t l, m, k;
      double plm,      /* eps * P(l,m) / u^m */
             pmm;      /* eps * P(m,m) / u^m */
      double plmp1;    /* P(l+1,m) */
      double rescalem; /* u^m / eps */
      size_t idxmm;    /* idx(m,m) for output array */
      size_t aidxmm;   /* aidx(m,m) for precomputed arrays */
      const double *al;

      /* initial values P(0,0) and P(1,0) */

      plm = alm[0];       /* P(0,0) */
      result_array[0] = plm;

      /* check for quick return */
      if (lmax == 0)
        return GSL_SUCCESS;

      plmp1 = alm[1] * x; /* P(1,0) */
      result_array[1] = plmp1;

      /* Compute P(l,0), l=2:lmax */

      k = 3;        /* idx(2,0) */
      al = &alm[4]; /* aidx(2,0) = 2 */
      for (l = 2; l < lmax; l += 2)
        {
          plm   = (al[0] * x) * plmp1 + al[1] * plm;
          plmp1 = (al[2] * x) * plm   + al[3] * plmp1;

          result_array[k]         = plm;
          result_array[k + l + 1] = plmp1;

          k += (l << 1) + 3; al += 4;
        }

      if (l == lmax)
        result_array[k] = (al[0] * x) * plmp1 + al[1] * plm;

      /* compute P(m,m), P(m+1,m) and P(l,m) */

      pmm = result_array[0] * eps; /* eps * P(0,0) */
      rescalem = 1.0 / eps;
      idxmm = 0;                   /* idx(0,0) for result_array */
      aidxmm = 0;                  /* idx(0,0) for alm */

      for (m = 1; m < lmax; ++m)
        {
          /* rescalem = u^m / eps */
          rescalem *= u;

          /* compute P(m,m) = d_m * u * P(m-1,m-1) */
          idxmm += m + 1; /* idx(m,m) = idx(m-1,m-1) + m + 1 */
          pmm *= dl[m];
          result_array[idxmm] = pmm * rescalem;
          plm = pmm;

          /* compute P(m+1,m) = c_m * x * P(m,m) */
          plmp1 = (cl[m] * x) * plm;
          result_array[idxmm + m + 1] = plmp1 * rescalem;

          /* compute P(l,m) for l=m+2:lmax */
          aidxmm += Lp2 - m;            /* aidx(m,m) = aidx(m-1,m-1) + L + 2 - m */
          al = &alm[(aidxmm + 2) << 1]; /* aidx(m+2,m) = aidx(m,m) + 2 */
          k = idxmm + 2 * m + 3;        /* idx(m+2,m) = idx(m,m) + 2 * m + 3 */
          for (l = m + 2; l < lmax; l += 2)
            {
              plm   = (al[0] * x) * plmp1 + al[1] * plm;
              plmp1 = (al[2] * x) * plm   + al[3] * plmp1;

              result_array[k]         = plm * rescalem;
              result_array[k + l + 1] = plmp1 * rescalem;

              k += (l << 1) + 3; al += 4;
            }

          if (l == lmax)
            {
              plm = (al[0] * x) * plmp1 + al[1] * plm;
              result_array[k] = plm * rescalem;
            }
        }

      /* compute P(lmax,lmax) */

      rescalem *= u;
      idxmm += m + 1; /* idx(lmax,lmax) */
      pmm *= dl[lmax];
      result_array[idxmm] = pmm * rescalem;

      return GSL_SUCCESS;
    }
}

/*
gsl_sf_legendre_deriv_alt_arrayx()
  Compute array of associated Legendre functions and their
derivatives with respect to theta at a given point x,

d/dtheta Plm(x)

Inputs: norm               - ALF normalization
        lmax               - maximum degree
        x                  - input point, cos(theta)
        result_array       - (input/output) output array of ALFs Plm(x)
        result_deriv_array - (output) d/dtheta Plm(x), length nlm

Notes:
1) result_array must be of length returned by gsl_sf_legendre_array_n()
and be initialized by gsl_sf_legendre_precompute()

2) result_deriv_array must be of length at least nlm

3) output arrays are indexed by gsl_sf_legendre_array_index(l,m)
*/

int
gsl_sf_legendre_deriv_alt_arrayx(const gsl_sf_legendre_t norm,
                                 const size_t lmax,
                                 const double x,
                                 double result_array[],
                                 double result_deriv_array[])
{
  int status;
  const size_t nlm = gsl_sf_legendre_nlm(lmax);
  const double * alm = &result_array[nlm];
  const double * cl = alm + 2 * nlm;
  const double * dl = cl + lmax + 1;
  const double * sqrts = dl + lmax + 1;
  const double csfac = dl[0];

  /* compute Plm(x) */
  status = gsl_sf_legendre_arrayx(norm, lmax, x, result_array);
  if (status)
    return status;

  /* compute d/dtheta Plm(x) */
  status = legendre_derivk_alt_array(norm, lmax, csfac, sqrts,
                                     result_array, result_deriv_array);
  if (status)
    return status;

  return GSL_SUCCESS;
}

/*
gsl_sf_legendre_deriv_arrayx()
  Compute array of associated Legendre functions and their
derivatives with at a given point x,

d/dx Plm(x)

Inputs: norm               - ALF normalization
        lmax               - maximum degree
        x                  - input point, cos(theta)
        result_array       - (input/output) output array of ALFs Plm(x)
        result_deriv_array - (output) d/dx Plm(x), length nlm

Notes:
1) result_array must be of length returned by gsl_sf_legendre_array_n()
and be initialized by gsl_sf_legendre_precompute()

2) result_deriv_array must be of length at least nlm

3) output arrays are indexed by gsl_sf_legendre_array_index(l,m)
*/

int
gsl_sf_legendre_deriv_arrayx(const gsl_sf_legendre_t norm,
                             const size_t lmax,
                             const double x,
                             double result_array[],
                             double result_deriv_array[])
{
  if (x == -1.0 || x == 1.0)
    {
      GSL_ERROR ("x cannot equal -1 or +1 for derivative computation", GSL_EDOM);
    }
  else
    {
      int status;
      const size_t nlm = gsl_sf_legendre_nlm(lmax);
      const double u = sqrt((1.0 - x) * (1.0 + x)); /* sin(theta) */
      const double uinv = 1.0 / u;
      size_t i;

      /* compute Plm and d/dtheta Plm */
      status = gsl_sf_legendre_deriv_alt_arrayx(norm, lmax, x, result_array, result_deriv_array);
      if (status)
        return status;

      /* d/dx Plm = -1/sin(theta) d/dtheta Plm */
      for (i = 0; i < nlm; ++i)
        result_deriv_array[i] *= -uinv;

      return GSL_SUCCESS;
    }
}

/*
gsl_sf_legendre_deriv2_alt_arrayx()
  Compute array of associated Legendre functions and their
first and second derivatives with respect to theta at a given point x

Inputs: norm                - ALF normalization
        lmax                - maximum degree
        x                   - input point, cos(theta)
        result_array        - (input/output) output array of ALFs Plm(x)
        result_deriv_array  - (output) d/dtheta Plm(x), length nlm
        result_deriv2_array - (output) d^2/dtheta^2 Plm(x), length nlm

Notes:
1) result_array must be of length returned by gsl_sf_legendre_array_n()
and be initialized by gsl_sf_legendre_precompute()

2) result_deriv_array and result_deriv2_array must be of length at least nlm

3) output arrays are indexed by gsl_sf_legendre_array_index(l,m)
*/

int
gsl_sf_legendre_deriv2_alt_arrayx(const gsl_sf_legendre_t norm,
                                  const size_t lmax,
                                  const double x,
                                  double result_array[],
                                  double result_deriv_array[],
                                  double result_deriv2_array[])
{
  int status;
  const size_t nlm = gsl_sf_legendre_nlm(lmax);
  const double * alm = &result_array[nlm];
  const double * cl = alm + 2 * nlm;
  const double * dl = cl + lmax + 1;
  const double * sqrts = dl + lmax + 1;
  const double csfac = dl[0];

  /* compute Plm(x) and d/dtheta Plm(x) */
  status = gsl_sf_legendre_deriv_alt_arrayx(norm, lmax, x, result_array, result_deriv_array);
  if (status)
    return status;

  /* compute d^2/dtheta^2 Plm(x) */
  status = legendre_derivk_alt_array(norm, lmax, csfac, sqrts,
                                     result_deriv_array, result_deriv2_array);
  if (status)
    return status;

  return GSL_SUCCESS;
}

/*
gsl_sf_legendre_deriv2_alt_arrayx()
  Compute array of associated Legendre functions and their
first and second derivatives at a given point x

Inputs: norm                - ALF normalization
        lmax                - maximum degree
        x                   - input point, cos(theta)
        result_array        - (input/output) output array of ALFs Plm(x)
        result_deriv_array  - (output) d/dx Plm(x), length nlm
        result_deriv2_array - (output) d^2/dx^2 Plm(x), length nlm

Notes:
1) result_array must be of length returned by gsl_sf_legendre_array_n()
and be initialized by gsl_sf_legendre_precompute()

2) result_deriv_array and result_deriv2_array must be of length at least nlm

3) output arrays are indexed by gsl_sf_legendre_array_index(l,m)
*/

int
gsl_sf_legendre_deriv2_arrayx(const gsl_sf_legendre_t norm,
                              const size_t lmax,
                              const double x,
                              double result_array[],
                              double result_deriv_array[],
                              double result_deriv2_array[])
{
  if (x == -1.0 || x == 1.0)
    {
      GSL_ERROR ("x cannot equal -1 or +1 for derivative computation", GSL_EDOM);
    }
  else
    {
      int status;
      const size_t nlm = gsl_sf_legendre_nlm(lmax);
      const double u = sqrt((1.0 - x) * (1.0 + x)); /* sin(theta) */
      const double uinv = 1.0 / u;
      const double uinv2 = uinv * uinv;
      size_t i;

      /* compute Plm, d/dtheta Plm and d^2/dtheta^2 Plm */
      status = gsl_sf_legendre_deriv2_alt_arrayx(norm, lmax, x,
                                                 result_array,
                                                 result_deriv_array,
                                                 result_deriv2_array);
      if (status)
        return status;

      for (i = 0; i < nlm; ++i)
        {
          /* d/dx Plm = -1/sin(theta) d/dtheta Plm */
          result_deriv_array[i] *= -uinv;

          /* d^2/dx^2 Plm = 1/sin^2(theta) [ d^2/dtheta^2 Plm + x d/dx Plm ] */
          result_deriv2_array[i] = uinv2 * (result_deriv2_array[i] + x * result_deriv_array[i]);
        }

      return GSL_SUCCESS;
    }
}

int
gsl_sf_legendre_array(const gsl_sf_legendre_t norm,
                      const size_t lmax,
                      const double x,
                      double result_array[])
{
  int status;

  status = gsl_sf_legendre_precompute(norm, lmax, 0, result_array);
  if (status)
    return status;

  status = gsl_sf_legendre_arrayx(norm, lmax, x, result_array);
  if (status)
    return status;

  return GSL_SUCCESS;
}

int
gsl_sf_legendre_deriv_alt_array(const gsl_sf_legendre_t norm, const size_t lmax,
                                const double x, double result_array[],
                                double result_deriv_array[])
{
  int status;

  status = gsl_sf_legendre_precompute(norm, lmax, 0, result_array);
  if (status)
    return status;

  status = gsl_sf_legendre_deriv_alt_arrayx(norm, lmax, x, result_array, result_deriv_array);
  if (status)
    return status;

  return GSL_SUCCESS;
}

int
gsl_sf_legendre_deriv_array(const gsl_sf_legendre_t norm, const size_t lmax,
                            const double x, double result_array[],
                            double result_deriv_array[])
{
  int status;

  status = gsl_sf_legendre_precompute(norm, lmax, 0, result_array);
  if (status)
    return status;

  status = gsl_sf_legendre_deriv_arrayx(norm, lmax, x, result_array, result_deriv_array);
  if (status)
    return status;

  return GSL_SUCCESS;
}

int
gsl_sf_legendre_deriv2_alt_array(const gsl_sf_legendre_t norm, const size_t lmax,
                                 const double x, double result_array[],
                                 double result_deriv_array[], double result_deriv2_array[])
{
  int status;

  status = gsl_sf_legendre_precompute(norm, lmax, 0, result_array);
  if (status)
    return status;

  status = gsl_sf_legendre_deriv2_alt_arrayx(norm, lmax, x, result_array,
                                             result_deriv_array, result_deriv2_array);
  if (status)
    return status;

  return GSL_SUCCESS;
}

int
gsl_sf_legendre_deriv2_array(const gsl_sf_legendre_t norm, const size_t lmax,
                             const double x, double result_array[],
                             double result_deriv_array[], double result_deriv2_array[])
{
  int status;

  status = gsl_sf_legendre_precompute(norm, lmax, 0, result_array);
  if (status)
    return status;

  status = gsl_sf_legendre_deriv2_arrayx(norm, lmax, x, result_array,
                                         result_deriv_array, result_deriv2_array);
  if (status)
    return status;

  return GSL_SUCCESS;
}

/*********************************************************
 *                 INTERNAL ROUTINES                     *
 *********************************************************/

/*
legendre_derivk_alt_array()
  Compute the kth derivative of Plm(x) with respect to theta
using Eq. 14 of [2]

Inputs: norm                - ALF normalization
        lmax                - maximum degree
        csfac               - -1 to include CS phase, +1 to omit
        sqrts               - sqrts[i] = sqrt(i), length 2*lmax+2
        input_array         - (input) array of d^{k-1}/dtheta^{k-1} Plm(x), length nlm
        output_array        - (output) d^k/dtheta^k Plm(x), length nlm
*/

static int
legendre_derivk_alt_array(const gsl_sf_legendre_t norm,
                          const size_t lmax,
                          const double csfac,
                          const double sqrts[],
                          const double input_array[],
                          double output_array[])
{
  size_t l, m, k;

  /* d^k/dtheta^k P(0,0) = 0 */
  output_array[0] = 0.0;

  if (lmax == 0)
    return GSL_SUCCESS;

  if (norm == GSL_SF_LEGENDRE_NONE)
    {
      k = 0; /* idx(0,0) */
      for (l = 1; l <= lmax; ++l)
        {
          k += l; /* idx(l,0) */

          /* Eq. 13: compute d^k/dtheta^k P(l,0) = -d^{k-1}/dtheta^{k-1} P(l,1) */
          output_array[k] = -csfac * input_array[k + 1];

          /* Eq. 12: compute d^k/dtheta^k P(l,l) = l d^{k-1}/dtheta^{k-1} P(l,l-1) */
          output_array[k + l] = csfac * l * input_array[k + l - 1];

          for (m = 1; m < l; ++m)
            {
              output_array[k + m] = 0.5 * csfac * ((l + m) * (l - m + 1.0) * input_array[k + m - 1] -
                                                   input_array[k + m + 1]);
            }
        }
    }
  else if (norm == GSL_SF_LEGENDRE_SCHMIDT)
    {
      /* d^k/dtheta^k P(1,0) = -d^{k-1}/dtheta^{k-1} P(1,1) */
      output_array[1] = -csfac * input_array[2];

      /* d^k/dtheta^k P(1,1) = d^{k-1}/dtheta^{k-1} P(1,0) */
      output_array[2] = csfac * input_array[1];

      k = 1; /* idx(1,0) */
      for (l = 2; l <= lmax; ++l)
        {
          k += l; /* idx(l,0) */

          /* d^k/dtheta^k P(l,0) = -sqrt(l(l+1)/2) d^{k-1}/dtheta^{k-1} P(l,1) */
          output_array[k] = -csfac * (sqrts[l] / M_SQRT2) * sqrts[l + 1] * input_array[k + 1];

          /* d^k/dtheta^k P(l,1) = sqrt(l(l+1)/2) d^{k-1}/dtheta^{k-1} P(l,0) -
           *                       1/2 sqrt((l-1)(l+2)) d^{k-1}/dtheta^{k-1} P(l,2) */
          output_array[k + 1] = csfac * ((sqrts[l] / M_SQRT2) * sqrts[l + 1] * input_array[k] -
                                         0.5 * sqrts[l - 1] * sqrts[l + 2] * input_array[k + 2]);

          /* d^k/dtheta^k P(l,l) = sqrt(l/2) d^{k-1}/dtheta^{k-1} P(l,l-1) */
          output_array[k + l] = csfac * (sqrts[l] / M_SQRT2) * input_array[k + l - 1];

          for (m = 2; m < l; ++m)
            {
              output_array[k + m] = 0.5 * csfac * (sqrts[l + m] * sqrts[l - m + 1] * input_array[k + m - 1] -
                                                   sqrts[l + m + 1] * sqrts[l - m] * input_array[k + m + 1]);
            }
        }
    }
  else
    {
      k = 0; /* idx(0,0) */
      for (l = 1; l <= lmax; ++l)
        {
          k += l; /* idx(l,0) */

          /* d^k/dtheta^k P(l,0) = -sqrt(l(l+1)) d^{k-1}/dtheta^{k-1} P(l,1) */
          output_array[k] = -csfac * sqrts[l] * sqrts[l + 1] * input_array[k + 1];

          /* d^k/dtheta^k P(l,l) = sqrt(l/2) d^{k-1}/dtheta^{k-1} P(l,l-1) */
          output_array[k + l] = csfac * (sqrts[l] / M_SQRT2) * input_array[k + l - 1];

          for (m = 1; m < l; ++m)
            {
              output_array[k + m] = 0.5 * csfac * (sqrts[l + m] * sqrts[l - m + 1] * input_array[k + m - 1] -
                                                   sqrts[l + m + 1] * sqrts[l - m] * input_array[k + m + 1]);
            }
        }
    }

  return GSL_SUCCESS;
}

/*
legendre_sqrts()
  Precompute square root factors needed for Legendre recurrence.
On output, array[i] = sqrt(i)
*/

static void
legendre_sqrts(const size_t lmax, double *array)
{
  size_t l;
  for (l = 0; l <= 2 * lmax + 1; ++l)
    array[l] = sqrt((double) l);
}
