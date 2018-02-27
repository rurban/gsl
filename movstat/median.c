/* movstat/median.c
 *
 * Routines related to a moving window median
 * 
 * Copyright (C) 2018 Patrick Alken
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
#include <gsl/gsl_movstat.h>

/*
gsl_movstat_median_alloc()
  Allocate a workspace for median filtering. The workspace
is set up to calculate a running median with a given window size.
The window around sample x_i is defined as:

W_i^{H,J} = {x_{i-H},...,x_i,...x_{i+J}}

The total window size is:

K = H + J + 1

Inputs: K - total samples in window (H = J = K / 2)

Return: pointer to workspace

Notes:
1) If K is even, it is rounded up to the next odd
*/

gsl_movstat_median_workspace *
gsl_movstat_median_alloc(const size_t K)
{
  const size_t H = K / 2;
  return gsl_movstat_median_alloc2(H, H);
}
 
/*
gsl_movstat_median_alloc2()
  Allocate a workspace for median filtering. The workspace
is set up to calculate a running median with a given window size.
The window around sample x_i is defined as:

W_i^{H,J} = {x_{i-H},...,x_i,...x_{i+J}}

The total window size is:

K = H + J + 1

Inputs: H - number of samples before current sample
        J - number of samples after current sample

Return: pointer to workspace
*/

gsl_movstat_median_workspace *
gsl_movstat_median_alloc2(const size_t H, const size_t J)
{
  gsl_movstat_median_workspace *w;

  w = calloc(1, sizeof(gsl_movstat_median_workspace));
  if (w == 0)
    {
      GSL_ERROR_NULL ("failed to allocate space for workspace", GSL_ENOMEM);
    }

  w->H = H;
  w->J = J;
  w->K = H + J + 1;

  w->medacc_workspace_p = gsl_movstat_medacc_alloc(w->K);
  if (w->medacc_workspace_p == 0)
    {
      gsl_movstat_median_free(w);
      GSL_ERROR_NULL ("failed to allocate space for mediator workspace", GSL_ENOMEM);
    }

  return w;
}

void
gsl_movstat_median_free(gsl_movstat_median_workspace * w)
{
  if (w->medacc_workspace_p)
    gsl_movstat_medacc_free(w->medacc_workspace_p);

  free(w);
}

/*
gsl_movstat_median()
  Apply median filter to input vector

Inputs: etype - edge handling criteria
        x     - input vector, size n
        y     - output vector, size n
        w     - workspace
*/

int
gsl_movstat_median(const gsl_movstat_edge_t etype, const gsl_vector * x, gsl_vector * y, gsl_movstat_median_workspace * w)
{
  const size_t n = x->size;

  if (n != y->size)
    {
      GSL_ERROR("input and output vectors must have same length", GSL_EBADLEN);
    }
  else
    {
      const int J = w->J; /* number of samples to right of current sample */
      size_t i;
      double x1 = 0.0;    /* pad values for data edges */
      double xN = 0.0;

      if (etype != GSL_MOVSTAT_EDGE_TRUNCATE)
        {
          if (etype == GSL_MOVSTAT_EDGE_PADZERO)
            {
              x1 = 0.0;
              xN = 0.0;
            }
          else if (etype == GSL_MOVSTAT_EDGE_PADVALUE)
            {
              x1 = gsl_vector_get(x, 0);
              xN = gsl_vector_get(x, n - 1);
            }

          /* pad initial windows with H values */
          for (i = 0; i < w->H; ++i)
            gsl_movstat_medacc_insert(x1, w->medacc_workspace_p);
        }

      /* process input vector and fill y(1:n - J) */
      for (i = 0; i < n; ++i)
        {
          double xi = gsl_vector_get(x, i);
          int idx = (int) i - J;

          gsl_movstat_medacc_insert(xi, w->medacc_workspace_p);

          if (idx >= 0)
            gsl_vector_set(y, idx, gsl_movstat_medacc_median(w->medacc_workspace_p));
        }

      if (etype == GSL_MOVSTAT_EDGE_TRUNCATE)
        {
          for (i = 0; i < w->J; ++i)
            {
            }
        }
      else
        {
          /* pad final windows and fill y(n - J + 1:n) */
          for (i = 0; i < w->J; ++i)
            {
              int idx = (int) n - J + (int) i;

              gsl_movstat_medacc_insert(xN, w->medacc_workspace_p);

              if (idx >= 0)
                gsl_vector_set(y, idx, gsl_movstat_medacc_median(w->medacc_workspace_p));
            }
        }

      return GSL_SUCCESS;
    }
}
