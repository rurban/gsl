/* linalg/trimult_complex.c
 *
 * Copyright (C) 2019 Patrick Alken
 *
 * This is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation; either version 3, or (at your option) any
 * later version.
 *
 * This source is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 * for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 *
 * This module contains code to compute L^T L where L is a lower triangular matrix
 */

#include <config.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>

#include "recurse.h"

static int triangular_multiply_L2(CBLAS_UPLO_t Uplo, gsl_matrix_complex * T);
static int triangular_multiply_L3(CBLAS_UPLO_t Uplo, gsl_matrix_complex * T);
static void complex_conj_vector(gsl_vector_complex * v);

int
gsl_linalg_complex_tri_LHL(gsl_matrix_complex * L)
{
  return triangular_multiply_L3(CblasLower, L);
}

/*
triangular_multiply_L2()
  Compute L^H L or U U^H

Inputs: Uplo - CblasUpper or CblasLower
        T    - on output the upper (or lower) part of T
               is replaced by L^H L or U U^H

Return: success/error

Notes:
1) Based on LAPACK routine ZLAUU2 using Level 2 BLAS
*/

static int
triangular_multiply_L2(CBLAS_UPLO_t Uplo, gsl_matrix_complex * T)
{
  const size_t N = T->size1;

  if (N != T->size2)
    {
      GSL_ERROR ("matrix must be square", GSL_ENOTSQR);
    }
  else
    {
      size_t i;

      if (Uplo == CblasUpper)
        {
        }
      else
        {
          for (i = 0; i < N; ++i)
            {
              gsl_complex * Tii = gsl_matrix_complex_ptr(T, i, i);
              gsl_complex z0 = *Tii;

              if (i < N - 1)
                {
                  gsl_vector_complex_view v = gsl_matrix_complex_subcolumn(T, i, i + 1, N - i - 1);
                  double norm = gsl_blas_dznrm2(&v.vector);

                  GSL_REAL(*Tii) = gsl_complex_abs2(*Tii) + norm * norm;

                  if (i > 0)
                    {
                      gsl_vector_complex_view w = gsl_matrix_complex_subrow(T, i, 0, i);
                      gsl_matrix_complex_view m = gsl_matrix_complex_submatrix(T, i + 1, 0, N - i - 1, i);

                      complex_conj_vector(&w.vector);
                      gsl_blas_zgemv(CblasConjTrans, GSL_COMPLEX_ONE, &m.matrix, &v.vector, z0, &w.vector);
                      complex_conj_vector(&w.vector);
                    }
                }
              else
                {
                  gsl_vector_complex_view w = gsl_matrix_complex_row(T, i);
                  gsl_blas_zdscal(GSL_REAL(z0), &w.vector);
                }

              GSL_IMAG(*Tii) = 0.0;
            }
        }

      return GSL_SUCCESS;
    }
}

/*
triangular_multiply_L3()
  Compute L^H L or U U^H

Inputs: Uplo - CblasUpper or CblasLower
        T    - on output the upper (or lower) part of T
               is replaced by L^H L or U U^H

Return: success/error

Notes:
1) Based on ReLAPACK using Level 3 BLAS
*/

static int
triangular_multiply_L3(CBLAS_UPLO_t Uplo, gsl_matrix_complex * T)
{
  const size_t N = T->size1;

  if (N != T->size2)
    {
      GSL_ERROR ("matrix must be square", GSL_ENOTSQR);
    }
  else if (N <= CROSSOVER_TRIMULT)
    {
      /* use Level 2 algorithm */
      return triangular_multiply_L2(Uplo, T);
    }
  else
    {
      /*
       * partition matrix:
       *
       * T11 T12
       * T21 T22
       *
       * where T11 is N1-by-N1
       */
      int status;
      const size_t N1 = GSL_LINALG_SPLIT_COMPLEX(N);
      const size_t N2 = N - N1;
      gsl_matrix_complex_view T11 = gsl_matrix_complex_submatrix(T, 0, 0, N1, N1);
      gsl_matrix_complex_view T12 = gsl_matrix_complex_submatrix(T, 0, N1, N1, N2);
      gsl_matrix_complex_view T21 = gsl_matrix_complex_submatrix(T, N1, 0, N2, N1);
      gsl_matrix_complex_view T22 = gsl_matrix_complex_submatrix(T, N1, N1, N2, N2);

      /* recursion on T11 */
      status = triangular_multiply_L3(Uplo, &T11.matrix);
      if (status)
        return status;

      if (Uplo == CblasLower)
        {
          /* T11 += T21^T T21 */
          gsl_blas_zherk(Uplo, CblasConjTrans, 1.0, &T21.matrix, 1.0, &T11.matrix);

          /* T21 = T22^T * T21 */
          gsl_blas_ztrmm(CblasLeft, Uplo, CblasConjTrans, CblasNonUnit, GSL_COMPLEX_ONE, &T22.matrix, &T21.matrix);
        }
      else
        {
          /* T11 += T12 T12^T */
          gsl_blas_zherk(Uplo, CblasNoTrans, 1.0, &T12.matrix, 1.0, &T11.matrix);

          /* T12 = T12 * T22^T */
          gsl_blas_ztrmm(CblasRight, Uplo, CblasConjTrans, CblasNonUnit, GSL_COMPLEX_ONE, &T22.matrix, &T12.matrix);
        }

      /* recursion on T22 */
      status = triangular_multiply_L3(Uplo, &T22.matrix);
      if (status)
        return status;

      return GSL_SUCCESS;
    }
}

/********************************************
 *           INTERNAL ROUTINES              *
 ********************************************/

static void
complex_conj_vector(gsl_vector_complex * v)
{
  size_t i;

  for (i = 0; i < v->size; ++i)
    {
      gsl_complex * vi = gsl_vector_complex_ptr(v, i);
      GSL_IMAG(*vi) = -GSL_IMAG(*vi);
    }
}
