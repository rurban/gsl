/* linalg/rqr.c
 * 
 * Copyright (C) 2019 Patrick Alken
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
#include <string.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

/*
 * this module contains routines for the QR factorization of a matrix
 * using the recursive Level 3 BLAS algorithm of Elmroth and Gustavson with
 * additional modifications courtesy of Julien Langou.
 */

static int unpack_Q1_r(gsl_matrix * Q);
static int aux_ULT(const gsl_matrix * L, gsl_matrix * U);
static int aux_mLU(gsl_matrix * A);
static int aux_ApUBT(const gsl_matrix * U, const gsl_matrix * B, gsl_matrix * A);

/*
gsl_linalg_QR_decomp_r()
  QR decomposition using Level 3 BLAS recursive algorithm of:
  
Elmroth, E. and Gustavson, F.G., 2000. Applying recursion to serial and parallel
  QR factorization leads to better performance. IBM Journal of Research and Development,
  44(4), pp.605-624.

Inputs: A - matrix to be factored, M-by-N with M >= N
        T - N-by-N upper triangular factor of block reflector

Return: success/error

Notes:
1) on output, diag(T) contains tau vector

2) on output, upper triangle of A contains R; elements below the diagonal
are columns of V, where the block reflector H is:

H = I - V T V^T

3) implementation provided by Julien Langou
*/

int
gsl_linalg_QR_decomp_r (gsl_matrix * A, gsl_matrix * T)
{
  const size_t M = A->size1;
  const size_t N = A->size2;

  if (M < N)
    {
      GSL_ERROR ("M must be >= N", GSL_EBADLEN);
    }
  else if (T->size1 != T->size2)
    {
      GSL_ERROR ("T matrix must be square", GSL_ENOTSQR);
    }
  else if (T->size1 != N)
    {
      GSL_ERROR ("T matrix does not match dimensions of A", GSL_EBADLEN);
    }
  else
    {
      if (N == 1)
        {
          /* base case, compute householder transform for single column matrix */

          double * T00 = gsl_matrix_ptr(T, 0, 0);
          gsl_vector_view v = gsl_matrix_column(A, 0);

          *T00 = gsl_linalg_householder_transform(&v.vector);
        }
      else
        {
          /*
           * partition matrices:
           *
           *       N1  N2              N1  N2
           * N1 [ A11 A12 ] and  N1 [ T11 T12 ]
           * M2 [ A21 A22 ]      N2 [  0  T22 ]
           */
          int status;
          const size_t N1 = N / 2;
          const size_t N2 = N - N1;
          const size_t M2 = M - N1;

          gsl_matrix_view A11 = gsl_matrix_submatrix(A, 0, 0, N1, N1);
          gsl_matrix_view A12 = gsl_matrix_submatrix(A, 0, N1, N1, N2);
          gsl_matrix_view A21 = gsl_matrix_submatrix(A, N1, 0, M2, N1);
          gsl_matrix_view A22 = gsl_matrix_submatrix(A, N1, N1, M2, N2);

          gsl_matrix_view T11 = gsl_matrix_submatrix(T, 0, 0, N1, N1);
          gsl_matrix_view T12 = gsl_matrix_submatrix(T, 0, N1, N1, N2);
          gsl_matrix_view T22 = gsl_matrix_submatrix(T, N1, N1, N2, N2);

          gsl_matrix_view m;

          /* recursion on (A(1:m,1:N1), T11) */
          m = gsl_matrix_submatrix(A, 0, 0, M, N1);
          status = gsl_linalg_QR_decomp_r(&m.matrix, &T11.matrix);
          if (status)
            return status;

          gsl_matrix_memcpy(&T12.matrix, &A12.matrix);

          gsl_blas_dtrmm(CblasLeft, CblasLower, CblasTrans, CblasUnit, 1.0, &A11.matrix, &T12.matrix);       /* T12 = lower(A11)' * T12 */
          gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &A21.matrix, &A22.matrix, 1.0, &T12.matrix);         /* T12 = T12 + A21' * A22 */
          gsl_blas_dtrmm(CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, 1.0, &T11.matrix, &T12.matrix);    /* T12 = T11' * T12 */
          gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, &A21.matrix, &T12.matrix, 1.0, &A22.matrix);      /* A22 = A22 - A21 * T12 */
          gsl_blas_dtrmm(CblasLeft, CblasLower, CblasNoTrans, CblasUnit, 1.0, &A11.matrix, &T12.matrix);     /* T12 = lower(A11) * T12 */

          gsl_matrix_sub(&A12.matrix, &T12.matrix);

          /* recursion on (A22, T22) */
          status = gsl_linalg_QR_decomp_r(&A22.matrix, &T22.matrix);
          if (status)
            return status;

          m = gsl_matrix_submatrix(&A21.matrix, 0, 0, N2, N1);
          gsl_matrix_transpose_memcpy(&T12.matrix, &m.matrix);

          A22 = gsl_matrix_submatrix(A, N1, N1, N2, N2);
          gsl_blas_dtrmm(CblasRight, CblasLower, CblasNoTrans, CblasUnit, 1.0, &A22.matrix, &T12.matrix);    /* T12 = T12 * lower(A22) */

          if (M > N)
            {
              gsl_matrix_view A31 = gsl_matrix_submatrix(A, N, 0, M - N, N1);
              gsl_matrix_view A32 = gsl_matrix_submatrix(A, N, N1, M - N, N2);

              gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, &A31.matrix, &A32.matrix, 1.0, &T12.matrix);     /* T12 = T12 + A31' * A32 */
            }

          gsl_blas_dtrmm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, -1.0, &T11.matrix, &T12.matrix); /* T12 = -T11 * T12 */
          gsl_blas_dtrmm(CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, &T22.matrix, &T12.matrix); /* T12 = T12 * T22 */
        }

      return GSL_SUCCESS;
    }
}

/*
gsl_linalg_QR_unpack_r()
  Unpack matrices Q and R

Inputs: QR - packed QR format, M-by-N
        T  - block reflector matrix, N-by-N
        Q  - (output) Q matrix, M-by-M
        R  - (output) R matrix, N-by-N

Return: success/error

Notes:
1) Implementation provided by Julien Langou

2) Lower triangular portion of R is used as temporary workspace
*/

int
gsl_linalg_QR_unpack_r(const gsl_matrix * QR, const gsl_matrix * T, gsl_matrix * Q, gsl_matrix * R)
{
  const size_t M = QR->size1;
  const size_t N = QR->size2;

  if (M < N)
    {
      GSL_ERROR ("M must be >= N", GSL_EBADLEN);
    }
  else if (Q->size1 != M || Q->size2 != M)
    {
      GSL_ERROR ("Q matrix must be M-by-M", GSL_EBADLEN);
    }
  else if (R->size1 != N || R->size2 != N)
    {
      GSL_ERROR ("R matrix must be N-by-N", GSL_EBADLEN);
    }
  else if (T->size1 != N || T->size2 != N)
    {
      GSL_ERROR ("T matrix must be N-by-N", GSL_EBADLEN);
    }
  else
    {
      gsl_matrix_const_view RV = gsl_matrix_const_submatrix(QR, 0, 0, N, N);
      gsl_matrix_view Q1 = gsl_matrix_submatrix(Q, 0, 0, M, N);
      gsl_matrix_view m;

      /*
       * set Q1 = [ T ]
       *          [ V ]
       */
      m = gsl_matrix_submatrix(Q, 0, 0, N, N);
      gsl_matrix_tricpy('U', 1, &m.matrix, T);
      gsl_matrix_tricpy('L', 0, &m.matrix, &RV.matrix);

      if (M > N)
        {
          gsl_matrix_const_view tmp = gsl_matrix_const_submatrix(QR, N, 0, M - N, N);
          m = gsl_matrix_submatrix(Q, N, 0, M - N, N);
          gsl_matrix_memcpy(&m.matrix, &tmp.matrix);
        }

      unpack_Q1_r(&Q1.matrix);

      /* copy R */
      gsl_matrix_tricpy('U', 1, R, &RV.matrix);

      return GSL_SUCCESS;
    }
}

/*
unpack_Q1_r()
  Compute Q_1

Inputs: Q  - on input, contains T in upper triangle and V in lower trapezoid
             on output, contains Q_1
             M-by-N
*/

static int
unpack_Q1_r(gsl_matrix * Q)
{
  int status;
  const size_t M = Q->size1;
  const size_t N = Q->size2;
  gsl_matrix_view T = gsl_matrix_submatrix(Q, 0, 0, N, N);
  size_t i;

  /* T := T V1^T */
  status = aux_ULT(&T.matrix, &T.matrix);
  if (status)
    return status;

  if (M > N)
    {
      gsl_matrix_view m = gsl_matrix_submatrix(Q, N, 0, M - N, N);
      gsl_blas_dtrmm(CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, -1.0, &T.matrix, &m.matrix);
    }

  status = aux_mLU(&T.matrix);
  if (status)
    return status;

  for (i = 0; i < N; ++i)
    {
      double * ptr = gsl_matrix_ptr(Q, i, i);
      *ptr += 1.0;
    }

  return GSL_SUCCESS;
}

/* U := U L^T for triangular matrices L and U; L is unit lower triangular */
static int
aux_ULT(const gsl_matrix * L, gsl_matrix * U)
{
  const size_t N = L->size1;

  if (N != L->size2)
    {
      GSL_ERROR ("L matrix must be square", GSL_ENOTSQR);
    }
  else if (U->size1 != N || U->size2 != N)
    {
      GSL_ERROR ("U matrix must be same size as L", GSL_EBADLEN);
    }
  else if (N == 1)
    {
      /* nothing to do */
      return GSL_SUCCESS;
    }
  else
    {
      int status;
      const size_t N1 = N / 2;
      const size_t N2 = N - N1;

      gsl_matrix_const_view L11 = gsl_matrix_const_submatrix(L, 0, 0, N1, N1);
      gsl_matrix_const_view L21 = gsl_matrix_const_submatrix(L, N1, 0, N2, N1);
      gsl_matrix_const_view L22 = gsl_matrix_const_submatrix(L, N1, N1, N2, N2);

      gsl_matrix_view U11 = gsl_matrix_submatrix(U, 0, 0, N1, N1);
      gsl_matrix_view U12 = gsl_matrix_submatrix(U, 0, N1, N1, N2);
      gsl_matrix_view U22 = gsl_matrix_submatrix(U, N1, N1, N2, N2);

      /* U12 = U12 * L22^T */
      gsl_blas_dtrmm(CblasRight, CblasLower, CblasTrans, CblasUnit, 1.0, &L22.matrix, &U12.matrix);

      /* U12 = U12 + U11 * L21^T */
      status = aux_ApUBT(&U11.matrix, &L21.matrix, &U12.matrix);
      if (status)
        return status;

      status = aux_ULT(&L11.matrix, &U11.matrix);
      if (status)
        return status;

      status = aux_ULT(&L22.matrix, &U22.matrix);
      if (status)
        return status;

      return GSL_SUCCESS;
    }
}

/* store -L*U in A */
static int
aux_mLU(gsl_matrix * A)
{
  const size_t N = A->size1;

  if (N != A->size2)
    {
      GSL_ERROR ("matrix must be square", GSL_ENOTSQR);
    }
  else if (N == 1)
    {
      double *A00 = gsl_matrix_ptr(A, 0, 0);
      *A00 = -(*A00);
      return GSL_SUCCESS;
    }
  else
    {
      int status;
      const size_t N1 = N / 2;
      const size_t N2 = N - N1;

      gsl_matrix_view A11 = gsl_matrix_submatrix(A, 0, 0, N1, N1);
      gsl_matrix_view A12 = gsl_matrix_submatrix(A, 0, N1, N1, N2);
      gsl_matrix_view A21 = gsl_matrix_submatrix(A, N1, 0, N2, N1);
      gsl_matrix_view A22 = gsl_matrix_submatrix(A, N1, N1, N2, N2);

      /* A22 = - L22 U22 */
      status = aux_mLU(&A22.matrix);
      if (status)
        return status;

      /* A22 = A22 - L21 U12 */
      gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, &A21.matrix, &A12.matrix, 1.0, &A22.matrix);

      /* A12 - -L11 U12 */
      gsl_blas_dtrmm(CblasLeft, CblasLower, CblasNoTrans, CblasUnit, -1.0, &A11.matrix, &A12.matrix);

      /* A21 = -L21 U11 */
      gsl_blas_dtrmm(CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, -1.0, &A11.matrix, &A21.matrix);

      /* A11 = - L11 U11 */
      status = aux_mLU(&A11.matrix);
      if (status)
        return status;

      return GSL_SUCCESS;
    }
}

/* A := A + U B^T where U is upper triangular */
static int
aux_ApUBT(const gsl_matrix * U, const gsl_matrix * B, gsl_matrix * A)
{
  const size_t M = A->size1;
  const size_t N = A->size2;

  if (U->size1 != M || U->size2 != M)
    {
      GSL_ERROR ("U matrix has wrong dimensions", GSL_EBADLEN);
    }
  else if (B->size1 != N || B->size2 != M)
    {
      GSL_ERROR ("B matrix has wrong dimensions", GSL_EBADLEN);
    }
  else if (M == 1 && N == 1)
    {
      double *aptr = gsl_matrix_ptr(A, 0, 0);
      const double *uptr = gsl_matrix_const_ptr(U, 0, 0);
      const double *bptr = gsl_matrix_const_ptr(B, 0, 0);
      *aptr += (*uptr) * (*bptr);
      return GSL_SUCCESS;
    }
  else if (M == 1)
    {
      double U00 = gsl_matrix_get(U, 0, 0);
      gsl_vector_view v = gsl_matrix_row(A, 0);
      gsl_vector_const_view w = gsl_matrix_const_column(B, 0);
      gsl_blas_daxpy(U00, &w.vector, &v.vector);
      return GSL_SUCCESS;
    }
  else if (N == 1)
    {
      size_t i;

      for (i = 0; i < M; ++i)
        {
          double * ai = gsl_matrix_ptr(A, i, 0);
          gsl_vector_const_view u = gsl_matrix_const_subrow(U, i, i, M - i);
          gsl_vector_const_view b = gsl_matrix_const_subrow(B, 0, i, M - i);
          double dot;

          gsl_blas_ddot(&u.vector, &b.vector, &dot);
          *ai += dot;
        }

      return GSL_SUCCESS;
    }
  else
    {
      int status;
      const size_t M1 = M / 2;
      const size_t M2 = M - M1;
      const size_t N1 = N / 2;
      const size_t N2 = N - N1;

      gsl_matrix_view A11 = gsl_matrix_submatrix(A, 0, 0, M1, N1);
      gsl_matrix_view A12 = gsl_matrix_submatrix(A, 0, N1, M1, N2);
      gsl_matrix_view A21 = gsl_matrix_submatrix(A, M1, 0, M2, N1);
      gsl_matrix_view A22 = gsl_matrix_submatrix(A, M1, N1, M2, N2);

      gsl_matrix_const_view U11 = gsl_matrix_const_submatrix(U, 0, 0, M1, M1);
      gsl_matrix_const_view U12 = gsl_matrix_const_submatrix(U, 0, M1, M1, M2);
      gsl_matrix_const_view U22 = gsl_matrix_const_submatrix(U, M1, M1, M2, M2);

      gsl_matrix_const_view B11 = gsl_matrix_const_submatrix(B, 0, 0, N1, M1);
      gsl_matrix_const_view B12 = gsl_matrix_const_submatrix(B, 0, M1, N1, M2);
      gsl_matrix_const_view B21 = gsl_matrix_const_submatrix(B, N1, 0, N2, M1);
      gsl_matrix_const_view B22 = gsl_matrix_const_submatrix(B, N1, M1, N2, M2);

      /* A11 := A11 + U11 B11^T */
      status = aux_ApUBT(&U11.matrix, &B11.matrix, &A11.matrix);
      if (status)
        return status;

      /* A11 := A11 + U12 B12^T */
      gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &U12.matrix, &B12.matrix, 1.0, &A11.matrix);

      /* A12 := A12 + U11 B21^T */
      status = aux_ApUBT(&U11.matrix, &B21.matrix, &A12.matrix);
      if (status)
        return status;

      /* A12 := A12 + U12 B22^T */
      gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, &U12.matrix, &B22.matrix, 1.0, &A12.matrix);

      /* A21 := A21 + U22 B12^T */
      status = aux_ApUBT(&U22.matrix, &B12.matrix, &A21.matrix);
      if (status)
        return status;

      /* A22 := A22 + U22 B22^T */
      status = aux_ApUBT(&U22.matrix, &B22.matrix, &A22.matrix);
      if (status)
        return status;

      return GSL_SUCCESS;
    }
}
