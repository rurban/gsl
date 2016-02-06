/* test.c
 * 
 * Copyright (C) 2012-2014, 2016 Patrick Alken
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
#include <math.h>
#include <unistd.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_spmatrix.h>

/*
create_random_sparse()
  Create a random sparse matrix with approximately
M*N*density non-zero entries

Inputs: M       - number of rows
        N       - number of columns
        density - sparse density \in [0,1]
                  0 = no non-zero entries
                  1 = all m*n entries are filled
        r       - random number generator

Return: pointer to sparse matrix in triplet format (must be freed by caller)

Notes:
1) non-zero matrix entries are uniformly distributed in [0,1]
*/

static gsl_spmatrix *
create_random_sparse(const size_t M, const size_t N, const double density,
                     const gsl_rng *r)
{
  size_t nnzwanted = (size_t) floor(M * N * GSL_MIN(density, 1.0));
  gsl_spmatrix *m = gsl_spmatrix_alloc_nzmax(M, N,
                                             nnzwanted,
                                             GSL_SPMATRIX_TRIPLET);

  while (gsl_spmatrix_nnz(m) < nnzwanted)
    {
      /* generate a random row and column */
      size_t i = gsl_rng_uniform(r) * M;
      size_t j = gsl_rng_uniform(r) * N;

      /* generate random m_{ij} and add it */
      double x = gsl_rng_uniform(r);
      gsl_spmatrix_set(m, i, j, x);
    }

  return m;
} /* create_random_sparse() */

static gsl_spmatrix *
create_random_sparse_int(const size_t M, const size_t N, const double density,
                         const gsl_rng *r)
{
  const double lower = 1.0;
  const double upper = 10.0;
  size_t nnzwanted = (size_t) floor(M * N * GSL_MIN(density, 1.0));
  gsl_spmatrix *m = gsl_spmatrix_alloc_nzmax(M, N,
                                             nnzwanted,
                                             GSL_SPMATRIX_TRIPLET);

  while (gsl_spmatrix_nnz(m) < nnzwanted)
    {
      /* generate a random row and column */
      size_t i = gsl_rng_uniform(r) * M;
      size_t j = gsl_rng_uniform(r) * N;

      /* generate random m_{ij} and add it */
      int x = (int) (gsl_rng_uniform(r) * (upper - lower) + lower);
      gsl_spmatrix_set(m, i, j, (double) x);
    }

  return m;
}

static void
test_getset(const size_t M, const size_t N,
            const double density, const gsl_rng *r)
{
  int status;
  size_t i, j;

  /* test triplet versions of _get and _set */
  {
    size_t k = 0;
    gsl_spmatrix *m = gsl_spmatrix_alloc(M, N);

    status = 0;
    for (i = 0; i < M; ++i)
      {
        for (j = 0; j < N; ++j)
          {
            double x = (double) ++k;
            double y;

            gsl_spmatrix_set(m, i, j, x);
            y = gsl_spmatrix_get(m, i, j);
            if (x != y)
              status = 1;
          }
      }

    gsl_test(status, "test_getset: M=%zu N=%zu _get != _set", M, N);

    /* test setting an element to 0 */
    gsl_spmatrix_set(m, 0, 0, 1.0);
    gsl_spmatrix_set(m, 0, 0, 0.0);

    status = gsl_spmatrix_get(m, 0, 0) != 0.0;
    gsl_test(status, "test_getset: M=%zu N=%zu m(0,0) = %f",
             M, N, gsl_spmatrix_get(m, 0, 0));

    /* test gsl_spmatrix_set_zero() */
    gsl_spmatrix_set(m, 0, 0, 1.0);
    gsl_spmatrix_set_zero(m);
    status = gsl_spmatrix_get(m, 0, 0) != 0.0;
    gsl_test(status, "test_getset: M=%zu N=%zu set_zero m(0,0) = %f",
             M, N, gsl_spmatrix_get(m, 0, 0));

    /* resassemble matrix to ensure nz is calculated correctly */
    k = 0;
    for (i = 0; i < M; ++i)
      {
        for (j = 0; j < N; ++j)
          {
            double x = (double) ++k;
            gsl_spmatrix_set(m, i, j, x);
          }
      }

    status = gsl_spmatrix_nnz(m) != M * N;
    gsl_test(status, "test_getset: M=%zu N=%zu set_zero nz = %zu",
             M, N, gsl_spmatrix_nnz(m));

    gsl_spmatrix_free(m);
  }

  /* test duplicate values are handled correctly */
  {
    size_t min = GSL_MIN(M, N);
    size_t expected_nnz = min;
    size_t nnz;
    size_t k = 0;
    gsl_spmatrix *m = gsl_spmatrix_alloc(M, N);

    status = 0;
    for (i = 0; i < min; ++i)
      {
        for (j = 0; j < 5; ++j)
          {
            double x = (double) ++k;
            double y;

            gsl_spmatrix_set(m, i, i, x);
            y = gsl_spmatrix_get(m, i, i);
            if (x != y)
              status = 1;
          }
      }

    gsl_test(status, "test_getset: duplicate test M=%zu N=%zu _get != _set", M, N);

    nnz = gsl_spmatrix_nnz(m);
    status = nnz != expected_nnz;
    gsl_test(status, "test_getset: duplicate test M=%zu N=%zu nnz=%zu, expected=%zu",
             M, N, nnz, expected_nnz);

    gsl_spmatrix_free(m);
  }

  /* test CCS version of gsl_spmatrix_get() */
  {
    gsl_spmatrix *T = create_random_sparse(M, N, density, r);
    gsl_spmatrix *C = gsl_spmatrix_ccs(T);

    status = 0;
    for (i = 0; i < M; ++i)
      {
        for (j = 0; j < N; ++j)
          {
            double Tij = gsl_spmatrix_get(T, i, j);
            double Cij = gsl_spmatrix_get(C, i, j);

            if (Tij != Cij)
              status = 1;
          }
      }

    gsl_test(status, "test_getset: M=%zu N=%zu CCS _get", M, N);

    gsl_spmatrix_free(T);
    gsl_spmatrix_free(C);
  }

  /* test CRS version of gsl_spmatrix_get() */
  {
    gsl_spmatrix *T = create_random_sparse(M, N, density, r);
    gsl_spmatrix *C = gsl_spmatrix_crs(T);

    status = 0;
    for (i = 0; i < M; ++i)
      {
        for (j = 0; j < N; ++j)
          {
            double Tij = gsl_spmatrix_get(T, i, j);
            double Cij = gsl_spmatrix_get(C, i, j);

            if (Tij != Cij)
              status = 1;
          }
      }

    gsl_test(status, "test_getset: M=%zu N=%zu CRS _get", M, N);

    gsl_spmatrix_free(T);
    gsl_spmatrix_free(C);
  }
} /* test_getset() */

static void
test_memcpy(const size_t M, const size_t N,
            const double density, const gsl_rng *r)
{
  int status;

  {
    gsl_spmatrix *A = create_random_sparse(M, N, density, r);
    gsl_spmatrix *A_ccs = gsl_spmatrix_ccs(A);
    gsl_spmatrix *A_crs = gsl_spmatrix_crs(A);
    gsl_spmatrix *B_t, *B_ccs, *B_crs;
  
    B_t = gsl_spmatrix_alloc(M, N);
    gsl_spmatrix_memcpy(B_t, A);

    status = gsl_spmatrix_equal(A, B_t) != 1;
    gsl_test(status, "test_memcpy: _memcpy M=%zu N=%zu triplet format", M, N);

    B_ccs = gsl_spmatrix_alloc_nzmax(M, N, A_ccs->nzmax, GSL_SPMATRIX_CCS);
    B_crs = gsl_spmatrix_alloc_nzmax(M, N, A_ccs->nzmax, GSL_SPMATRIX_CRS);

    gsl_spmatrix_memcpy(B_ccs, A_ccs);
    gsl_spmatrix_memcpy(B_crs, A_crs);

    status = gsl_spmatrix_equal(A_ccs, B_ccs) != 1;
    gsl_test(status, "test_memcpy: _memcpy M=%zu N=%zu CCS", M, N);

    status = gsl_spmatrix_equal(A_crs, B_crs) != 1;
    gsl_test(status, "test_memcpy: _memcpy M=%zu N=%zu CRS", M, N);

    gsl_spmatrix_free(A);
    gsl_spmatrix_free(A_ccs);
    gsl_spmatrix_free(A_crs);
    gsl_spmatrix_free(B_t);
    gsl_spmatrix_free(B_ccs);
    gsl_spmatrix_free(B_crs);
  }

  /* test transpose_memcpy */
  {
    gsl_spmatrix *A = create_random_sparse(M, N, density, r);
    gsl_spmatrix *AT = gsl_spmatrix_alloc(N, M);
    gsl_spmatrix *B = gsl_spmatrix_ccs(A);
    gsl_spmatrix *BT = gsl_spmatrix_alloc_nzmax(N, M, 1, GSL_SPMATRIX_CCS);
    gsl_spmatrix *C = gsl_spmatrix_crs(A);
    gsl_spmatrix *CT = gsl_spmatrix_alloc_nzmax(N, M, 1, GSL_SPMATRIX_CRS);
    size_t i, j;

    gsl_spmatrix_transpose_memcpy(AT, A);
    gsl_spmatrix_transpose_memcpy(BT, B);
    gsl_spmatrix_transpose_memcpy(CT, C);

    status = 0;
    for (i = 0; i < M; ++i)
      {
        for (j = 0; j < N; ++j)
          {
            double Aij = gsl_spmatrix_get(A, i, j);
            double ATji = gsl_spmatrix_get(AT, j, i);

            if (Aij != ATji)
              status = 1;
          }
      }

    gsl_test(status, "test_memcpy: _transpose_memcpy M=%zu N=%zu triplet format", M, N);

    status = 0;
    for (i = 0; i < M; ++i)
      {
        for (j = 0; j < N; ++j)
          {
            double Aij = gsl_spmatrix_get(A, i, j);
            double Bij = gsl_spmatrix_get(B, i, j);
            double BTji = gsl_spmatrix_get(BT, j, i);

            if ((Bij != BTji) || (Aij != Bij))
              status = 1;
          }
      }

    gsl_test(status, "test_memcpy: _transpose_memcpy M=%zu N=%zu CCS format", M, N);

    status = 0;
    for (i = 0; i < M; ++i)
      {
        for (j = 0; j < N; ++j)
          {
            double Aij = gsl_spmatrix_get(A, i, j);
            double Cij = gsl_spmatrix_get(C, i, j);
            double CTji = gsl_spmatrix_get(CT, j, i);

            if ((Cij != CTji) || (Aij != Cij))
              status = 1;
          }
      }

    gsl_test(status, "test_memcpy: _transpose_memcpy M=%zu N=%zu CRS format", M, N);

    gsl_spmatrix_free(A);
    gsl_spmatrix_free(AT);
    gsl_spmatrix_free(B);
    gsl_spmatrix_free(BT);
    gsl_spmatrix_free(C);
    gsl_spmatrix_free(CT);
  }
} /* test_memcpy() */

static void
test_transpose(const size_t M, const size_t N,
               const double density, const gsl_rng *r)
{
  int status;
  gsl_spmatrix *A = create_random_sparse(M, N, density, r);
  gsl_spmatrix *AT = gsl_spmatrix_alloc_nzmax(M, N, A->nz, A->sptype);
  gsl_spmatrix *AT_ccs, *AT_crs;
  size_t i, j;

  /* test triplet transpose */

  gsl_spmatrix_memcpy(AT, A);
  gsl_spmatrix_transpose(AT);

  status = 0;
  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          double Aij = gsl_spmatrix_get(A, i, j);
          double ATji = gsl_spmatrix_get(AT, j, i);

          if (Aij != ATji)
            status = 1;
        }
    }

  gsl_test(status, "test_transpose: _transpose M=%zu N=%zu triplet format",
           M, N);

  /* test CCS transpose */

  AT_ccs = gsl_spmatrix_ccs(A);
  gsl_spmatrix_transpose(AT_ccs);

  status = 0;
  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          double Aij = gsl_spmatrix_get(A, i, j);
          double ATji = gsl_spmatrix_get(AT_ccs, j, i);

          if (Aij != ATji)
            status = 1;
        }
    }

  gsl_test(status, "test_transpose: _transpose M=%zu N=%zu CCS format",
           M, N);

  /* test CRS transpose */

  AT_crs = gsl_spmatrix_crs(A);
  gsl_spmatrix_transpose(AT_crs);

  status = 0;
  for (i = 0; i < M; ++i)
    {
      for (j = 0; j < N; ++j)
        {
          double Aij = gsl_spmatrix_get(A, i, j);
          double ATji = gsl_spmatrix_get(AT_crs, j, i);

          if (Aij != ATji)
            status = 1;
        }
    }

  gsl_test(status, "test_transpose: _transpose M=%zu N=%zu CRS format",
           M, N);

  gsl_spmatrix_free(A);
  gsl_spmatrix_free(AT);
  gsl_spmatrix_free(AT_ccs);
  gsl_spmatrix_free(AT_crs);
}

static void
test_ops(const size_t M, const size_t N,
         const double density, const gsl_rng *r)
{
  size_t i, j;
  int status;

  /* test gsl_spmatrix_add */
  {
    gsl_spmatrix *A = create_random_sparse(M, N, density, r);
    gsl_spmatrix *B = create_random_sparse(M, N, density, r);
    gsl_spmatrix *A_ccs = gsl_spmatrix_ccs(A);
    gsl_spmatrix *B_ccs = gsl_spmatrix_ccs(B);
    gsl_spmatrix *C_ccs = gsl_spmatrix_alloc_nzmax(M, N, 1, GSL_SPMATRIX_CCS);
    
    gsl_spmatrix_add(C_ccs, A_ccs, B_ccs);

    status = 0;
    for (i = 0; i < M; ++i)
      {
        for (j = 0; j < N; ++j)
          {
            double aij = gsl_spmatrix_get(A_ccs, i, j);
            double bij = gsl_spmatrix_get(B_ccs, i, j);
            double cij = gsl_spmatrix_get(C_ccs, i, j);

            if (aij + bij != cij)
              status = 1;
          }
      }

    gsl_test(status, "test_ops: _add M=%zu N=%zu CCS", M, N);

    gsl_spmatrix_free(A);
    gsl_spmatrix_free(B);
    gsl_spmatrix_free(A_ccs);
    gsl_spmatrix_free(B_ccs);
    gsl_spmatrix_free(C_ccs);
  }
} /* test_ops() */

static void
test_io(const size_t M, const size_t N,
        const double density, const gsl_rng *r)
{
  int status;
  gsl_spmatrix *A = create_random_sparse_int(M, N, density, r);

  char filename[] = "test.XXXXXX";
#if !defined( _WIN32 )
  int fd = mkstemp(filename);
#else
  char * fd = _mktemp(filename);
# define fdopen fopen
#endif

  /* test triplet I/O */
  {
    FILE *f = fdopen(fd, "w");

    gsl_spmatrix_fprintf(f, A, "%lg");

    fclose(f);
  }

  {
    FILE *f = fopen(filename, "r");
    gsl_spmatrix *B = gsl_spmatrix_alloc_nzmax(M, N, A->nz, A->sptype);

    gsl_spmatrix_fscanf(f, B);

    status = gsl_spmatrix_equal(A, B) != 1;
    gsl_test(status, "test_io: fprintf/fscanf M=%zu N=%zu triplet format", M, N);

    fclose(f);
    gsl_spmatrix_free(B);
  }

  unlink(filename);

  gsl_spmatrix_free(A);
}

int
main()
{
  gsl_rng *r = gsl_rng_alloc(gsl_rng_default);

  test_memcpy(10, 10, 0.2, r);
  test_memcpy(10, 15, 0.3, r);
  test_memcpy(53, 213, 0.4, r);
  test_memcpy(920, 2, 0.2, r);
  test_memcpy(2, 920, 0.3, r);

  test_getset(20, 20, 0.3, r);
  test_getset(30, 20, 0.3, r);
  test_getset(15, 210, 0.3, r);

  test_transpose(50, 50, 0.5, r);
  test_transpose(10, 40, 0.3, r);
  test_transpose(40, 10, 0.3, r);
  test_transpose(57, 13, 0.2, r);

  test_ops(20, 20, 0.2, r);
  test_ops(50, 20, 0.3, r);
  test_ops(20, 50, 0.3, r);
  test_ops(76, 43, 0.4, r);

  test_io(30, 30, 0.3, r);
  test_io(20, 10, 0.2, r);
  test_io(10, 20, 0.2, r);
  test_io(34, 78, 0.3, r);

  gsl_rng_free(r);

  exit (gsl_test_summary());
} /* main() */
