/* steffenson.c -- steffenson root finding algorithm 

   This is Newton's method with an Aitken "delta-squared"
   acceleration of the iterates. This can improve the convergence on
   multiple roots where the ordinary Newton algorithm is slow.

   x[i+1] = x[i] - f(x[i]) / f'(x[i])

   x_accelerated[i] = x[i] - (x[i+1] - x[i])**2 / (x[i+2] - 2*x[i+1] - x[i])

   We can only use the accelerated estimate after three iterations,
   and use the unaccelerated value until then.

 */

#include <config.h>

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include <gsl_math.h>
#include <gsl_errno.h>
#include <gsl_roots.h>

#include "roots.h"

typedef struct
  {
    double f, df;
    double x;
    double x_1;
    double x_2;
    int count;
  }
steffenson_state_t;

int steffenson_init (void * vstate, gsl_function_fdf * fdf, double * root);
int steffenson_iterate (void * vstate, gsl_function_fdf * fdf, double * root);

int
steffenson_init (void * vstate, gsl_function_fdf * fdf, double * root)
{
  steffenson_state_t * state = (steffenson_state_t *) vstate;

  const double x = *root ;

  state->f = GSL_FN_FDF_EVAL_F (fdf, x);
  state->df = GSL_FN_FDF_EVAL_DF (fdf, x) ;

  state->x = x;
  state->x_1 = 0.0;
  state->x_2 = 0.0;

  state->count = 1;

  return GSL_SUCCESS;

}

int
steffenson_iterate (void * vstate, gsl_function_fdf * fdf, double * root)
{
  steffenson_state_t * state = (steffenson_state_t *) vstate;
  
  double x_new, f_new, df_new;

  double x_2 = state->x_2 ;
  double x_1 = state->x_1 ;
  double x = state->x ;

  if (state->df == 0.0)
    {
      GSL_ERROR("derivative is zero", GSL_EZERODIV);
    }

  x_new = x - (state->f / state->df);
  
  GSL_FN_FDF_EVAL_F_DF(fdf, x_new, &f_new, &df_new);

  state->x_2 = x_1 ;
  state->x_1 = x ;
  state->x = x_new;

  state->f = f_new ;
  state->df = df_new ;

  if (!GSL_IS_REAL (f_new))
    {
      GSL_ERROR ("function not continuous", GSL_EBADFUNC);
    }

  if (state->count < 3)
    {
      *root = x_new ;
      state->count++ ;
    }
  else 
    {
      double u = (x - x_1) ;
      *root = x_1 - u * u / (x_new - 2 * x + x_1) ;
    }

  if (!GSL_IS_REAL (df_new))
    {
      GSL_ERROR ("function not differentiable", GSL_EBADFUNC);
    }
      
  return GSL_SUCCESS;
}


static const gsl_root_fdfsolver_type steffenson_type =
{"steffenson",				/* name */
 sizeof (steffenson_state_t),
 &steffenson_init,
 &steffenson_iterate};

const gsl_root_fdfsolver_type  * gsl_root_fdfsolver_steffenson = &steffenson_type;
