/* multimin/fminimizer.c
 */

#include <config.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_multimin.h>

gsl_multimin_fminimizer *
gsl_multimin_fminimizer_alloc (const gsl_multimin_fminimizer_type * T,
			       size_t n)
{
  int status;

  gsl_multimin_fminimizer *s =
    (gsl_multimin_fminimizer *) malloc (sizeof (gsl_multimin_fminimizer));

  if (s == 0)
    {
      GSL_ERROR_VAL ("failed to allocate space for minimizer struct",
		     GSL_ENOMEM, 0);
    }

  s->type = T;

  s->x = gsl_vector_calloc (n);

  if (s->x == 0) 
    {
      free (s);
      GSL_ERROR_VAL ("failed to allocate space for x", GSL_ENOMEM, 0);
    }

  s->state = malloc (T->size);

  if (s->state == 0)
    {
      gsl_vector_free (s->x);
      free (s);
      GSL_ERROR_VAL ("failed to allocate space for minimizer state",
		     GSL_ENOMEM, 0);
    }

  status = (T->alloc) (s->state, n);

  if (status != GSL_SUCCESS)
    {
      free (s->state);
      gsl_vector_free (s->x);
      free (s);

      GSL_ERROR_VAL ("failed to initialize minimizer state", GSL_ENOMEM, 0);
    }

  return s;
}

int
gsl_multimin_fminimizer_set (gsl_multimin_fminimizer * s,
			     gsl_multimin_function * f,
			     const gsl_vector * x,
			     const gsl_vector * step_size)
{
  if (s->x->size != f->n)
    {
      GSL_ERROR ("function incompatible with solver size", GSL_EBADLEN);
    }
  
  if (x->size != f->n || step_size->size != f->n) 
    {
      GSL_ERROR ("vector length not compatible with function", GSL_EBADLEN);
    }  
    
  s->f = f;

  gsl_vector_memcpy (s->x,x);

  return (s->type->set) (s->state, s->f, s->x, step_size);
}

void
gsl_multimin_fminimizer_free (gsl_multimin_fminimizer * s)
{
  (s->type->free) (s->state);
  free (s->state);
  gsl_vector_free (s->x);
  free (s);
}

int
gsl_multimin_fminimizer_iterate (gsl_multimin_fminimizer * s)
{
  return (s->type->iterate) (s->state, s->f, s->x, &(s->fval));
}

const char * 
gsl_multimin_fminimizer_name (const gsl_multimin_fminimizer * s)
{
  return s->type->name;
}


gsl_vector * 
gsl_multimin_fminimizer_x (gsl_multimin_fminimizer * s)
{
  return s->x;
}

double 
gsl_multimin_fminimizer_minimum (gsl_multimin_fminimizer * s)
{
  return s->fval;
}