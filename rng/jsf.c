/*
 *  rng/jsf.c
 * 
 *  Copyright (C) 2007 Bob Jenkins, placed into the Public domain.
 *  Copyright (C) 2020 Reini Urban, GSL packaging
 *
 *  Jenkin's Small Fast pseudo-random number generator 2007
 *  From https://burtleburtle.net/bob/rand/smallprng.html
 *
 *  The 32bit variant passes all dieharder tests, the 64bit variant fails some.
 */

#include <config.h>
#include <stdlib.h>
#include <stdint.h>
#include <gsl/gsl_rng.h>

static unsigned long int jsf_get (void *vstate);
static double jsf_get_double (void *vstate);
static void jsf_set (void *vstate, unsigned long int s);
static unsigned long int jsf64_get (void *vstate);
static double jsf64_get_double (void *vstate);
static void jsf64_set (void *vstate, unsigned long int s);

typedef struct {
 uint32_t a;
 uint32_t b;
 uint32_t c;
 uint32_t d;
} jsf_state_t;

static inline uint32_t rotl32(const uint32_t x, int k) {
  return (x << k) | (x >> (32 - k));
}

static unsigned long int jsf_get (void *vstate)
{

 jsf_state_t *x = vstate;
 uint32_t e;
 e    = x->a - rotl32(x->b, 27);
 x->a = x->b ^ rotl32(x->c, 17);
 x->b = x->c + x->d;
 x->c = x->d + e;
 x->d = e + x->a;
 return (unsigned long int)x->d;
}

static double jsf_get_double (void *vstate)
{
  return (double) jsf_get (vstate) / (double) UINT32_MAX;
}

/*
 (April 2009) Elias Yarrkov used a black-box solver to find some fixed
  points of this generator (which produce cycles of length 1):

  { 0x00000000, 0x00000000, 0x00000000, 0x00000000 },
  { 0x77777777, 0x55555555, 0x11111111, 0x44444444 },
  { 0x5591F2E3, 0x69EBA6CD, 0x2A171E3D, 0x3FD48890 },
  { 0x47CB8D56, 0xAE9B35A7, 0x5C78F4A8, 0x522240FF }

 (January 2016) David Blackman found two more, and claims that these
 six solutions are all that there are. I (Bob Jekins) have not verified.

  { 0x71AAC8F9, 0x66B4F5D3, 0x1E950B8F, 0x481FEA44 },
  { 0xAB23E5C6, 0xD3D74D9A, 0x542E3C7A, 0x7FA91120 }
 */
static void
jsf_set (void *vstate, unsigned long int s)
{
 jsf_state_t *x = (jsf_state_t *) vstate;
 int i;
 x->a = 0xf1ea5eed, x->b = x->c = x->d = s;
 for (i=0; i<20; ++i) {
   (void)jsf_get(x);
 } 
 return;
}

static const gsl_rng_type jsf_type =
{"jsf",			/* name */
 UINT32_MAX,			/* RAND_MAX */
 0,				/* RAND_MIN */
 sizeof (jsf_state_t),
 &jsf_set,
 &jsf_get,
 &jsf_get_double};

/* ============ 64bit variant ============= */

typedef struct {
 uint64_t a;
 uint64_t b;
 uint64_t c;
 uint64_t d;
} jsf64_state_t;

static inline uint64_t rotl64(const uint64_t x, int k) {
  return (x << k) | (x >> (64 - k));
}

static unsigned long int jsf64_get (void *vstate)
{

 jsf64_state_t *x = vstate;
 uint64_t e;
 e    = x->a - rotl64(x->b, 7);
 x->a = x->b ^ rotl64(x->c, 13);
 x->b = x->c + rotl64(x->d, 37);
 x->c = x->d + e;
 x->d = e + x->a;
 return (unsigned long int)x->d;
}

static double jsf64_get_double (void *vstate)
{
  return (double) ((jsf64_get (vstate) >> 11) * 0x1.0p-53);
}

static void
jsf64_set (void *vstate, unsigned long int s)
{
 /* Initialize automaton using specified seed. */
 jsf64_state_t *x = (jsf64_state_t *) vstate;
 int i;
 x->a = 0xf1ea5eed, x->b = x->c = x->d = s;
 for (i=0; i<20; ++i) {
   (void)jsf64_get(x);
 } 
 return;
}

static const gsl_rng_type jsf64_type =
{"jsf64",			/* name */
 UINT64_MAX,			/* RAND_MAX */
 0,				/* RAND_MIN */
 sizeof (jsf64_state_t),
 &jsf64_set,
 &jsf64_get,
 &jsf64_get_double};

const gsl_rng_type *gsl_rng_jsf = &jsf_type;
const gsl_rng_type *gsl_rng_jsf64 = &jsf64_type;
