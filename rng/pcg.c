/*
 *  rng/pcg.c
 *  From https://www.pcg-random.org/download.html
 *
 *  Copyright (C) 2014 M.E. O'Neill
 *  Copyright (C) 2020 Reini Urban, GSL packaging
 *
 * Licensed under Apache License 2.0 (NO WARRANTY, etc.
 * See http://www.apache.org/licenses/LICEMSE-2.0
 */

#include <config.h>
#include <stdlib.h>
#include <stdint.h>
#include <gsl/gsl_rng.h>

// *Really* minimal PCG32 code

typedef struct { uint64_t state;  uint64_t inc; } pcg32_random_t;

static unsigned long int pcg32_get(void *vstate)
{
  pcg32_random_t* rng = (pcg32_random_t*) vstate;
  uint64_t oldstate = rng->state;
  // Advance internal state
  rng->state = oldstate * 6364136223846793005ULL + (rng->inc|1);
  // Calculate output function (XSH RR), uses old state for max ILP
  uint32_t xorshifted = ((oldstate >> 18u) ^ oldstate) >> 27u;
  uint32_t rot = oldstate >> 59u;
  return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

// state for global RNGs
#define PCG32_INITIALIZER   { 0x853c49e6748fea9bULL, 0xda3e39cb94b95bdbULL }
static pcg32_random_t pcg32_global = PCG32_INITIALIZER;

// pcg32_srandom(initstate, initseq)
// pcg32_srandom_r(rng, initstate, initseq):
//     Seed the rng.  Specified in two parts, state initializer and a
//     sequence selection constant (a.k.a. stream id)

static void pcg32_srandom_r(pcg32_random_t* rng, uint64_t initstate, uint64_t initseq)
{
  rng->state = 0U;
  rng->inc = (initseq << 1u) | 1u;
  (void)pcg32_get(rng);
  rng->state += initstate;
  (void)pcg32_get(rng);
}

static void pcg32_set(void *vstate, unsigned long int seed)
{
  pcg32_random_t* rng = (pcg32_random_t*) vstate;
  pcg32_srandom_r(&pcg32_global, seed, (intptr_t)&rng);
}

static double
pcg32_get_double (void *vstate)
{
  return pcg32_get(vstate) / (double) UINT32_MAX;
}

static const gsl_rng_type pcg32_type =
{"pcg32",                       /* name */
 UINT32_MAX,			/* RAND_MAX */
 0,				/* RAND_MIN */
 sizeof (pcg32_random_t),
 &pcg32_set,
 &pcg32_get,
 &pcg32_get_double};

// pcg64
// pcg64_cmdxsm

const gsl_rng_type *gsl_rng_pcg32 = &pcg32_type;
//const gsl_rng_type *gsl_rng_pcg64 = &pcg64_type;
//const gsl_rng_type *gsl_rng_pcg64_cmdxsm = &pcg64_cmdxsm_type;
