.. index:: random number generators

************************
Random Number Generation
************************

.. include:: include.rst

The library provides a large collection of random number generators
which can be accessed through a uniform interface.  Environment
variables allow you to select different generators and seeds at runtime,
so that you can easily switch between generators without needing to
recompile your program.  Each instance of a generator keeps track of its
own state, allowing the generators to be used in multi-threaded
programs.  Additional functions are available for transforming uniform
random numbers into samples from continuous or discrete probability
distributions such as the Gaussian, log-normal or Poisson distributions.

These functions are declared in the header file :file:`gsl_rng.h`.

.. Need to explain the difference between SERIAL and PARALLEL random 
.. number generators here

General comments on random numbers
==================================

In 1988, Park and Miller wrote a paper entitled "Random number
generators: good ones are hard to find." [Commun.: ACM, 31, 1192--1201].
Fortunately, some excellent random number generators are available,
though poor ones are still in common use.  You may be happy with the
system-supplied random number generator on your computer, but you should
be aware that as computers get faster, requirements on random number
generators increase.  Nowadays, a simulation that calls a random number
generator millions of times can often finish before you can make it down
the hall to the coffee machine and back.

A very nice review of random number generators was written by Pierre
L'Ecuyer, as Chapter 4 of the book: Handbook on Simulation, Jerry Banks,
ed. (Wiley, 1997).  The chapter is available in postscript from
L'Ecuyer's ftp site (see references).  Knuth's volume on Seminumerical
Algorithms (originally published in 1968) devotes 170 pages to random
number generators, and has recently been updated in its 3rd edition
(1997).
It is brilliant, a classic.  If you don't own it, you should stop reading
right now, run to the nearest bookstore, and buy it.

A good random number generator will satisfy both theoretical and
statistical properties.  Theoretical properties are often hard to obtain
(they require real math!), but one prefers a random number generator
with a long period, low serial correlation, and a tendency **not** to
"fall mainly on the planes."  Statistical tests are performed with
numerical simulations.  Generally, a random number generator is used to
estimate some quantity for which the theory of probability provides an
exact answer.  Comparison to this exact answer provides a measure of
"randomness".

The Random Number Generator Interface
=====================================

It is important to remember that a random number generator is not a
"real" function like sine or cosine.  Unlike real functions, successive
calls to a random number generator yield different return values.  Of
course that is just what you want for a random number generator, but to
achieve this effect, the generator must keep track of some kind of
"state" variable.  Sometimes this state is just an integer (sometimes
just the value of the previously generated random number), but often it
is more complicated than that and may involve a whole array of numbers,
possibly with some indices thrown in.  To use the random number
generators, you do not need to know the details of what comprises the
state, and besides that varies from algorithm to algorithm.

.. type:: gsl_rng_type
          gsl_rng

   The random number generator library uses two special structs,
   :type:`gsl_rng_type` which holds static information about each type of
   generator and :type:`gsl_rng` which describes an instance of a generator
   created from a given :type:`gsl_rng_type`.

The functions described in this section are declared in the header file
:file:`gsl_rng.h`.

Random number generator initialization
======================================

.. function:: gsl_rng * gsl_rng_alloc (const gsl_rng_type * T)

   This function returns a pointer to a newly-created
   instance of a random number generator of type :data:`T`.
   For example, the following code creates an instance of the Tausworthe
   generator::

      gsl_rng * r = gsl_rng_alloc (gsl_rng_taus);

   If there is insufficient memory to create the generator then the
   function returns a null pointer and the error handler is invoked with an
   error code of :macro:`GSL_ENOMEM`.

   The generator is automatically initialized with the default seed,
   :data:`gsl_rng_default_seed`.  This is zero by default but can be changed
   either directly or by using the environment variable :macro:`GSL_RNG_SEED`.

   The details of the available generator types are
   described later in this chapter.

.. function:: void gsl_rng_set (const gsl_rng * r, unsigned long int s)

   This function initializes (or "seeds") the random number generator.  If
   the generator is seeded with the same value of :data:`s` on two different
   runs, the same stream of random numbers will be generated by successive
   calls to the routines below.  If different values of :math:`s \geq 1`
   are supplied, then the generated streams of random
   numbers should be completely different.  If the seed :data:`s` is zero
   then the standard seed from the original implementation is used
   instead.  For example, the original Fortran source code for the
   :code:`ranlux` generator used a seed of 314159265, and so choosing
   :data:`s` equal to zero reproduces this when using
   :data:`gsl_rng_ranlux`.

   When using multiple seeds with the same generator, choose seed values
   greater than zero to avoid collisions with the default setting.  

   Note that the most generators only accept 32-bit seeds, with higher
   values being reduced modulo :math:`2^{32}`.
   For generators with smaller ranges the maximum seed value will typically be lower.

.. function:: void gsl_rng_free (gsl_rng * r)

   This function frees all the memory associated with the generator
   :data:`r`.

Sampling from a random number generator
=======================================

The following functions return uniformly distributed random numbers,
either as integers or double precision floating point numbers.  |inlinefns|
To obtain non-uniform distributions, see :ref:`chap_random-number-distributions`.

.. function:: unsigned long int gsl_rng_get (const gsl_rng * r)

   This function returns a random integer from the generator :data:`r`.  The
   minimum and maximum values depend on the algorithm used, but all
   integers in the range [:data:`min`, :data:`max`] are equally likely.  The
   values of :data:`min` and :data:`max` can be determined using the auxiliary
   functions :func:`gsl_rng_max` and :func:`gsl_rng_min`.

.. function:: double gsl_rng_uniform (const gsl_rng * r)

   This function returns a double precision floating point number uniformly
   distributed in the range [0,1).  The range includes 0.0 but excludes 1.0.
   The value is typically obtained by dividing the result of
   :code:`gsl_rng_get(r)` by :code:`gsl_rng_max(r) + 1.0` in double
   precision.  Some generators compute this ratio internally so that they
   can provide floating point numbers with more than 32 bits of randomness
   (the maximum number of bits that can be portably represented in a single
   :code:`unsigned long int`).

.. function:: double gsl_rng_uniform_pos (const gsl_rng * r)

   This function returns a positive double precision floating point number
   uniformly distributed in the range (0,1), excluding both 0.0 and 1.0.
   The number is obtained by sampling the generator with the algorithm of
   :func:`gsl_rng_uniform` until a non-zero value is obtained.  You can use
   this function if you need to avoid a singularity at 0.0.

.. function:: unsigned long int gsl_rng_uniform_int (const gsl_rng * r, unsigned long int n)

   This function returns a random integer from 0 to :math:`n-1` inclusive
   by scaling down and/or discarding samples from the generator :data:`r`.
   All integers in the range :math:`[0,n-1]` are produced with equal
   probability.  For generators with a non-zero minimum value an offset
   is applied so that zero is returned with the correct probability.

   Note that this function is designed for sampling from ranges smaller
   than the range of the underlying generator.  The parameter :data:`n`
   must be less than or equal to the range of the generator :data:`r`.
   If :data:`n` is larger than the range of the generator then the function
   calls the error handler with an error code of :macro:`GSL_EINVAL` and
   returns zero.

   In particular, this function is not intended for generating the full range of
   unsigned integer values :math:`[0,2^{32}-1]`.
   Instead choose a generator with the maximal integer range and zero minimum
   value, such as :data:`gsl_rng_ranlxd1`, :data:`gsl_rng_mt19937` or
   :data:`gsl_rng_taus`, and sample it directly using
   :func:`gsl_rng_get`.  The range of each generator can be found using
   the auxiliary functions described in the next section.

Auxiliary random number generator functions
===========================================

The following functions provide information about an existing
generator.  You should use them in preference to hard-coding the generator
parameters into your own code.

.. function:: const char * gsl_rng_name (const gsl_rng * r)

   This function returns a pointer to the name of the generator.
   For example::

      printf ("r is a '%s' generator\n", gsl_rng_name (r));

   would print something like::
   
      r is a 'taus' generator

.. function:: unsigned long int gsl_rng_max (const gsl_rng * r)

   This function returns the largest value that :func:`gsl_rng_get`
   can return.

.. function:: unsigned long int gsl_rng_min (const gsl_rng * r)

   This function returns the smallest value that :func:`gsl_rng_get`
   can return.  Usually this value is zero.  There are some generators with
   algorithms that cannot return zero, and for these generators the minimum
   value is 1.

.. function:: void * gsl_rng_state (const gsl_rng * r)
              size_t gsl_rng_size (const gsl_rng * r)

   These functions return a pointer to the state of generator :data:`r` and
   its size.  You can use this information to access the state directly.  For
   example, the following code will write the state of a generator to a
   stream::

      void * state = gsl_rng_state (r);
      size_t n = gsl_rng_size (r);
      fwrite (state, n, 1, stream);

.. function:: const gsl_rng_type ** gsl_rng_types_setup (void)

   This function returns a pointer to an array of all the available
   generator types, terminated by a null pointer. The function should be
   called once at the start of the program, if needed.  The following code
   fragment shows how to iterate over the array of generator types to print
   the names of the available algorithms::

      const gsl_rng_type **t, **t0;

      t0 = gsl_rng_types_setup ();

      printf ("Available generators:\n");

      for (t = t0; *t != 0; t++)
        {
          printf ("%s\n", (*t)->name);
        }

Random number environment variables
===================================

The library allows you to choose a default generator and seed from the
environment variables :macro:`GSL_RNG_TYPE` and :macro:`GSL_RNG_SEED` and
the function :func:`gsl_rng_env_setup`.  This makes it easy try out
different generators and seeds without having to recompile your program.

.. macro:: GSL_RNG_TYPE

   This environment variable specifies the default random number generator.
   It should be the name of a generator, such as :code:`taus` or :code:`mt19937`.

.. macro:: GSL_RNG_SEED

   This environment variable specifies the default seed for the random
   number generator

.. var:: gsl_rng_type * gsl_rng_default

   This global library variable specifies the default random number generator,
   and can be initialized from :macro:`GSL_RNG_TYPE` using :func:`gsl_rng_env_setup`.
   It is defined as follows::

      extern const gsl_rng_type *gsl_rng_default

.. var:: unsigned long int gsl_rng_default_seed

   This global library variable specifies the seed for the default random number generator,
   and can be initialized from :macro:`GSL_RNG_SEED` using :func:`gsl_rng_env_setup`.
   It is set to zero by default and is defined as follows::

      extern unsigned long int gsl_rng_default_seed

.. function:: const gsl_rng_type * gsl_rng_env_setup (void)

   This function reads the environment variables :macro:`GSL_RNG_TYPE` and
   :macro:`GSL_RNG_SEED` and uses their values to set the corresponding
   library variables :data:`gsl_rng_default` and
   :data:`gsl_rng_default_seed`.

   The value of :macro:`GSL_RNG_SEED` is converted to an :code:`unsigned long int`
   using the C library function :func:`strtoul`.

   If you don't specify a generator for :macro:`GSL_RNG_TYPE` then
   :data:`gsl_rng_mt19937` is used as the default.  The initial value of
   :data:`gsl_rng_default_seed` is zero.

Here is a short program which shows how to create a global
generator using the environment variables :macro:`GSL_RNG_TYPE` and
:macro:`GSL_RNG_SEED`,

.. include:: examples/rng.c
   :code:

Running the program without any environment variables uses the initial
defaults, an :code:`mt19937` generator with a seed of 0,

.. include:: examples/rng.txt
   :code:

By setting the two variables on the command line we can
change the default generator and the seed::

  $ GSL_RNG_TYPE="taus" GSL_RNG_SEED=123 ./a.out 
  GSL_RNG_TYPE=taus
  GSL_RNG_SEED=123
  generator type: taus
  seed = 123
  first value = 2720986350

Copying random number generator state
=====================================

The above methods do not expose the random number state which changes
from call to call.  It is often useful to be able to save and restore
the state.  To permit these practices, a few somewhat more advanced
functions are supplied.  These include:

.. function:: int gsl_rng_memcpy (gsl_rng * dest, const gsl_rng * src)

   This function copies the random number generator :data:`src` into the
   pre-existing generator :data:`dest`, making :data:`dest` into an exact copy
   of :data:`src`.  The two generators must be of the same type.

.. function:: gsl_rng * gsl_rng_clone (const gsl_rng * r)

   This function returns a pointer to a newly created generator which is an
   exact copy of the generator :data:`r`.

Reading and writing random number generator state
=================================================

The library provides functions for reading and writing the random
number state to a file as binary data.

.. function:: int gsl_rng_fwrite (FILE * stream, const gsl_rng * r)

   This function writes the random number state of the random number
   generator :data:`r` to the stream :data:`stream` in binary format.  The
   return value is 0 for success and :macro:`GSL_EFAILED` if there was a
   problem writing to the file.  Since the data is written in the native
   binary format it may not be portable between different architectures.

.. function:: int gsl_rng_fread (FILE * stream, gsl_rng * r)

   This function reads the random number state into the random number
   generator :data:`r` from the open stream :data:`stream` in binary format.
   The random number generator :data:`r` must be preinitialized with the
   correct random number generator type since type information is not
   saved.  The return value is 0 for success and :macro:`GSL_EFAILED` if
   there was a problem reading from the file.  The data is assumed to
   have been written in the native binary format on the same
   architecture.

Random number generator algorithms
==================================

The functions described above make no reference to the actual algorithm
used.  This is deliberate so that you can switch algorithms without
having to change any of your application source code.  The library
provides a large number of generators of different types, including
simulation quality generators, generators provided for compatibility
with other libraries and historical generators from the past.

The following generators are recommended for use in simulation.  They
have extremely long periods, low correlation and pass most statistical
tests.  For the most reliable source of uncorrelated numbers, the
second-generation RANLUX generators have the strongest proof of
randomness.

.. index:: MT19937 random number generator

.. var:: gsl_rng_type * gsl_rng_mt19937

   The MT19937 generator of Makoto Matsumoto and Takuji Nishimura is a
   variant of the twisted generalized feedback shift-register algorithm,
   and is known as the "Mersenne Twister" generator.  It has a Mersenne
   prime period of :math:`2^{19937} - 1`
   (about :math:`10^{6000}`) and is
   equi-distributed in 623 dimensions.  It has passed the DIEHARD
   statistical tests.  It uses 624 words of state per generator and is
   comparable in speed to the other generators.  The original generator used
   a default seed of 4357 and choosing :data:`s` equal to zero in
   :func:`gsl_rng_set` reproduces this.  Later versions switched to 5489
   as the default seed, you can choose this explicitly via :func:`gsl_rng_set`
   instead if you require it.

   For more information see,

   * Makoto Matsumoto and Takuji Nishimura, "Mersenne Twister: A
     623-dimensionally equidistributed uniform pseudorandom number
     generator". ACM Transactions on Modeling and Computer
     Simulation, Vol.: 8, No.: 1 (Jan. 1998), Pages 3--30

   The generator :data:`gsl_rng_mt19937` uses the second revision of the
   seeding procedure published by the two authors above in 2002.  The
   original seeding procedures could cause spurious artifacts for some seed
   values. They are still available through the alternative generators
   :data:`gsl_rng_mt19937_1999` and :data:`gsl_rng_mt19937_1998`.

.. index:: RANLXS random number generator

.. var:: gsl_rng_type * gsl_rng_ranlxs0
         gsl_rng_type * gsl_rng_ranlxs1
         gsl_rng_type * gsl_rng_ranlxs2

   The generator :code:`ranlxs0` is a second-generation version of the
   RANLUX algorithm of Luscher, which produces "luxury random
   numbers".  This generator provides single precision output (24 bits) at
   three luxury levels :code:`ranlxs0`, :code:`ranlxs1` and :code:`ranlxs2`,
   in increasing order of strength.  
   It uses double-precision floating point arithmetic internally and can be
   significantly faster than the integer version of :code:`ranlux`,
   particularly on 64-bit architectures.  The period of the generator is
   about :math:`10^{171}`.
   The algorithm has mathematically proven properties and
   can provide truly decorrelated numbers at a known level of randomness.
   The higher luxury levels provide increased decorrelation between samples
   as an additional safety margin.

   Note that the range of allowed seeds for this generator is :math:`[0,2^{31}-1]`.
   Higher seed values are wrapped modulo :math:`2^{31}`.

.. index:: RANLXD random number generator

.. var:: gsl_rng_type * gsl_rng_ranlxd1
         gsl_rng_type * gsl_rng_ranlxd2

   These generators produce double precision output (48 bits) from the
   RANLXS generator.  The library provides two luxury levels
   :code:`ranlxd1` and :code:`ranlxd2`, in increasing order of strength.

.. index:: RANLUX random number generator

.. var:: gsl_rng_type * gsl_rng_ranlux
         gsl_rng_type * gsl_rng_ranlux389

   The :code:`ranlux` generator is an implementation of the original
   algorithm developed by Luscher.  It uses a
   lagged-fibonacci-with-skipping algorithm to produce "luxury random
   numbers".  It is a 24-bit generator, originally designed for
   single-precision IEEE floating point numbers.  This implementation is
   based on integer arithmetic, while the second-generation versions
   RANLXS and RANLXD described above provide floating-point
   implementations which will be faster on many platforms.
   The period of the generator is about :math:`10^{171}`.
   The algorithm has mathematically proven properties and
   it can provide truly decorrelated numbers at a known level of
   randomness.  The default level of decorrelation recommended by Luscher
   is provided by :data:`gsl_rng_ranlux`, while :data:`gsl_rng_ranlux389`
   gives the highest level of randomness, with all 24 bits decorrelated.
   Both types of generator use 24 words of state per generator.

   For more information see,

   * M. Luscher, "A portable high-quality random number generator for
     lattice field theory calculations", Computer Physics
     Communications, 79 (1994) 100--110.

   * F. James, "RANLUX: A Fortran implementation of the high-quality
     pseudo-random number generator of Luscher", Computer Physics
     Communications, 79 (1994) 111--114

.. index::
   single: CMRG, combined multiple recursive random number generator

.. var:: gsl_rng_type * gsl_rng_cmrg

   This is a combined multiple recursive generator by L'Ecuyer. 
   Its sequence is,

   .. math:: z_n = (x_n - y_n) \mod m_1

   where the two underlying generators :math:`x_n` and :math:`y_n` are,

   .. only:: not texinfo

      .. math::

         x_n & = (a_1 x_{n-1} + a_2 x_{n-2} + a_3 x_{n-3}) \mod m_1 \\
         y_n & = (b_1 y_{n-1} + b_2 y_{n-2} + b_3 y_{n-3}) \mod m_2

   .. only:: texinfo

      ::

         x_n = (a_1 x_{n-1} + a_2 x_{n-2} + a_3 x_{n-3}) mod m_1
         y_n = (b_1 y_{n-1} + b_2 y_{n-2} + b_3 y_{n-3}) mod m_2

   with coefficients 
   :math:`a_1 = 0`, 
   :math:`a_2 = 63308`, 
   :math:`a_3 = -183326`,
   :math:`b_1 = 86098`, 
   :math:`b_2 = 0`,
   :math:`b_3 = -539608`,
   and moduli 
   :math:`m_1 = 2^{31} - 1 = 2147483647`
   and 
   :math:`m_2 = 2145483479`.

   The period of this generator is  
   :math:`\hbox{lcm}(m_1^3-1, m_2^3-1)`,
   which is approximately
   :math:`2^{185}`
   (about :math:`10^{56}`).
   It uses 6 words of state per generator.  For more information see,

   * P. L'Ecuyer, "Combined Multiple Recursive Random Number
     Generators", Operations Research, 44, 5 (1996), 816--822.

.. index::
   single: MRG, multiple recursive random number generator

.. var:: gsl_rng_type * gsl_rng_mrg

   This is a fifth-order multiple recursive generator by L'Ecuyer, Blouin
   and Coutre.  Its sequence is,

   .. math:: x_n = (a_1 x_{n-1} + a_5 x_{n-5}) \mod m

   with 
   :math:`a_1 = 107374182`,
   :math:`a_2 = a_3 = a_4 = 0`, 
   :math:`a_5 = 104480`
   and 
   :math:`m = 2^{31}-1`.

   The period of this generator is about 
   :math:`10^{46}`.
   It uses 5 words
   of state per generator.  More information can be found in the following
   paper,

   * P. L'Ecuyer, F. Blouin, and R. Coutre, "A search for good multiple
     recursive random number generators", ACM Transactions on Modeling and
     Computer Simulation 3, 87--98 (1993).

.. index:: Tausworthe random number generator

.. var:: gsl_rng_type * gsl_rng_taus
         gsl_rng_type * gsl_rng_taus2

   This is a maximally equidistributed combined Tausworthe generator by
   L'Ecuyer.  The sequence is,

   .. only:: not texinfo

      .. math:: x_n = (s^1_n \oplus s^2_n \oplus s^3_n) 

   .. only:: texinfo

      ::

         x_n = (s1_n ^^ s2_n ^^ s3_n) 

   where,

   .. only:: not texinfo

      .. math::

         s^1_{n+1} &= (((s^1_n \& 4294967294)\ll 12) \oplus (((s^1_n\ll 13) \oplus s^1_n)\gg 19)) \\
         s^2_{n+1} &= (((s^2_n \& 4294967288)\ll 4) \oplus (((s^2_n\ll 2) \oplus s^2_n)\gg 25)) \\
         s^3_{n+1} &= (((s^3_n \& 4294967280)\ll 17) \oplus (((s^3_n\ll 3) \oplus s^3_n)\gg 11))

   .. only:: texinfo

      ::

         s1_{n+1} = (((s1_n&4294967294)<<12)^^(((s1_n<<13)^^s1_n)>>19))
         s2_{n+1} = (((s2_n&4294967288)<< 4)^^(((s2_n<< 2)^^s2_n)>>25))
         s3_{n+1} = (((s3_n&4294967280)<<17)^^(((s3_n<< 3)^^s3_n)>>11))

   computed modulo 
   :math:`2^{32}`.
   In the formulas above 
   :math:`\oplus`
   denotes *exclusive-or*.  Note that the algorithm relies on the properties
   of 32-bit unsigned integers and has been implemented using a bitmask
   of :code:`0xFFFFFFFF` to make it work on 64 bit machines.

   The period of this generator is :math:`2^{88}`
   (about :math:`10^{26}`).
   It uses 3 words of state per generator.  For more
   information see,

   * P. L'Ecuyer, "Maximally Equidistributed Combined Tausworthe
     Generators", Mathematics of Computation, 65, 213 (1996), 203--213.

   The generator :data:`gsl_rng_taus2` uses the same algorithm as
   :data:`gsl_rng_taus` but with an improved seeding procedure described in
   the paper,

   * P. L'Ecuyer, "Tables of Maximally Equidistributed Combined LFSR
     Generators", Mathematics of Computation, 68, 225 (1999), 261--269

   The generator :data:`gsl_rng_taus2` should now be used in preference to
   :data:`gsl_rng_taus`.

.. index:: Four-tap Generalized Feedback Shift Register

.. var:: gsl_rng_type * gsl_rng_gfsr4

   The :code:`gfsr4` generator is like a lagged-fibonacci generator, and 
   produces each number as an :code:`xor`'d sum of four previous values.

   .. only:: not texinfo

      .. math:: r_n = r_{n-A} \oplus r_{n-B} \oplus r_{n-C} \oplus r_{n-D}

   .. only:: texinfo

      ::

         r_n = r_{n-A} ^^ r_{n-B} ^^ r_{n-C} ^^ r_{n-D}

   Ziff (ref below) notes that "it is now widely known" that two-tap
   registers (such as R250, which is described below)
   have serious flaws, the most obvious one being the three-point
   correlation that comes from the definition of the generator.  Nice
   mathematical properties can be derived for GFSR's, and numerics bears
   out the claim that 4-tap GFSR's with appropriately chosen offsets are as
   random as can be measured, using the author's test.

   This implementation uses the values suggested the example on p392 of
   Ziff's article: :math:`A=471`, :math:`B=1586`, :math:`C=6988`, :math:`D=9689`.

   If the offsets are appropriately chosen (such as the one ones in this
   implementation), then the sequence is said to be maximal; that means
   that the period is :math:`2^D - 1`, where :math:`D` is the longest lag.
   (It is one less than :math:`2^D` because it is not permitted to have all
   zeros in the :code:`ra[]` array.)  For this implementation with
   :math:`D=9689` that works out to about :math:`10^{2917}`.

   Note that the implementation of this generator using a 32-bit
   integer amounts to 32 parallel implementations of one-bit
   generators.  One consequence of this is that the period of this
   32-bit generator is the same as for the one-bit generator.
   Moreover, this independence means that all 32-bit patterns are
   equally likely, and in particular that 0 is an allowed random
   value.  (We are grateful to Heiko Bauke for clarifying for us these
   properties of GFSR random number generators.)

   For more information see,

   * Robert M. Ziff, "Four-tap shift-register-sequence random-number 
     generators", Computers in Physics, 12(4), Jul/Aug
     1998, pp 385--392.

Unix random number generators
=============================

The standard Unix random number generators :code:`rand`, :code:`random`
and :code:`rand48` are provided as part of GSL. Although these
generators are widely available individually often they aren't all
available on the same platform.  This makes it difficult to write
portable code using them and so we have included the complete set of
Unix generators in GSL for convenience.  Note that these generators
don't produce high-quality randomness and aren't suitable for work
requiring accurate statistics.  However, if you won't be measuring
statistical quantities and just want to introduce some variation into
your program then these generators are quite acceptable.

.. index::
   single: rand, BSD random number generator
   single: Unix random number generators, rand
   single: Unix random number generators, rand48

.. index:: BSD random number generator

.. var:: gsl_rng_type * gsl_rng_rand

   This is the BSD :code:`rand` generator.  Its sequence is

   .. math:: x_{n+1} = (a x_n + c) \mod m

   with 
   :math:`a = 1103515245`, 
   :math:`c = 12345` and 
   :math:`m = 2^{31}`.
   The seed specifies the initial value, 
   :math:`x_1`.  The period of this
   generator is 
   :math:`2^{31}`,
   and it uses 1 word of storage per generator.

.. var:: gsl_rng_type * gsl_rng_random_bsd
         gsl_rng_type * gsl_rng_random_libc5
         gsl_rng_type * gsl_rng_random_glibc2

   These generators implement the :code:`random` family of functions, a
   set of linear feedback shift register generators originally used in BSD
   Unix.  There are several versions of :code:`random` in use today: the
   original BSD version (e.g. on SunOS4), a libc5 version (found on
   older GNU/Linux systems) and a glibc2 version.  Each version uses a
   different seeding procedure, and thus produces different sequences.

   The original BSD routines accepted a variable length buffer for the
   generator state, with longer buffers providing higher-quality
   randomness.  The :code:`random` function implemented algorithms for
   buffer lengths of 8, 32, 64, 128 and 256 bytes, and the algorithm with
   the largest length that would fit into the user-supplied buffer was
   used.  To support these algorithms additional generators are available
   with the following names::

      gsl_rng_random8_bsd
      gsl_rng_random32_bsd
      gsl_rng_random64_bsd
      gsl_rng_random128_bsd
      gsl_rng_random256_bsd

   where the numeric suffix indicates the buffer length.  The original BSD
   :code:`random` function used a 128-byte default buffer and so
   :data:`gsl_rng_random_bsd` has been made equivalent to
   :data:`gsl_rng_random128_bsd`.  Corresponding versions of the :code:`libc5`
   and :code:`glibc2` generators are also available, with the names
   :data:`gsl_rng_random8_libc5`, :data:`gsl_rng_random8_glibc2`, etc.

.. index:: rand48 random number generator

.. var:: gsl_rng_type * gsl_rng_rand48

   This is the Unix :code:`rand48` generator.  Its sequence is

   .. math:: x_{n+1} = (a x_n + c) \mod m

   defined on 48-bit unsigned integers with 
   :math:`a = 25214903917`, 
   :math:`c = 11` and 
   :math:`m = 2^{48}`.
   The seed specifies the upper 32 bits of the initial value, :math:`x_1`,
   with the lower 16 bits set to :code:`0x330E`.  The function
   :func:`gsl_rng_get` returns the upper 32 bits from each term of the
   sequence.  This does not have a direct parallel in the original
   :code:`rand48` functions, but forcing the result to type :code:`long int`
   reproduces the output of :code:`mrand48`.  The function
   :func:`gsl_rng_uniform` uses the full 48 bits of internal state to return
   the double precision number :math:`x_n/m`, which is equivalent to the
   function :code:`drand48`.  Note that some versions of the GNU C Library
   contained a bug in :code:`mrand48` function which caused it to produce
   different results (only the lower 16-bits of the return value were set).

Other random number generators
==============================

The generators in this section are provided for compatibility with
existing libraries.  If you are converting an existing program to use GSL
then you can select these generators to check your new implementation
against the original one, using the same random number generator.  After
verifying that your new program reproduces the original results you can
then switch to a higher-quality generator.

Note that most of the generators in this section are based on single
linear congruence relations, which are the least sophisticated type of
generator.  In particular, linear congruences have poor properties when
used with a non-prime modulus, as several of these routines do (e.g.
with a power of two modulus, 
:math:`2^{31}` or
:math:`2^{32}`).
This leads to periodicity in the least significant bits of each number,
with only the higher bits having any randomness.  Thus if you want to
produce a random bitstream it is best to avoid using the least
significant bits.

.. index::
   single: RANF random number generator
   single: CRAY random number generator, RANF

.. var:: gsl_rng_type * gsl_rng_ranf

   This is the CRAY random number generator :code:`RANF`.  Its sequence is

   .. math:: x_{n+1} = (a x_n) \mod m

   defined on 48-bit unsigned integers with :math:`a = 44485709377909` and
   :math:`m = 2^{48}`.
   The seed specifies the lower 32 bits of the initial value, 
   :math:`x_1`, with the lowest bit set to
   prevent the seed taking an even value.  The upper 16 bits of 
   :math:`x_1`
   are set to 0. A consequence of this procedure is that the pairs of seeds
   2 and 3, 4 and 5, etc.: produce the same sequences.

   The generator compatible with the CRAY MATHLIB routine RANF. It
   produces double precision floating point numbers which should be
   identical to those from the original RANF.

   There is a subtlety in the implementation of the seeding.  The initial
   state is reversed through one step, by multiplying by the modular
   inverse of :math:`a` mod :math:`m`.  This is done for compatibility with
   the original CRAY implementation.

   Note that you can only seed the generator with integers up to
   :math:`2^{32}`,
   while the original CRAY implementation uses
   non-portable wide integers which can cover all 
   :math:`2^{48}`
   states of the generator.

   The function :func:`gsl_rng_get` returns the upper 32 bits from each term
   of the sequence.  The function :func:`gsl_rng_uniform` uses the full 48
   bits to return the double precision number :math:`x_n/m`.

   The period of this generator is :math:`2^{46}`.

.. index:: RANMAR random number generator

.. var:: gsl_rng_type * gsl_rng_ranmar

   This is the RANMAR lagged-fibonacci generator of Marsaglia, Zaman and
   Tsang.  It is a 24-bit generator, originally designed for
   single-precision IEEE floating point numbers.  It was included in the
   CERNLIB high-energy physics library.

.. index::
   single: shift-register random number generator
   single: R250 shift-register random number generator

.. var:: gsl_rng_type * gsl_rng_r250

   This is the shift-register generator of Kirkpatrick and Stoll.  The
   sequence is based on the recurrence

   .. only:: not texinfo

      .. math:: x_n = x_{n-103} \oplus x_{n-250}

   .. only:: texinfo

      ::

         x_n = x_{n-103} ^^ x_{n-250}

   where 
   :math:`\oplus`
   denotes *exclusive-or*, defined on
   32-bit words.  The period of this generator is about :math:`2^{250}` and it
   uses 250 words of state per generator.

   For more information see,

   * S. Kirkpatrick and E. Stoll, "A very fast shift-register sequence random
     number generator", Journal of Computational Physics, 40, 517--526
     (1981)

.. index:: TT800 random number generator

.. var:: gsl_rng_type * gsl_rng_tt800

   This is an earlier version of the twisted generalized feedback
   shift-register generator, and has been superseded by the development of
   MT19937.  However, it is still an acceptable generator in its own
   right.  It has a period of 
   :math:`2^{800}`
   and uses 33 words of storage per generator.

   For more information see,

   * Makoto Matsumoto and Yoshiharu Kurita, "Twisted GFSR Generators
     II", ACM Transactions on Modelling and Computer Simulation,
     Vol.: 4, No.: 3, 1994, pages 254--266.

.. The following generators are included only for historical reasons, so
.. that you can reproduce results from old programs which might have used
.. them.  These generators should not be used for real simulations since
.. they have poor statistical properties by modern standards.

.. index:: VAX random number generator

.. var:: gsl_rng_type * gsl_rng_vax

   This is the VAX generator :code:`MTH$RANDOM`.  Its sequence is,

   .. math:: x_{n+1} = (a x_n + c) \mod m

   with 
   :math:`a = 69069`, :math:`c = 1` and 
   :math:`m = 2^{32}`.
   The seed specifies the initial value, 
   :math:`x_1`.  The
   period of this generator is 
   :math:`2^{32}`
   and it uses 1 word of storage per
   generator.

.. var:: gsl_rng_type * gsl_rng_transputer

   This is the random number generator from the INMOS Transputer
   Development system.  Its sequence is,

   .. math:: x_{n+1} = (a x_n) \mod m

   with :math:`a = 1664525` and 
   :math:`m = 2^{32}`.
   The seed specifies the initial value, 
   :math:`x_1`.

.. index:: RANDU random number generator

.. var:: gsl_rng_type * gsl_rng_randu

   This is the IBM :code:`RANDU` generator.  Its sequence is

   .. math:: x_{n+1} = (a x_n) \mod m

   with :math:`a = 65539` and 
   :math:`m = 2^{31}`. The
   seed specifies the initial value, 
   :math:`x_1`.  The period of this
   generator was only 
   :math:`2^{29}`.
   It has become a textbook example of a poor generator.

.. index:: RANMAR random number generator

.. var:: gsl_rng_type * gsl_rng_minstd

   This is Park and Miller's "minimal standard" MINSTD generator, a
   simple linear congruence which takes care to avoid the major pitfalls of
   such algorithms.  Its sequence is,

   .. math:: x_{n+1} = (a x_n) \mod m

   with :math:`a = 16807` and 
   :math:`m = 2^{31} - 1 = 2147483647`.
   The seed specifies the initial value, 
   :math:`x_1`.  The period of this
   generator is about 
   :math:`2^{31}`.

   This generator was used in the IMSL Library (subroutine RNUN) and in
   MATLAB (the RAND function) in the past.  It is also sometimes known by
   the acronym "GGL" (I'm not sure what that stands for).

   For more information see,

   * Park and Miller, "Random Number Generators: Good ones are hard to find",
     Communications of the ACM, October 1988, Volume 31, No 10, pages
     1192--1201.

.. var:: gsl_rng_type * gsl_rng_uni
         gsl_rng_type * gsl_rng_uni32

   This is a reimplementation of the 16-bit SLATEC random number generator
   RUNIF. A generalization of the generator to 32 bits is provided by
   :data:`gsl_rng_uni32`.  The original source code is available from NETLIB.

.. var:: gsl_rng_type * gsl_rng_slatec

   This is the SLATEC random number generator RAND. It is ancient.  The
   original source code is available from NETLIB.

.. var:: gsl_rng_type * gsl_rng_zuf

   This is the ZUFALL lagged Fibonacci series generator of Peterson.  Its
   sequence is,

   .. only:: not texinfo

      .. math::

         t &= u_{n-273} + u_{n-607} \\
         u_n  &= t - \hbox{floor}(t)

   .. only:: texinfo

      ::

         t = u_{n-273} + u_{n-607}
         u_n  = t - floor(t)

   The original source code is available from NETLIB.  For more information
   see,

   * W. Petersen, "Lagged Fibonacci Random Number Generators for the NEC
     SX-3", International Journal of High Speed Computing (1994).

.. var:: gsl_rng_type * gsl_rng_knuthran2

   This is a second-order multiple recursive generator described by Knuth
   in Seminumerical Algorithms, 3rd Ed., page 108.  Its sequence is,

   .. math:: x_n = (a_1 x_{n-1} + a_2 x_{n-2}) \mod m

   with 
   :math:`a_1 = 271828183`, 
   :math:`a_2 = 314159269`, 
   and 
   :math:`m = 2^{31}-1`.

.. var:: gsl_rng_type * gsl_rng_knuthran2002
         gsl_rng_type * gsl_rng_knuthran

   This is a second-order multiple recursive generator described by Knuth
   in Seminumerical Algorithms, 3rd Ed., Section 3.6.  Knuth
   provides its C code.  The updated routine :data:`gsl_rng_knuthran2002`
   is from the revised 9th printing and corrects some weaknesses in the
   earlier version, which is implemented as :data:`gsl_rng_knuthran`.

.. var:: gsl_rng_type * gsl_rng_borosh13
         gsl_rng_type * gsl_rng_fishman18
         gsl_rng_type * gsl_rng_fishman20
         gsl_rng_type * gsl_rng_lecuyer21
         gsl_rng_type * gsl_rng_waterman14

   These multiplicative generators are taken from Knuth's
   Seminumerical Algorithms, 3rd Ed., pages 106--108. Their sequence
   is,

   .. math:: x_{n+1} = (a x_n) \mod m

   where the seed specifies the initial value,
   :math:`x_1`.
   The parameters :math:`a` and :math:`m` are as follows,
   Borosh-Niederreiter: 
   :math:`a = 1812433253`, :math:`m = 2^{32}`,
   Fishman18:
   :math:`a = 62089911`,
   :math:`m = 2^{31}-1`,
   Fishman20:
   :math:`a = 48271`,
   :math:`m = 2^{31}-1`,
   L'Ecuyer:
   :math:`a = 40692`,
   :math:`m = 2^{31}-249`,
   Waterman:
   :math:`a = 1566083941`,
   :math:`m = 2^{32}`.

.. var:: gsl_rng_type * gsl_rng_fishman2x

   This is the L'Ecuyer--Fishman random number generator. It is taken from
   Knuth's Seminumerical Algorithms, 3rd Ed., page 108. Its sequence
   is,

   .. math:: z_{n+1} = (x_n - y_n) \mod m

   with :math:`m = 2^{31}-1`.
   :math:`x_n` and :math:`y_n` are given by the :code:`fishman20` 
   and :code:`lecuyer21` algorithms.
   The seed specifies the initial value, 
   :math:`x_1`.

.. var:: gsl_rng_type * gsl_rng_coveyou

   This is the Coveyou random number generator. It is taken from Knuth's
   Seminumerical Algorithms, 3rd Ed., Section 3.2.2. Its sequence
   is,

   .. math:: x_{n+1} = (x_n (x_n + 1)) \mod m

   with :math:`m = 2^{32}`.
   The seed specifies the initial value, 
   :math:`x_1`.

.. index:: JSF random number generator

.. var:: gsl_rng_jsf

   This is Bob Jenkin's small pseudo-random number generator (2007)
   from https://burtleburtle.net/bob/rand/smallprng.html

   The state needs to contain 4 random 32-bit words, with those
   forbidden fixpoints:

       { 0x00000000, 0x00000000, 0x00000000, 0x00000000 },
       { 0x77777777, 0x55555555, 0x11111111, 0x44444444 },
       { 0x5591F2E3, 0x69EBA6CD, 0x2A171E3D, 0x3FD48890 },
       { 0x47CB8D56, 0xAE9B35A7, 0x5C78F4A8, 0x522240FF }
       { 0x71AAC8F9, 0x66B4F5D3, 0x1E950B8F, 0x481FEA44 },
       { 0xAB23E5C6, 0xD3D74D9A, 0x542E3C7A, 0x7FA91120 }

   The given set function avoids forbidden states for all given seeds.
   
   The period of this generator is about :math:`2^{126}`, with a shortest
   cycle of :math:`2^{94}` and it uses 4 words of state per generator.
   The 32bit variant passes all DIEHARDER tests, the 64bit variant
   `gsl_rng_jsf64` fails some.


Performance
===========

.. I made the original plot like this
.. ./benchmark > tmp; cat tmp | perl -n -e '($n,$s) = split(" ",$_); printf("%17s ",$n); print "-" x ($s/1e5), "\n";'

.. The large number of generators based on single linear congruences are
.. represented by the :code:`random` generator below.  These generators are
.. fast but have the lowest statistical quality.

The following table shows the relative performance of a all
available random number generators.  The fastest simulation quality
generators are :code:`taus`, :code:`gfsr4` and :code:`mt19937`.  The
generators which offer the best mathematically-proven quality are those
based on the RANLUX algorithm::

  328964 k ints/sec, 309224 k doubles/sec, transputer
  328636 k ints/sec, 309094 k doubles/sec, waterman14
  319741 k ints/sec, 297242 k doubles/sec, borosh13
  319059 k ints/sec, 271855 k doubles/sec, ran3
  314252 k ints/sec, 281049 k doubles/sec, r250
  310345 k ints/sec, 293818 k doubles/sec, jsf
  306351 k ints/sec, 295357 k doubles/sec, xoroshiro64**
  306348 k ints/sec, 294652 k doubles/sec, xoshiro256+
  298077 k ints/sec, 298469 k doubles/sec, jsf64
  296997 k ints/sec, 268763 k doubles/sec, xoroshiro64*
  295887 k ints/sec, 299276 k doubles/sec, coveyou
  295622 k ints/sec, 297637 k doubles/sec, xoroshiro128+
  295145 k ints/sec, 298620 k doubles/sec, vax
  293844 k ints/sec, 136665 k doubles/sec, ranmar
  293525 k ints/sec, 289265 k doubles/sec, random-glibc2
  293315 k ints/sec, 297062 k doubles/sec, randu
  293195 k ints/sec, 288070 k doubles/sec, random64-glibc2
  293123 k ints/sec, 288626 k doubles/sec, random-bsd
  292932 k ints/sec, 289939 k doubles/sec, random128-bsd
  291894 k ints/sec, 288411 k doubles/sec, random-libc5
  291792 k ints/sec, 287893 k doubles/sec, random64-libc5
  291670 k ints/sec, 272732 k doubles/sec, random32-libc5
  291389 k ints/sec, 288216 k doubles/sec, random32-glibc2
  290627 k ints/sec, 277519 k doubles/sec, uni
  289957 k ints/sec, 276028 k doubles/sec, random256-glibc2
  289719 k ints/sec, 278012 k doubles/sec, random256-bsd
  286423 k ints/sec, 274391 k doubles/sec, xoshiro256++
  286205 k ints/sec, 278385 k doubles/sec, random256-libc5
  285164 k ints/sec, 279405 k doubles/sec, xoshiro128+
  284937 k ints/sec, 282090 k doubles/sec, random128-glibc2
  284853 k ints/sec, 267730 k doubles/sec, random64-bsd
  284797 k ints/sec, 269101 k doubles/sec, uni32
  283980 k ints/sec, 285506 k doubles/sec, random8-glibc2
  283802 k ints/sec, 280524 k doubles/sec, random8-bsd
  283488 k ints/sec, 283380 k doubles/sec, xoroshiro128**
  283015 k ints/sec, 285540 k doubles/sec, random8-libc5
  282902 k ints/sec, 285138 k doubles/sec, rand
  282296 k ints/sec, 282613 k doubles/sec, random32-bsd
  280662 k ints/sec, 289208 k doubles/sec, xoshiro256**
  278943 k ints/sec, 261021 k doubles/sec, mt19937
  272474 k ints/sec, 271866 k doubles/sec, xoroshiro128++
  271239 k ints/sec,  55285 k doubles/sec, ranf
  269604 k ints/sec, 247058 k doubles/sec, xoshiro128++
  265568 k ints/sec, 263181 k doubles/sec, gfsr4
  264719 k ints/sec, 251332 k doubles/sec, taus
  262505 k ints/sec,  52593 k doubles/sec, rand48
  261220 k ints/sec, 245894 k doubles/sec, tt800
  260970 k ints/sec, 254503 k doubles/sec, xoshiro128**
  247667 k ints/sec, 257275 k doubles/sec, random128-libc5
  233522 k ints/sec, 225526 k doubles/sec, taus113
  198365 k ints/sec, 199142 k doubles/sec, knuthran
  177301 k ints/sec, 175731 k doubles/sec, slatec
  155473 k ints/sec, 155889 k doubles/sec, fishman20
  152459 k ints/sec, 153368 k doubles/sec, minstd
  150635 k ints/sec, 150672 k doubles/sec, ran1
  148195 k ints/sec, 139760 k doubles/sec, fishman2x
  146463 k ints/sec, 151964 k doubles/sec, ran0
  137318 k ints/sec, 140146 k doubles/sec, lecuyer21
  119559 k ints/sec, 118447 k doubles/sec, ran2
  106604 k ints/sec,  97263 k doubles/sec, zuf
   85861 k ints/sec,  82312 k doubles/sec, mrg
   69997 k ints/sec,  71964 k doubles/sec, fishman18
   53427 k ints/sec,  53351 k doubles/sec, cmrg
   37676 k ints/sec,  40992 k doubles/sec, knuthran2
   35773 k ints/sec,  35799 k doubles/sec, ranlxs0
   23627 k ints/sec,  23611 k doubles/sec, ranlxs1
   20209 k ints/sec,  20051 k doubles/sec, ranlux
   13674 k ints/sec,  13634 k doubles/sec, ranlxs2
   13671 k ints/sec,  13614 k doubles/sec, ranlxd1
   11852 k ints/sec,  11536 k doubles/sec, ranlux389
    7430 k ints/sec,   7427 k doubles/sec, ranlxd2

Examples
========

The following program demonstrates the use of a random number generator
to produce uniform random numbers in the range [0.0, 1.0),

.. include:: examples/rngunif.c
   :code:

Here is the output of the program,

.. include:: examples/rngunif.txt
   :code:

The numbers depend on the seed used by the generator.  The default seed
can be changed with the :macro:`GSL_RNG_SEED` environment variable to
produce a different stream of numbers.  The generator itself can be
changed using the environment variable :macro:`GSL_RNG_TYPE`.  Here is the
output of the program using a seed value of 123 and the
multiple-recursive generator :code:`mrg`::

  $ GSL_RNG_SEED=123 GSL_RNG_TYPE=mrg ./a.out 

.. include:: examples/rngunif2.txt
   :code:

References and Further Reading
==============================

The subject of random number generation and testing is reviewed
extensively in Knuth's *Seminumerical Algorithms*.

* Donald E. Knuth, The Art of Computer Programming: Seminumerical
  Algorithms (Vol 2, 3rd Ed, 1997), Addison-Wesley, ISBN 0201896842.

Further information is available in the review paper written by Pierre
L'Ecuyer,

* P. L'Ecuyer, "Random Number Generation", Chapter 4 of the
  Handbook on Simulation, Jerry Banks Ed., Wiley, 1998, 93--137.

* http://www.iro.umontreal.ca/~lecuyer/papers.html in the file :file:`handsim.ps`.

The source code for the DIEHARD random number generator tests is also
available online,

* DIEHARD source code, G. Marsaglia, http://stat.fsu.edu/pub/diehard/

The source code for the DIEHARDER random number generator tests is
available online,
  
* DIEHARDER source code, Robert G. Brown, https://webhome.phy.duke.edu/~rgb/General/dieharder.php
  
A comprehensive set of random number generator tests is available from
NIST,

* NIST Special Publication 800-22, "A Statistical Test Suite for the
  Validation of Random Number Generators and Pseudo Random Number
  Generators for Cryptographic Applications".

* http://csrc.nist.gov/rng/

Acknowledgements
================

Thanks to Makoto Matsumoto, Takuji Nishimura and Yoshiharu Kurita for
making the source code to their generators (MT19937, MM&TN; TT800,
MM&YK) available under the GNU General Public License.  Thanks to Martin
Luscher for providing notes and source code for the RANLXS and
RANLXD generators.

.. lcg
.. [ LCG(n) := n * 69069 mod (2^32) ]
.. First 6: [69069, 475559465, 2801775573, 1790562961, 3104832285, 4238970681]
.. %2^31-1   69069, 475559465, 654291926, 1790562961, 957348638, 2091487034
.. mrg
.. [q([x1, x2, x3, x4, x5]) := [107374182 mod 2147483647 * x1 + 104480 mod 2147483647 * x5, x1, x2, x3, x4]]
..
.. cmrg
.. [q1([x1,x2,x3]) := [63308 mod 2147483647 * x2 -183326 mod 2147483647 * x3, x1, x2],
..  q2([x1,x2,x3]) := [86098 mod 2145483479 * x1 -539608 mod 2145483479 * x3, x1, x2] ]
..  initial for q1 is [69069, 475559465, 654291926]
..  initial for q2 is  [1790562961, 959348806, 2093487202]

.. tausworthe
..    [ b1(x) := rsh(xor(lsh(x, 13), x), 19),
..      q1(x) := xor(lsh(and(x, 4294967294), 12), b1(x)),
..      b2(x) := rsh(xor(lsh(x, 2), x), 25),
..      q2(x) := xor(lsh(and(x, 4294967288), 4), b2(x)),
..      b3(x) := rsh(xor(lsh(x, 3), x), 11),
..      q3(x) := xor(lsh(and(x, 4294967280), 17), b3(x)) ]
..      [s1, s2, s3] = [600098857, 1131373026, 1223067536] 
.. [2948905028, 441213979, 394017882]
