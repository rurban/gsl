.. index::
   single: associated Legendre functions
   single: Legendre functions, associated
   single: spherical harmonics

The associated Legendre functions are described in
Abramowitz & Stegun, Chapter 8.  These functions are declared in 
the header file :file:`gsl_sf_alf.h`.

The associated Legendre functions :math:`P_l^m(x)` are solutions
of the differential equation

.. only:: not texinfo

   .. math:: (1 - x^2) {d^2 \over dx^2} P_l^m(x) - 2x {d \over dx} P_l^m(x) +
             \left( l(l+1) - {m^2 \over 1 - x^2} \right) P_l^m(x) = 0

.. only:: texinfo

   ::

      (1 - x^2) d^2 P_l^m(x) / dx^2 P_l^m(x) - 2x d/dx P_l^m(x) +
      ( l(l+1) - m^2 / (1 - x^2) ) P_l^m(x) = 0

where the degree :math:`l` and order :math:`m` satisfy :math:`0 \le l` and
:math:`0 \le m \le l`.
The functions :math:`P_l^m(x)` grow combinatorially with
:math:`l` and can overflow in double precision for :math:`l` larger than about 150.
Alternatively, one may calculate normalized associated Legendre
polynomials. There are a number of different normalization conventions,
and these functions can be stably computed up to degree and order 2700. The
following normalizations are provided:

* **Schmidt semi-normalization**

  Schmidt semi-normalized associated Legendre polynomials are often
  used in the geophysical community and are defined as (Winch et al, 2005)

  .. only:: not texinfo

     .. math::

        S_l^0(x) &= P_l^0(x) \\
        S_l^m(x) &= (-1)^m \sqrt{2 {(l-m)! \over (l+m)!}} P_l^m(x), m > 0 

  .. only:: texinfo

     ::

        S_l^0(x) = P_l^0(x)
        S_l^m(x) = (-1)^m \sqrt((2(l-m)! / (l+m)!)) P_l^m(x), m > 0 

  The factor of :math:`(-1)^m` is called the Condon-Shortley phase
  factor and can be included if desired by setting the flag
  :macro:`GSL_SF_ALF_FLG_CSPHASE` in the function :func:`gsl_sf_alf_precompute`
  below. These functions satisfy the normalization condition,

  .. only:: not texinfo

     .. math::

        \int_{-1}^1 S_k^m(x) S_l^m(x) dx =
        \left\{
          \begin{array}{ll}
            \frac{2}{2l+1} \delta_{kl}, & m = 0 \\
            \frac{4}{2l+1} \delta_{kl}, & m > 0
          \end{array}
        \right.

  .. only:: texinfo

     ::

        \int_{-1}^1 S_k^m(x) S_l^m(x) dx = { 2/(2l+1) \delta_{kl}, m = 0
                                           { 4/(2l+1) \delta_{kl}, m > 0

* **Spherical Harmonic Normalization**

  The associated Legendre polynomials suitable for calculating spherical
  harmonics are defined as

  .. only:: not texinfo

     .. math:: Y_l^m(x) = (-1)^m \sqrt{{2l + 1 \over 4 \pi} {(l-m)! \over (l+m)!}} P_l^m(x)

  .. only:: texinfo

     ::

        Y_l^m(x) = (-1)^m \sqrt((2l + 1) * (l-m)! / (4 \pi) / (l+m)!) P_l^m(x)

  where again the phase factor :math:`(-1)^m` can be included or excluded
  if desired. These functions satisfy the normalization condition,

  .. only:: not texinfo

     .. math::

        \int_{-1}^1 Y_k^m(x) Y_l^m(x) dx = \frac{\delta_{kl}}{2\pi}

  .. only:: texinfo

     ::

        \int_{-1}^1 Y_k^m(x) Y_l^m(x) dx = \delta_{kl} / (2 \pi)

  Note that these functions, when coupled with the factor
  :math:`e^{i m \phi}` produce the orthonormalized complex spherical
  harmonics.

* **Full Normalization**

  The fully normalized associated Legendre polynomials are defined as

  .. only:: not texinfo

     .. math:: N_l^m(x) = (-1)^m \sqrt{(l + {1 \over 2}) {(l-m)! \over (l+m)!}} P_l^m(x)

  .. only:: texinfo
  
     ::
     
        N_l^m(x) = (-1)^m \sqrt((l + 1/2) (l-m)! / (l+m)!) P_l^m(x)

  and satisfy the normalization condition,

  .. math:: \int_{-1}^1 N_k^m(x) N_l^m(x) dx = \delta_{kl}

* :math:`4 \pi` **Normalization**

  The :math:`4 \pi` normalized associated Legendre polynomials are often used in geodesy and are defined as

  .. only:: not texinfo

     .. math:: R_l^m(x) = (-1)^m \sqrt{(2 - \delta_{m0}) (2 l + 1) {(l-m)! \over (l+m)!}} P_l^m(x)

  .. only:: texinfo
  
     ::
     
        R_l^m(x) = (-1)^m \sqrt((2 - \delta_{m0}) (2 l + 1) (l-m)! / (l+m)!) P_l^m(x)

  These functions satisfy the normalization condition,

  .. only:: not texinfo

     .. math:: \int_{-1}^1 R_k^m(x) R_l^m(x) dx = 2 \left( 2 - \delta_{m0} \right) \delta_{kl}

  .. only:: texinfo

     ::

        \int_{-1}^1 R_k^m(x) R_l^m(x) dx = 2 (2 - \delta_{m0}) \delta_{kl}

  When used in the definition of real spherical harmonics, they satisfy a
  :math:`4\pi` normalization condition when integrated over the unit sphere.
  More information on these functions can be found in Hofmann-Wellenhof and Moritz, 2006.

The normalized associated Legendre routines below use a recurrence
relation which is stable up to a degree and order of about 2700.
Beyond this, the computed functions could suffer from underflow
leading to incorrect results. Routines are provided to compute
first and second derivatives
:math:`dP_l^m(x)/dx` and :math:`d^2 P_l^m(x)/dx^2` as well as their alternate
versions :math:`d P_l^m(\cos{\theta})/d\theta` and
:math:`d^2 P_l^m(\cos{\theta})/d\theta^2`. While there is a simple
scaling relationship between the two forms, the derivatives
involving :math:`\theta` are heavily used in spherical harmonic
expansions, and also do not suffer from singularities at the poles,
:math:`x = \pm 1`, and so these routines are also provided.

In the functions below, a parameter of type :type:`gsl_sf_alf_t`
specifies the type of normalization to use. The possible values are

.. type:: gsl_sf_alf_t

   ============================= ===============================================================================
   Value                         Description
   ============================= ===============================================================================
   :code:`GSL_SF_ALF_NONE`       The unnormalized associated Legendre polynomials :math:`P_l^m(x)`
   :code:`GSL_SF_ALF_SCHMIDT`    The Schmidt semi-normalized associated Legendre polynomials :math:`S_l^m(x)`
   :code:`GSL_SF_ALF_SPHARM`     The spherical harmonic associated Legendre polynomials :math:`Y_l^m(x)`
   :code:`GSL_SF_ALF_FULL`       The fully normalized associated Legendre polynomials :math:`N_l^m(x)`
   :code:`GSL_SF_ALF_FOURPI`     The :math:`4\pi` normalized associated Legendre polynomials :math:`R_l^m(x)`
   ============================= ===============================================================================

The routines below which return an array of ALFs organize the outputs into blocks
of the same order :math:`m`.  This scheme uses the following indexing function to
locate an ALF of degree :math:`l` and order :math:`m`,

.. math:: \mathcal{I}_m(l,m,L) = m L - \frac{m(m-1)}{2} + l

where :math:`L \geq 0` is the maximum degree, :code:`lmax`. This corresponds to the following memory layout,

.. math:: l \quad \overbrace{0 \; 1 \; 2 \; \cdots \; L}^{m = 0} \quad \overbrace{1 \; 2 \; \cdots \; L}^{m = 1} \quad \overbrace{2 \; 3 \; \cdots \; L}^{m = 2} \quad \cdots \quad \overbrace{M \; \cdots \; L}^{m = M}

Here, :math:`M` is the maximum order, :code:`mmax`, which satisfies :math:`0 \leq M \leq L`.
This memory layout corresponds to the order that ALFs are computed in their
recurrence relations, and is therefore cache-efficient. The following
code demonstrates how to access the array elements in order,

.. code::

   idx = 0;
   for (m = 0; m <= mmax; ++m) {
     for (l = m; l <= lmax; ++l) {
       double value = Plm[idx]; /* (l,m) element */
       ++idx;
     }
   }

.. function:: int gsl_sf_alf_precompute (const gsl_sf_alf_t norm, const size_t lmax, const size_t mmax, const size_t flags, double result_array[])

   This function precomputes the multiplicative factors needed for the
   associated Legendre recurrence relations. The input :data:`norm`
   specifies the ALF normalization. The input :data:`lmax` specifies
   the maximum ALF degree. The input :data:`mmax` specifies the maximum
   ALF order. The input :data:`flags` is a bitmask which
   specifies how the ALFs are computed and stored. It can contain the
   following values,

   .. macro:: GSL_SF_ALF_FLG_CSPHASE

      This flag will include the Condon-Shortley phase factor in the calculation
      of the ALFs

   The output array :data:`result_array` should have a length
   as returned by the function :func:`gsl_sf_alf_array_size`. The computed
   recurrence factors are stored at the end of :data:`result_array`, leaving
   room at the front for the calculation of the ALFs.

   This routine must be called prior to calling ALF array functions.

.. function:: int gsl_sf_alf_array (const size_t lmax, const size_t mmax, const double x, double result_array[])

   This function calculates all associated Legendre polynomials
   for :math:`0 \le l \le lmax` and :math:`0 \le m \le \min{(l,mmax)}` for :math:`|x| \le 1`.
   The :data:`norm` parameter specifies which normalization is used.
   The normalized :math:`P_l^m(x)` values are stored in :data:`result_array`, whose
   minimum size can be obtained from calling :func:`gsl_sf_alf_array_size`.
   The array index of :math:`P_l^m(x)` is obtained from calling
   :code:`gsl_sf_alf_array_index(l, m, lmax)`.

   The function :func:`gsl_sf_alf_precompute` must be called first
   using the same :data:`lmax` and :data:`mmax` inputs
   to initialize :data:`result_array` with the multiplicative factors
   used in the recurrence relations.

.. function:: int gsl_sf_alf_deriv_array (const size_t lmax, const size_t mmax, const double x, double result_array[], double result_deriv_array[])

   This function calculates all associated Legendre
   functions and their first derivatives for :math:`0 \leq l \leq lmax`
   and :math:`0 \leq m \leq \min{(l,mmax)}` for :math:`|x| < 1`.
   The :math:`P_l^m(x)` values and their derivatives
   :math:`dP_l^m(x)/dx` are stored in :data:`result_array` and
   :data:`result_deriv_array` respectively.

   Note that for some orders :math:`m`, the derivatives :math:`dP_l^m(x)/dx` have
   singularities at the end points :math:`x = \pm 1`, and so this function only
   accepts interior points as input, :math:`x \in (-1,1)`.

   The function :func:`gsl_sf_alf_precompute` must be called first
   using the same :data:`lmax` and :data:`mmax` inputs
   to initialize :data:`result_array` with the multiplicative factors
   used in the recurrence relations.

.. function:: int gsl_sf_alf_vsh_array (const size_t lmax, const size_t mmax, const double x, double result_array[], double result_deriv_array[])

   In vector spherical harmonic expansions, it is often necessary to simultaneously
   compute terms of the form

   .. math:: \frac{d}{d\theta} P_l^m(\cos{\theta}) e^{i m \phi}

   and

   .. math:: \frac{i m}{\sin{\theta}} P_l^m(\cos{\theta}) e^{i m \phi}

   Analogous expressions to those above also arise when using real-valued spherical
   harmonics. The factor :math:`\frac{1}{\sin{\theta}}` in the second term could
   present problems at the poles, :math:`\theta = 0,\pi` (e.g. :math:`x = \pm 1`). However,
   for :math:`m \neq 0`, :math:`P_l^m(\cos{\theta}) / \sin{\theta}` is well defined
   and can be computed stably for all :math:`\theta`.

   This function computes and stores the following quantities in the output arrays,

   .. math::

      \textrm{result\_array[index(l,0,lmax)]} &= P_l^0(x) \\
      \textrm{result\_array[index(l,m,lmax)]} &= \frac{P_l^m(x)}{\sin{\theta}}, \qquad m > 0 \\
      \textrm{result\_deriv\_array[index(l,m,lmax)]} &= \frac{d}{d\theta} P_l^m(x), \qquad m \geq 0
   
   The :code:`index` notation above refers to the :func:`gsl_sf_alf_array_index` function.
   All associated Legendre functions and their first :math:`\theta` derivatives are computed
   for :math:`0 \leq l \leq lmax` and :math:`0 \leq m \leq \min{(l,mmax)}` for :math:`|x| \leq 1`.

   The function :func:`gsl_sf_alf_precompute` must be called first
   using the same :data:`lmax` and :data:`mmax` inputs
   to initialize :data:`result_array` with the multiplicative factors
   used in the recurrence relations.

.. function:: size_t gsl_sf_alf_nlm(const size_t lmax, const size_t mmax)

   This function returns the total number of associated Legendre
   functions :math:`P_l^m(x)` for a given :data:`lmax` and :data:`mmax`.

   An inline version of this function is used if :macro:`HAVE_INLINE` is
   defined.

.. function:: size_t gsl_sf_alf_array_size (const size_t lmax, const size_t mmax)

   This function returns the minimum array size for maximum degree :data:`lmax`
   and maximum order :data:`mmax` needed for the array versions of the associated Legendre functions.
   Size is calculated as the total number of :math:`P_l^m(x)` functions
   (see :func:`gsl_sf_alf_nlm`),
   plus extra space for precomputing multiplicative factors used in the
   recurrence relations.

.. function:: size_t gsl_sf_alf_array_index (const size_t l, const size_t m, const size_t lmax)

   This function returns the index into the ALF arrays corresponding
   to the function of degree :math:`l` and order :math:`m`. The maximum
   degree :math:`L` must also be specified in the parameter :data:`lmax`.

   An inline version of this function is used if :macro:`HAVE_INLINE` is
   defined.
