/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <config.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include "gsl_sf_clausen.h"
#include "gsl_sf_trig.h"
#include "gsl_sf_log.h"
#include "gsl_sf_dilog.h"


/* Evaluate series for real dilog(x)
 * Sum[ x^k / k^2, {k,1,Infinity}]
 *
 * Converges rapidly for |x| < 1/2.
 */
static
int
dilog_series(const double x, gsl_sf_result * result)
{
  const int kmax = 1000;
  double sum  = x;
  double term = x;
  int k;
  for(k=2; k<kmax; k++) {
    double rk = (k-1.0)/k;
    term *= x;
    term *= rk*rk;
    sum += term;
    if(fabs(term/sum) < GSL_DBL_EPSILON) break;
  }

  result->val = sum;
  result->err = 2.0 * fabs(term);

  if(k == kmax)
    return GSL_EMAXITER;
  else
    return GSL_SUCCESS;
}


/* Assumes x >= 0.0
 */
static
int
dilog_xge0(const double x, gsl_sf_result * result)
{
  if(x > 2.0) {
    const double log_x = log(x);
    gsl_sf_result ser;
    int stat_ser = dilog_series(1.0/x, &ser);
    result->val  = M_PI*M_PI/3.0 - ser.val - 0.5*log_x*log_x;
    result->err  = GSL_DBL_EPSILON * fabs(log_x) + ser.err;
    result->val += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat_ser;
  }
  else if(x > 1.01) {
    const double log_x    = log(x);
    const double log_term = log_x * (log(1.0-1.0/x) + 0.5*log_x);
    gsl_sf_result ser;
    int stat_ser = dilog_series(1.0 - 1.0/x, &ser);
    result->val  = M_PI*M_PI/6.0 + ser.val - log_term;
    result->err  = GSL_DBL_EPSILON * fabs(log_x) + ser.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat_ser;
  }
  else if(x > 1.0) {
    /* series around x = 1.0 */
    const double eps = x - 1.0;
    const double lne = log(eps);
    const double c0 = M_PI*M_PI/6.0;
    const double c1 =   1.0 - lne;
    const double c2 = -(1.0 - 2.0*lne)/4.0;
    const double c3 =  (1.0 - 3.0*lne)/9.0;
    const double c4 = -(1.0 - 4.0*lne)/16.0;
    const double c5 =  (1.0 - 5.0*lne)/25.0;
    const double c6 = -(1.0 - 6.0*lne)/36.0;
    const double c7 =  (1.0 - 7.0*lne)/49.0;
    const double c8 = -(1.0 - 8.0*lne)/64.0;
    result->val = c0+eps*(c1+eps*(c2+eps*(c3+eps*(c4+eps*(c5+eps*(c6+eps*(c7+eps*c8)))))));
    result->err = 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x == 1.0) {
    result->val = M_PI*M_PI/6.0;
    result->err = 2.0 * GSL_DBL_EPSILON * M_PI*M_PI/6.0;
    return GSL_SUCCESS;
  }
  else if(x > 0.5) {
    const double log_x = log(x);
    gsl_sf_result ser;
    int stat_ser = dilog_series(1.0-x, &ser);
    result->val  = M_PI*M_PI/6.0 - ser.val - log_x*log(1.0-x);
    result->err  = GSL_DBL_EPSILON * fabs(log_x) + ser.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat_ser;
  }
  else if(x > 0.0) {
    return dilog_series(x, result);
  }
  else {
    /* x == 0.0 */
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
}


/* Evaluate the series representation for Li2(z):
 *
 *   Li2(z) = Sum[ |z|^k / k^2 Exp[i k arg(z)], {k,1,Infinity}]
 *   |z|    = r
 *   arg(z) = theta
 *   
 * Assumes 0 < r < 1. 
 */
static
int
dilogc_series_1(double r, double cos_theta, double sin_theta,
                double * real_result, double * imag_result)
{
  double alpha = 1.0 - cos_theta;
  double beta  = sin_theta;
  double ck = cos_theta;
  double sk = sin_theta;
  double rk = r;
  double real_sum = r*ck;
  double imag_sum = r*sk;
  int kmax = 50 + (int)(22.0/(-log(r))); /* tuned for double-precision */
  int k;
  for(k=2; k<kmax; k++) {
    double ck_tmp = ck;
    ck = ck - (alpha*ck + beta*sk);
    sk = sk - (alpha*sk - beta*ck_tmp);
    rk *= r;
    real_sum += rk/((double)k*k) * ck;
    imag_sum += rk/((double)k*k) * sk;
  }
  
  *real_result = real_sum;
  *imag_result = imag_sum;
  return GSL_SUCCESS;
}


/* Evaluate a series for Li_2(z) when |z| is near 1.
 * This is uniformly good away from z=1.
 *
 *   Li_2(z) = Sum[ a^n/n! H_n(theta), {n, 0, Infinity}]
 *
 * where
 *   H_n(theta) = Li_2(Exp(I theta))
 *   a = ln(r)
 *
 *  H_0(t) = Gl_2(t) + i Cl_2(t)
 *  H_1(t) = 1/2 ln(2(1-c)) + I atan2(-s, 1-c)
 *  H_2(t) = -1/2 + I/2 s/(1-c)
 *  H_3(t) = -1/2 /(1-c)
 *  H_4(t) = -I/2 s/(1-c)^2
 *  H_5(t) = 1/2 (2 + c)/(1-c)^2
 *  H_6(t) = I/2 s/(1-c)^5 (8(1-c) - s^2 (3 + c))
 *
 *  assumes: 0 <= theta <= 2Pi
 */
static
int
dilogc_series_2(double r, double theta, double cos_theta, double sin_theta,
                double * real_result, double * imag_result)
{
  double a = log(r);
  double omc = 1.0 - cos_theta;
  double H_re[7];
  double H_im[7];
  double an, nfact;
  double sum_re, sum_im;
  gsl_sf_result Him0;
  int n;

  H_re[0] = M_PI*M_PI/6.0 + 0.25*(theta*theta - 2.0*M_PI*fabs(theta));
  gsl_sf_clausen_impl(theta, &Him0);
  H_im[0] = Him0.val;
  
  H_re[1] = 0.5*log(2.0*omc);
  H_im[1] = atan2(-sin_theta, omc);
  
  H_re[2] = -0.5;
  H_im[2] = 0.5 * sin_theta/omc;
  
  H_re[3] = -0.5/omc;
  H_im[3] = 0.0;
  
  H_re[4] = 0.0;
  H_im[4] = -0.5*sin_theta/(omc*omc);
  
  H_re[5] = 0.5 * (2.0 + cos_theta)/(omc*omc);
  H_im[5] = 0.0;
  
  H_re[6] = 0.0;
  H_im[6] = 0.5 * sin_theta/(omc*omc*omc*omc*omc)
            * (8*omc - sin_theta*sin_theta*(3 + cos_theta));
 
  sum_re = H_re[0];
  sum_im = H_im[0];
  an = 1.0;
  nfact = 1.0;
  for(n=1; n<=6; n++) {
    double t;
    an *= -a;
    nfact *= n;
    t = an/nfact;
    sum_re += t * H_re[n];
    sum_im += t * H_im[n];
  }
  
  *real_result = sum_re;
  *imag_result = sum_im;
  return GSL_SUCCESS;
}


/* complex dilogarithm in the unit disk
 * assumes:  r < 1  and  0 <= theta <= 2Pi
 */
static
int
dilogc_unitdisk(double r, double theta, double * real_dl, double * imag_dl)
{
  const double zeta2 = M_PI*M_PI/6.0;
  int stat_dilog;
  double cos_theta = cos(theta);
  double sin_theta = sin(theta);
  double x = r * cos_theta;
  double y = r * sin_theta;
  double x_tmp, y_tmp, r_tmp;
  double result_re_tmp, result_im_tmp;

  /* Reflect away from z = 1 if
   * we are too close.
   */
  if(x > 0.5) {
    x_tmp = 1.0 - x;
    y_tmp = -y;
    r_tmp = sqrt(x_tmp*x_tmp + y_tmp*y_tmp);
  }
  else {
    x_tmp = x;
    y_tmp = y;
    r_tmp = r;
  }

  /* Calculate dilog of the transformed variable.
   */
  if(r_tmp < 0.98) {
    double cos_theta_tmp = x_tmp / r_tmp;
    double sin_theta_tmp = y_tmp / r_tmp;
    stat_dilog = dilogc_series_1(r_tmp, cos_theta_tmp, sin_theta_tmp,
				 &result_re_tmp, &result_im_tmp
				 );
  }
  else {
    double cos_theta_tmp = x_tmp / r_tmp;
    double sin_theta_tmp = y_tmp / r_tmp;
    double theta_tmp = atan2(y_tmp, x_tmp);
    stat_dilog = dilogc_series_2(r_tmp, theta_tmp, cos_theta_tmp, sin_theta_tmp,
                                 &result_re_tmp, &result_im_tmp
				 );
  }

  /* Unwind reflection if necessary.
   *
   * Li2(z) = -Li2(1-z) + zeta(2) - ln(z) ln(1-z)
   */
  if(x > 0.5) {
    double lnz    =  log(r);                 /*  log(|z|)   */
    double lnomz  =  log(r_tmp);             /*  log(|1-z|) */
    double argz   =  theta;                  /*  arg(z)     */
    double argomz =  atan2(y_tmp, x_tmp);    /*  arg(1-z)   */
    *real_dl = -result_re_tmp + zeta2 - lnz*lnomz + argz*argomz;
    *imag_dl = -result_im_tmp - argz*lnomz - argomz*lnz;
  }
  else {
    *real_dl = result_re_tmp;
    *imag_dl = result_im_tmp;
  }

  return stat_dilog;
}



/*-*-*-*-*-*-*-*-*-*-*-* (semi)Private Implementations *-*-*-*-*-*-*-*-*-*-*-*/


int
gsl_sf_dilog_impl(const double x, gsl_sf_result * result)
{
  if(x >= 0.0) {
    return dilog_xge0(x, result);
  }
  else {
    gsl_sf_result d1, d2;
    int stat_d1 = dilog_xge0( -x, &d1);
    int stat_d2 = dilog_xge0(x*x, &d2);
    result->val = -d1.val + 0.5 * d2.val;
    result->err =  d1.err + 0.5 * d2.err;
    return GSL_ERROR_SELECT_2(stat_d1, stat_d2);
  }
}


int
gsl_sf_complex_dilog_impl(const double r, double theta,
                          gsl_sf_result * real_dl, gsl_sf_result * imag_dl)
{
  if(theta < 0.0 || theta > 2.0*M_PI) {
    gsl_sf_angle_restrict_pos_impl(&theta);
  }

  if(r == 0.0) {
    real_dl->val = 0.0;
    real_dl->err = 0.0;
    imag_dl->val = 0.0;
    imag_dl->err = 0.0;
    return GSL_SUCCESS;
  }

  /* Trap cases of real-valued argument.
   */
  if(theta == 0.0) {
    imag_dl->val = ( r > 1.0 ? -M_PI * log(r) : 0.0 );
    imag_dl->err = GSL_DBL_EPSILON * fabs(imag_dl->val);
    return gsl_sf_dilog_impl(r, real_dl);
  }
  if(theta == M_PI) {
    imag_dl->val = 0.0;
    imag_dl->err = 0.0;
    return gsl_sf_dilog_impl(-r, real_dl);
  }

  /* Trap unit circle case.
   */
  if(r == 1.0) {
    real_dl->val = M_PI*M_PI/6.0 + 0.25*(theta*theta - 2.0*M_PI*fabs(theta));
    real_dl->err = 3.0 * GSL_DBL_EPSILON * fabs(real_dl->val);
    return gsl_sf_clausen_impl(theta, imag_dl);
  }

  /* Generic case.
   */
  {
    int stat_dilog;
    double r_tmp, theta_tmp;
    double result_re_tmp, result_im_tmp;

    /* Reduce argument to unit disk.
     */
    if(r > 1.0) {
      r_tmp     = 1.0 / r;
      theta_tmp = 2.0*M_PI - theta;
    }
    else {
      r_tmp     = r;
      theta_tmp = theta;
    }

    /* Calculate in the unit disk.
     */
    stat_dilog = dilogc_unitdisk(r_tmp, theta_tmp,
                                 &result_re_tmp, &result_im_tmp
				 );

    /* Unwind the inversion if necessary. We calculate
     * the imaginary part explicitly if using the inversion
     * because there is no simple relationship between
     * arg(1-z) and arg(1 - 1/z), which is the "omega"
     * term in [Lewin A.2.5 (1)].
     */
    if(r > 1.0) {
      const double zeta2 = M_PI*M_PI/6.0;
      double x = r * cos(theta);
      double y = r * sin(theta);
      double omega = atan2(y, 1.0-x);
      double lnr = log(r);
      double pmt = M_PI - theta;
      gsl_sf_result Cl_a, Cl_b, Cl_c;
      gsl_sf_clausen_impl(2.0*omega, &Cl_a);
      gsl_sf_clausen_impl(2.0*theta, &Cl_b);
      gsl_sf_clausen_impl(2.0*(omega+theta), &Cl_c);
      real_dl->val = -result_re_tmp - 0.5*lnr*lnr + 0.5*pmt*pmt - zeta2;
      real_dl->err = 4.0 * GSL_DBL_EPSILON * fabs(real_dl->val);
      imag_dl->val = omega*log(r) + 0.5*(Cl_a.val + Cl_b.val - Cl_c.val);
      imag_dl->err = GSL_DBL_EPSILON * fabs(imag_dl->val) + 0.5*(Cl_a.err + Cl_b.err + Cl_c.err);
    }
    else {
      real_dl->val = result_re_tmp;
      real_dl->err = GSL_DBL_EPSILON * fabs(real_dl->val);
      imag_dl->val = result_im_tmp;
      imag_dl->err = GSL_DBL_EPSILON * fabs(imag_dl->val);
    }

    return stat_dilog;
  }
}


/*-*-*-*-*-*-*-*-*-*-*-* Functions w/ Error Handling *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_dilog_e(const double x, gsl_sf_result * result)
{
  int status = gsl_sf_dilog_impl(x, result);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_dilog_e", status);
  }
  return status;
}


int
gsl_sf_complex_dilog_e(const double x, const double y,
                       gsl_sf_result * result_re,
                       gsl_sf_result * result_im)
{
  int status = gsl_sf_complex_dilog_impl(x, y, result_re, result_im);
  if(status != GSL_SUCCESS) {
    GSL_ERROR("gsl_sf_complex_dilog_e", status);
  }
  return status;
}
