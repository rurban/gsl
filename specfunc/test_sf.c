/* Author:  G. Jungman
 * RCS:     $Id$
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl_math.h>
#include <gsl_errno.h>
#include <gsl_test.h>
#include <gsl_sf.h>


double frac_diff(double x1, double x2)
{
  return fabs((x1-x2)/(x1+x2));
}


int check_airy(void)
{
  int status = 0;
  int s;

  /** functions */

  s = 0;
  s += ( frac_diff(gsl_sf_airy_Ai(-5),                   0.3507610090241142    ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai(-0.3000000000000094),  0.4309030952855831    ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai(0.6999999999999907),   0.1891624003981519    ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai(1.649999999999991),    0.05831058618720882   ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai(2.54999999999999),     0.01446149513295428   ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai(3.499999999999987),    0.002584098786989702  ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai(5.39999999999998),     4.272986169411866e-05 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_airy_Ai");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_airy_Ai_scaled(-5),            gsl_sf_airy_Ai(-5)         ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_scaled(-0.5),          gsl_sf_airy_Ai(-0.5)       ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_scaled(0.6999999999999907),   0.2795125667681217  ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_scaled(1.649999999999991),    0.2395493001442741  ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_scaled(2.54999999999999),     0.2183658595899388  ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_scaled(3.499999999999987),    0.2032920808163519  ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_scaled(5.39999999999998),     0.1836050093282229  ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_airy_Ai_scaled");
  status += s;
  
  s = 0;
  s += ( frac_diff(gsl_sf_airy_Bi(-5),                   -0.1383691349016005    ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi(-0.3000000000000094),   0.4779778401098885    ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi(0.6999999999999907),    0.9733286558781599    ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi(1.649999999999991),     2.196407956850028     ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi(2.54999999999999),      6.973628612493443     ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi(3.499999999999987),     33.05550675461069     ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi(5.39999999999998),      1604.476078241272     ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_airy_Bi");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_airy_Bi_scaled(-5),            gsl_sf_airy_Bi(-5)        ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_scaled(-0.5),          gsl_sf_airy_Bi(-0.5)      ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_scaled(0.6999999999999907),  0.6587080754582302  ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_scaled(1.649999999999991),   0.5346449995597539  ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_scaled(2.54999999999999),    0.461835455542297   ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_scaled(3.499999999999987),   0.4201771882353061  ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_scaled(5.39999999999998),    0.3734050675720473  ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_airy_Bi_scaled");
  status += s;
  
  /** derivatives */

  s = 0;
  s += ( frac_diff(gsl_sf_airy_Ai_deriv(-5),                   0.3271928185544435)  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_deriv(-0.5500000000000094), -0.1914604987143629)  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_deriv(0.4999999999999906),  -0.2249105326646850)  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_deriv(1.899999999999992),   -0.06043678178575718) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_deriv(3.249999999999988),   -0.007792687926790889) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_deriv(5.199999999999981),   -0.0001589434526459543) > 1.e-14 );
  gsl_test(s, "  gsl_sf_airy_Ai_deriv");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_airy_Ai_deriv_scaled(-5),          gsl_sf_airy_Ai_deriv(-5)   )  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_deriv_scaled(-0.5),        gsl_sf_airy_Ai_deriv(-0.5) )  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_deriv_scaled(0.5499999999999906), -0.2874057279170166 )  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_deriv_scaled(1.499999999999991),  -0.3314199796863637 )  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_deriv_scaled(2.49999999999999),   -0.366108938475162  )  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_deriv_scaled(3.649999999999986),  -0.3974033831453963 )  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Ai_deriv_scaled(6.299999999999977),  -0.4508799189585947 )  > 1.e-14 );
  gsl_test(s, "  gsl_sf_airy_Ai_deriv_scaled");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_airy_Bi_deriv(-5),                   0.778411773001899 )  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_deriv(-0.5500000000000094),  0.5155785358765014)  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_deriv(0.4999999999999906),   0.5445725641405883)  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_deriv(1.899999999999992),    3.495165862891568 ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_deriv(3.249999999999988),    36.55485149250338 ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_deriv(5.199999999999981),    2279.748293583233 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_airy_Bi_deriv");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_airy_Bi_deriv_scaled(-5),          gsl_sf_airy_Bi_deriv(-5)   )  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_deriv_scaled(-0.5),        gsl_sf_airy_Bi_deriv(-0.5) )  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_deriv_scaled(0.5499999999999906),  0.4322811281817566 )  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_deriv_scaled(1.499999999999991),   0.5542307563918037 )  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_deriv_scaled(2.49999999999999),    0.6755384441644985 )  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_deriv_scaled(3.649999999999986),   0.7613959373000228 )  > 1.e-14 );
  s += ( frac_diff(gsl_sf_airy_Bi_deriv_scaled(6.299999999999977),   0.8852064139737571 )  > 1.e-14 );
  gsl_test(s, "  gsl_sf_airy_Bi_deriv_scaled");
  status += s;
  return status;
}


int check_bessel(void)
{
  int status = 0;
  int s;

  s = 0;
  s += ( frac_diff(gsl_sf_bessel_Jnu(0.0001,10.0),  -0.2459270166445205       ) > 1.e-12 );
  s += ( frac_diff(gsl_sf_bessel_Jnu( 1.0, 0.001),  0.0004999999375000026     ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_bessel_Jnu( 1.0,   1.0),  0.4400505857449335160     ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_bessel_Jnu(30.0,   1.0),  3.482869794251482902e-42  ) > 1.e-13 );
  s += ( frac_diff(gsl_sf_bessel_Jnu(30.0, 100.0),  0.08146012958117222297    ) > 1.e-13 );
  s += ( frac_diff(gsl_sf_bessel_Jnu(10.0,   1.0),  2.6306151236874532070e-10 ) > 1.e-13 );
  s += ( frac_diff(gsl_sf_bessel_Jnu(10.0, 100.0),  -0.05473217693547201474   ) > 1.e-13 );
  gsl_test(s, "  gsl_sf_bessel_Jnu");
  status += s;


  return status;
}


int check_clausen(void)
{
  double y;
  int status = 0;
  int s;

  s = 0;
  s += ( frac_diff( gsl_sf_clausen(M_PI/20.0), 0.4478882448133546 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_clausen(M_PI/6.0),  0.8643791310538927 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_clausen(M_PI/3.0),  1.0149416064096535 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_clausen(  2.0*M_PI + M_PI/3.0),  1.0149416064096535 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_clausen(100.0*M_PI + M_PI/3.0),  1.0149416064096535 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_clausen");
  status += s;
  
  s = 0;
  s += ( gsl_sf_clausen_impl(1.0e+6*M_PI + M_PI/3.0, &y) != GSL_ELOSS);
  gsl_test(s, "  gsl_sf_clausen: trap accuracy loss from large argument");
  status += s;

  return status;
}


#define PRINT(n) printf("%20.16g  %20.16g  %20.16g  %20.16g    %20.16g\n", F[n], Fp[n], G[n], Gp[n], Fp[n]*G[n]-Gp[n]*F[n])

int check_coulomb(void)
{
  int status = 0;
  int s;
  
  const int kmax = 20;
  double F[kmax+1], Fp[kmax+1], G[kmax+1], Gp[kmax+1];
  double Fe, Ge;
  double lam_min;
  double eta, x;

  lam_min = 0.0;
  eta = -1000.0;
  x = 1.0;
  gsl_sf_coulomb_wave_FGp_impl(lam_min, kmax, eta, x, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( frac_diff(  F[0],  9.68222518991e-02 ) > 1.e-10 );
  s += ( frac_diff( Fp[0],  5.12063396274e+00 ) > 1.e-10 );
  s += ( frac_diff(  G[0],  1.13936784380e-01 ) > 1.e-10 );
  s += ( frac_diff( Gp[0], -4.30243486522e+00 ) > 1.e-10 );
  gsl_test(s, "  coulomb(0, -1000, 1)");
  PRINT(0);
  status += s;

  lam_min = 0.0;
  eta = -50.0;
  x = 5.0;
  gsl_sf_coulomb_wave_FGp_impl(lam_min, kmax, eta, x, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( frac_diff(  F[0],  1.52236975714e-01 ) > 1.e-10 );
  s += ( frac_diff( Fp[0],  2.03091041166e+00 ) > 1.e-10 );
  s += ( frac_diff(  G[0],  4.41680690236e-01 ) > 1.e-10 );
  s += ( frac_diff( Gp[0], -6.76485374767e-01 ) > 1.e-10 );
  gsl_test(s, "  coulomb(0, -50, 5)");
  PRINT(0);
  status += s;

  s = 0;
  s += ( frac_diff(  F[10], -3.68114360218e-01 ) > 1.e-10 );
  s += ( frac_diff( Fp[10],  1.33846751032e+00 ) > 1.e-10 );
  s += ( frac_diff(  G[10],  3.31588324611e-01 ) > 1.e-10 );
  s += ( frac_diff( Gp[10],  1.51088862814e+00 ) > 1.e-10 );
  gsl_test(s, "  coulomb(10, -50, 5)");
  PRINT(10);
  status += s;

  lam_min = 0.0;
  eta = -4.0;
  x = 5.0;
  gsl_sf_coulomb_wave_FGp_impl(lam_min, kmax, eta, x, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( frac_diff(  F[0],  4.07862723006e-01 ) > 1.e-10 );
  s += ( frac_diff( Fp[0],  1.09821233636e+00 ) > 1.e-10 );
  s += ( frac_diff(  G[0],  6.74327035383e-01 ) > 1.e-10 );
  s += ( frac_diff( Gp[0], -6.36110427280e-01 ) > 1.e-10 );
  gsl_test(s, "  coulomb(0, -4, 5)");
  PRINT(0);
  status += s;
  s = 0;
  s += ( frac_diff(  F[3], -2.56863093558e-01 ) > 1.e-10 );
  s += ( frac_diff( Fp[3],  1.14322942201e+00 ) > 1.e-10 );
  s += ( frac_diff(  G[3],  7.87989922393e-01 ) > 1.e-10 );
  s += ( frac_diff( Gp[3],  3.85990587811e-01 ) > 1.e-10 );
  gsl_test(s, "  coulomb(3, -4, 5)");
  PRINT(3);
  status += s;

  lam_min = 0.0;
  eta = 1.0;
  x = 2.0;
  gsl_sf_coulomb_wave_FGp_impl(lam_min, kmax, eta, x, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( frac_diff(  F[0],  6.61781613833e-01 ) > 1.e-10 );
  s += ( frac_diff( Fp[0],  4.81557455710e-01 ) > 1.e-10 );
  s += ( frac_diff(  G[0],  1.27577878477e+00 ) > 1.e-10 );
  s += ( frac_diff( Gp[0], -5.82728813097e-01 ) > 1.e-10 );
  gsl_test(s, "  coulomb(0, 1, 2)");
  PRINT(0);
  status += s;

  lam_min = 0.0;
  eta = 8.0;
  x = 1.05;
  gsl_sf_coulomb_wave_FGp_impl(lam_min, kmax, eta, x, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( frac_diff(  F[0],  9.88270778102e-09 ) > 1.e-6 );
  s += ( frac_diff( Fp[0],  4.00516771647e-08 ) > 1.e-6 );
  s += ( frac_diff(  G[0],  1.33312774446e+07 ) > 1.e-6 );
  s += ( frac_diff( Gp[0], -4.71591379533e+07 ) > 1.e-6 );
  gsl_test(s, "  coulomb(0, 8, 1.05)");
  PRINT(0);
  status += s;

  lam_min = 0.0;
  eta = 10.0;
  x = 5.0;
  gsl_sf_coulomb_wave_FGp_impl(lam_min, kmax, eta, x, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( frac_diff(  F[0],  1.72074542808e-06 ) > 1.e-10 );
  s += ( frac_diff( Fp[0],  3.09759950467e-06 ) > 1.e-10 );
  s += ( frac_diff(  G[0],  1.67637564193e+05 ) > 1.e-10 );
  s += ( frac_diff( Gp[0], -2.79370763593e+05 ) > 1.e-10 );
  gsl_test(s, "  coulomb(0, 10, 5)");
  PRINT(0);
  status += s;

  lam_min = 0.0;
  eta = 25.0;
  x = 10.0;
  gsl_sf_coulomb_wave_FGp_impl(lam_min, kmax, eta, x, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( frac_diff(  F[0],  1.54838713210e-16 ) > 1.e-10 );
  s += ( frac_diff( Fp[0],  3.14570932195e-16 ) > 1.e-10 );
  s += ( frac_diff(  G[0],  1.61423768840e+15 ) > 1.e-10 );
  s += ( frac_diff( Gp[0], -3.17884161766e+15 ) > 1.e-10 );
  gsl_test(s, "  coulomb(0, 25, 10)");
  PRINT(0);
  status += s;

  lam_min = 0.0;
  eta = 1.0;
  x = 9.2;
  gsl_sf_coulomb_wave_FGp_impl(lam_min, kmax, eta, x, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( frac_diff(  F[0], -2.56320123198e-01 ) > 1.e-10 );
  s += ( frac_diff( Fp[0],  9.15187922867e-01 ) > 1.e-10 );
  s += ( frac_diff(  G[0],  1.03120585919e+00 ) > 1.e-10 );
  s += ( frac_diff( Gp[0],  2.19463267175e-01 ) > 1.e-10 );
  gsl_test(s, "  coulomb(0, 1, 9.2)");
  PRINT(0);
  status += s;

  lam_min = 0.0;
  eta = 10.0;
  x = 10.0;
  gsl_sf_coulomb_wave_FGp_impl(lam_min, kmax, eta, x, F, Fp, G, Gp, &Fe, &Ge);
  s = 0;
  s += ( frac_diff(  F[0],  1.62627112503e-03 ) > 1.e-10 );
  s += ( frac_diff( Fp[0],  1.70604763209e-03 ) > 1.e-10 );
  s += ( frac_diff(  G[0],  3.07873216608e+02 ) > 1.e-10 );
  s += ( frac_diff( Gp[0], -2.91927723806e+02 ) > 1.e-10 );
  gsl_test(s, "  coulomb(0, 10, 10)");
  PRINT(0);
  status += s;

  return status;
}

int check_coupling(void)
{
  double y;
  int status = 0;
  int s;

  gsl_sf_coupling_3j_impl(0, 1, 1, 0, 1, -1, &y);
  s = 0;
  s += ( frac_diff( y, sqrt(1.0/2.0) ) > 1.0e-14 );
  gsl_test(s, "  3j(ja=0    jb=1/2  j=1/2  ma=0    mb=1/2  m=-1/2)");
  status += s;

  gsl_sf_coupling_3j_impl(1, 1, 2, 1, -1, 0, &y);
  s = 0;
  s += ( frac_diff( y, sqrt(1.0/6.0) ) > 1.0e-14 );
  gsl_test(s, "  3j(ja=1/2  jb=1/2  j=1    ma=1/2  mb=-1/2 m=0)");
  status += s;

  gsl_sf_coupling_3j_impl(2, 4, 6, 0, 2, -2, &y);
  s = 0;
  s += ( frac_diff( y, sqrt(8.0/105.0) ) > 1.0e-14 );
  gsl_test(s, "  3j(ja=1    jb=2    j=3    ma=0    mb=1    m=-1)");
  status += s;

  gsl_sf_coupling_6j_impl(2, 2, 4, 2, 2, 2, &y);
  s = 0;
  s += ( frac_diff( y, 1.0/6.0 ) > 1.0e-14 );
  gsl_test(s, "  6j(ja=1   jb=1   jc=2    jd=1     je=1     jf=1)");
  status += s;

  gsl_sf_coupling_9j_impl(4, 2, 4, 3, 3, 2, 1, 1, 2, &y);
  s = 0;
  s += ( frac_diff( y, -0.040824829046386 ) > 1.0e-13 );
  gsl_test(s, "  9j(ja=2   jb=1   jc=2    jd=3/2   je=3/2   jf=1   jg=1/2  jh=1/2  ji=1)");
  status += s;

  gsl_sf_coupling_9j_impl(8, 4, 10, 7, 3, 8, 1, 1, 2, &y);
  s = 0;
  s += ( frac_diff( y, 0.025458753860866 ) > 1.0e-13 );
  gsl_test(s, "  9j(ja=4   jb=2   jc=5    jd=7/2   je=3/2   jf=4   jg=1/2  jh=1/2  ji=1)");
  status += s;

  return status;
}

int check_gamma(void)
{
  double zr, zi, lg_r, lg_i;
  int status = 0;
  int s;
  
  s = 0;
  s += ( frac_diff( gsl_sf_lngamma(0.1),    2.252712651734205 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_lngamma(100.0),  359.1342053695753 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_lngamma");
  status += s;
  
  s = 0;
  zr = 5.0;
  zi = 2.0;
  gsl_sf_lngamma_complex_impl(zr, zi, &lg_r, &lg_i);
  s += ( frac_diff( lg_r, 2.7487017561338026749 ) > 1.0e-14 );
  s += ( frac_diff( lg_i, 3.0738434100497007915 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_lngamma_complex_impl(5 + 2 I)");
  status += s;
  
  s = 0;
  zr = 100.0;
  zi = 100.0;
  gsl_sf_lngamma_complex_impl(zr, zi, &lg_r, &lg_i);
  s += ( frac_diff( lg_r, 315.07804459949331323 ) > 1.0e-14 );
  s += ( frac_diff( lg_i, 2.0821801804113110099 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_lngamma_complex_impl(100 + 100 I)");
  status += s;
  
  s = 0;
  zr =   100.0;
  zi = -1000.0;
  gsl_sf_lngamma_complex_impl(zr, zi, &lg_r, &lg_i);
  s += ( frac_diff( lg_r, -882.3920483010362817000 ) > 1.0e-14 );
  s += ( frac_diff( lg_i,   -2.1169293725678813270 ) > 1.0e-12 );
  gsl_test(s, "  gsl_sf_lngamma_complex_impl(100 - 1000 I)");
  status += s;

  s = 0;
  zr = -100.0;
  zi =   -1.0;
  gsl_sf_lngamma_complex_impl(zr, zi, &lg_r, &lg_i);
  s += ( frac_diff( lg_r, -365.0362469529239516000 ) > 1.0e-14 );
  s += ( frac_diff( lg_i,   -3.0393820262864361140 ) > 1.0e-12 );
  gsl_test(s, "  gsl_sf_lngamma_complex_impl(-1000 - I)");
  status += s;

  return status;
}

int check_gegen(void)
{
  double y;
  int status = 0;
  int s;

  s = 0;
  s += ( frac_diff( gsl_sf_gegenpoly_1(-0.2,   1.0),  -0.4 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_1( 0.0,   1.0),   2.0 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_1( 1.0,   1.0),   2.0 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_1( 1.0,   0.5),   1.0 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_1( 5.0,   1.0),  10.0 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_1( 100.0, 0.5), 100.0 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_gegenpoly_1");
  status += s;
  
  s = 0;
  s += ( frac_diff( gsl_sf_gegenpoly_2(-0.2,   0.5),   0.12 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_2( 0.0,   1.0),   1.00 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_2( 1.0,   1.0),   3.00 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_2( 1.0,   0.1),  -0.96 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_2( 5.0,   1.0),   55.0 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_2( 100.0, 0.5), 4950.0 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_gegenpoly_2");
  status += s;
  
  s = 0;
  s += ( frac_diff( gsl_sf_gegenpoly_3(-0.2,   0.5),      0.112 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_3( 0.0,   1.0),   -2.0/3.0 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_3( 1.0,   1.0),      4.000 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_3( 1.0,   0.1),     -0.392 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_3( 5.0,   1.0),    220.000 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_3( 100.0, 0.5), 161600.000 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_gegenpoly_3");
  status += s;

  s = 0;
  s += ( frac_diff( gsl_sf_gegenpoly_n(1,       1.0, 1.0),  2.000		    ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_n(10,      1.0, 1.0), 11.000		    ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_n(10,      1.0, 0.1), -0.4542309376	    ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_n(10,      5.0, 1.0),  9.23780e+4  	    ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_n(10,    100.0, 0.5),  1.5729338392690000e+13  ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_n(1000,  100.0, 1.0),  3.3353666135627322e+232 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_n(100,  2000.0, 1.0),  5.8753432034937579e+202 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_n(103,   207.0, 2.0),  1.4210272202235983e+145 ) > 1.e-14 );
  s += ( frac_diff( gsl_sf_gegenpoly_n(103,    -0.4, 0.3), -1.64527498094522e-04    ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_gegenpoly_n");
  status += s;
  
  s += 0;
  s += (gsl_sf_gegenpoly_n_impl(103, -0.5,  1.0, &y) != GSL_EDOM);
  s += (gsl_sf_gegenpoly_n_impl(103, -0.51, 1.0, &y) != GSL_EDOM);
  gsl_test(s, "  gsl_sf_gegenpoly_n_impl: trap lambda <= -1/2");
  status += s;

  s += 0;
  s += (gsl_sf_gegenpoly_n_impl(-2, 0.0,  1.0, &y) != GSL_EDOM);
  gsl_test(s, "  gsl_sf_gegenpoly_n_impl: trap n < 0");
  status += s;
  
  return status;
}


int check_hyperg(void)
{
  int status = 0;
  int s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_1F1(1, 1, 0.5),  1.6487212707001281468    ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_hyperg_1F1(8, 1, 0.5),  13.108875178030540372    ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_hyperg_1F1(8, 1, 8),    5.481422453671217135e+7  ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_hyperg_1F1(1, 8, 8),    4.918996932358889820     ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_1F1");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_2F1(1, 1, 1, 0.5),  2.0000000000000000000    ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_hyperg_2F1(8, 8, 1, 0.5),  1.2451584000000000000e+7 ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_hyperg_2F1(8, 8, 5, 0.5),  4205.714285714285714     ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_2F1");
  status += s;

  /* FIXME: the "true" values here may not be so good */
  s = 0;
  s += ( frac_diff(gsl_sf_hyperg_2F0(0.01, 1.0, -0.02), 0.999803886708565    ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_hyperg_2F0(0.1,  0.5, -0.02), 0.999015947934831    ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_hyperg_2F0(1,   1, -0.02),   0.980755496569062     ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_hyperg_2F0(8,   8, -0.02),   0.3299059284994299    ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_hyperg_2F0(50, 50, -0.02),   2.688995263773233e-13 ) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hyperg_2F0");
  status += s;

  return status;
}

int check_pow_int(void)
{
  int status = 0;
  int s;
  
  s = 0;
  s += ( frac_diff( gsl_sf_pow_int(2.0, 3), 8.0 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_pow_int(2.0, 3)");
  status += s;
  
  s = 0;
  s += ( frac_diff( gsl_sf_pow_int(-2.0, 3), -8.0 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_pow_int(-2.0, 3)");
  status += s;
  
  s = 0;
  s += ( frac_diff( gsl_sf_pow_int(2.0, -3), 1.0/8.0 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_pow_int(2.0, -3)");
  status += s;
  
  s = 0;
  s += ( frac_diff( gsl_sf_pow_int(-2.0, -3), -1.0/8.0 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_pow_int(-2.0, -3)");
  status += s;
  
  return status;
}

int check_psi(void)
{
  int status = 0;
  int s;
  
  s = 0;
  s += ( frac_diff(gsl_sf_psi_int(5), 1.5061176684318004727 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_psi_int(5)");
  status += s;
  
  s = 0;
  s += ( frac_diff(gsl_sf_psi_int(5000), 8.517093188082904107 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_psi_int(5000)");
  status += s;
  
  s = 0;
  s += ( frac_diff(gsl_sf_psi(5.0), 1.5061176684318004727 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_psi(5.0)");
  status += s;
  
  s = 0;
  s += ( frac_diff(gsl_sf_psi(5000.0), 8.517093188082904107 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_psi(5000.0)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_psi(-100.5), 4.615124601338064117 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_psi(-100.5)");
  status += s;

  s = 0;
  s += ( frac_diff(gsl_sf_psi(-1.0e+5-0.5), 11.512935464924395337 ) > 1.0e-10 );
  gsl_test(s, "  gsl_sf_psi(-1.0e+5-0.5)");
  status += s;

  return status;
}

int check_trig(void)
{
  double zr, zi;
  double yr, yi;
  double x, y;
  double theta;
  int status = 0;
  int s;

  yr = 1.0;
  yi = 5.0;
  gsl_sf_complex_sin_impl(yr, yi, &zr, &zi);
  s = 0;
  s += ( frac_diff( zr, 62.44551846769653403 ) > 1.0e-14 );
  s += ( frac_diff( zi, 40.09216577799840254 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_complex_sin_impl(1 + 5 I)");
  status += s;
  
  yr = 1.0;
  yi = 5.0;
  gsl_sf_complex_cos_impl(yr, yi, &zr, &zi);
  s = 0;
  s += ( frac_diff( zr,  40.09580630629882573 ) > 1.0e-14 );
  s += ( frac_diff( zi, -62.43984868079963017 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_complex_cos_impl(1 + 5 I)");
  status += s;

  yr =   1.0;
  yi = 100.0;
  gsl_sf_complex_logsin_impl(yr, yi, &zr, &zi);
  s = 0;
  s += ( frac_diff( zr, 99.3068528194400546900 ) > 1.0e-14 );
  s += ( frac_diff( zi,  0.5707963267948966192 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_complex_logsin_impl(1 + 100 I)");
  status += s;
  
  yr =    1.0;
  yi = -100.0;
  gsl_sf_complex_logsin_impl(yr, yi, &zr, &zi);
  s = 0;
  s += ( frac_diff( zr,  99.3068528194400546900 ) > 1.0e-14 );
  s += ( frac_diff( zi,  -0.5707963267948966192 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_complex_logsin_impl(1 - 100 I)");
  status += s;

  yr = 5.0;
  yi = 5.0;
  gsl_sf_complex_logsin_impl(yr, yi, &zr, &zi);
  s = 0;
  s += ( frac_diff( zr, 4.3068909128079757420 ) > 1.0e-14 );
  s += ( frac_diff( zi, 2.8540063315538773952 ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_complex_logsin_impl(5 + 5 I)");
  status += s;

  gsl_sf_polar_to_rect_impl(10.0, M_PI/6.0, &x, &y);
  s = 0;
  s += ( frac_diff( x, 10.0 * sqrt(3) / 2.0 ) > 1.0e-14 );
  s += ( frac_diff( y, 10.0 * 0.5           ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_polar_to_rect_impl(10, Pi/6)");
  status += s;
  
  gsl_sf_polar_to_rect_impl(10.0, -2.0/3.0*M_PI, &x, &y);
  s = 0;
  s += ( frac_diff( x, 10.0 * (-0.5)           ) > 1.0e-14 );
  s += ( frac_diff( y, 10.0 * (-sqrt(3) / 2.0) ) > 1.0e-14 );
  gsl_test(s, "  gsl_sf_polar_to_rect_impl(10, -2/3 Pi)");
  status += s;
  
  theta = 5.0*M_PI + M_PI/2.0;
  gsl_sf_angle_restrict_pos_impl(&theta, GSL_SQRT_MACH_EPS);
  s = 0;
  s += ( frac_diff( theta, 3.0/2.0*M_PI ) > 1.e-14 );
  gsl_test(s, "  gsl_angle_restrict_pos_impl: theta =  11/2 Pi");
  status += s;

  theta = -5.0*M_PI - M_PI/2.0;
  gsl_sf_angle_restrict_pos_impl(&theta, GSL_SQRT_MACH_EPS);
  s = 0;
  s += ( frac_diff( theta, M_PI/2.0 ) > 1.e-14 );
  gsl_test(s, "  gsl_angle_restrict_pos_impl: theta = -11/2 Pi");
  status += s;

  theta = 5.0*M_PI + M_PI/2.0;
  gsl_sf_angle_restrict_symm_impl(&theta, GSL_SQRT_MACH_EPS);
  s = 0;
  s += ( frac_diff( theta, -M_PI/2.0 ) > 1.e-14 );
  gsl_test(s, "  gsl_angle_restrict_symm_impl: theta =  11/2 Pi");
  status += s;

  theta = -5.0*M_PI - M_PI/2.0;
  gsl_sf_angle_restrict_symm_impl(&theta, GSL_SQRT_MACH_EPS);
  s = 0;
  s += ( frac_diff( theta, M_PI/2.0 ) > 1.e-14 );
  gsl_test(s, "  gsl_angle_restrict_symm_impl: theta = -11/2 Pi");
  status += s;

  theta =  5.0*M_PI - M_PI/2.0;
  gsl_sf_angle_restrict_symm_impl(&theta, GSL_SQRT_MACH_EPS);
  s = 0;
  s += ( frac_diff( theta, M_PI/2.0 ) > 1.e-14 );
  gsl_test(s, "  gsl_angle_restrict_symm_impl: theta = -9/2 Pi");
  status += s;

  theta =  3.0/2.0*M_PI;
  gsl_sf_angle_restrict_symm_impl(&theta, GSL_SQRT_MACH_EPS);
  s = 0;
  s += ( frac_diff( theta, -M_PI/2.0 ) > 1.e-14 );
  gsl_test(s, "  gsl_angle_restrict_symm_impl: theta =  3/2 Pi");
  status += s;

  theta = -3.0/2.0*M_PI;
  gsl_sf_angle_restrict_symm_impl(&theta, GSL_SQRT_MACH_EPS);
  s = 0;
  s += ( frac_diff( theta, M_PI/2.0 ) > 1.e-14 );
  gsl_test(s, "  gsl_angle_restrict_symm_impl: theta = -3/2 Pi");
  status += s;

  for(theta=1000.0; theta < 1.e20 ; theta *= 10.0) {
    double nyak = 2.0*M_PI*theta + 0.5*M_PI;
    gsl_sf_angle_restrict_symm_impl(&nyak, 1.0);
    printf("%24.20g   %24.20g   %24.20g\n", theta, nyak, (nyak-0.5*M_PI)/(nyak+0.5*M_PI));
  }
  return status;
}

int check_zeta(void)
{
  int status = 0;
  int s;

  /** functions */

  s = 0;
  s += ( frac_diff(gsl_sf_zeta_int(-5),  -0.003968253968253968253968) > 1.e-14 );
  s += ( frac_diff(gsl_sf_zeta_int(-61), -3.30660898765775767257e+34) > 1.e-14 );
  s += ( frac_diff(gsl_sf_zeta_int(5),   1.0369277551433699263313655) > 1.e-14 );
  s += ( frac_diff(gsl_sf_zeta_int(31),  1.0000000004656629065033784) > 1.e-14 );
  gsl_test(s, "  gsl_sf_zeta_int");
  status += s;
  
  s = 0;
  s += ( frac_diff(gsl_sf_zeta(-51),   9.68995788746359406565e+24) > 1.e-14 );
  s += ( frac_diff(gsl_sf_zeta(-5),   -0.003968253968253968253968) > 1.e-14 );
  s += ( frac_diff(gsl_sf_zeta(-0.5), -0.207886224977354566017307) > 1.e-14 );
  s += ( frac_diff(gsl_sf_zeta(0.5),  -1.460354508809586812889499) > 1.e-14 );
  s += ( frac_diff(gsl_sf_zeta(5),     1.036927755143369926331365) > 1.e-14 );
  s += ( frac_diff(gsl_sf_zeta(25.5),  1.000000021074106110269959) > 1.e-14 );
  gsl_test(s, "  gsl_sf_zeta");
  status += s;
  
  s = 0;
  s += ( frac_diff(gsl_sf_hzeta(2,  1.0),  1.6449340668482264365   ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_hzeta(5,  1.0),  1.0369277551433699263   ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_hzeta(9,  0.1),  1.0000000004253980e+9   ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_hzeta(2, 10.0),  0.1051663356816857461   ) > 1.e-14 );
  s += ( frac_diff(gsl_sf_hzeta(5, 10.0),  0.0000304137986764702764) > 1.e-14 );
  gsl_test(s, "  gsl_sf_hzeta");
  status += s;
  
  return status;
}


int main(int argc, char * argv[])
{
  gsl_test(check_airy(),       "Airy Functions");
  gsl_test(check_bessel(),     "Bessel Functions");
  gsl_test(check_clausen(),    "Clausen Integral");
  gsl_test(check_coulomb(),    "Coulomb Wave Functions");
  gsl_test(check_coupling(),   "Coupling Coefficients");

  gsl_test(check_gamma(),      "Gamma Functions");
  gsl_test(check_gegen(),      "Gegenbauer Polynomials");
  gsl_test(check_hyperg(),     "Hypergeometric Functions");

  gsl_test(check_pow_int(),    "Integer Powers");
  gsl_test(check_psi(),        "Psi Functions");
  gsl_test(check_trig(),       "Trigonometric and Related Functions");
  gsl_test(check_zeta(),       "Zeta Functions");

  gsl_test_summary();

  return 0;  
}
