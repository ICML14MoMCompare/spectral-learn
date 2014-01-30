/* sys/test.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

#include <config.h>
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_test.h>

#include <gsl/gsl_ieee_utils.h>

int
main (void)
{
  double y, y_expected;

  gsl_ieee_env_setup ();

  /* Test for expm1 */

  y = gsl_expm1 (0.0); y_expected = 0.0;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_expm1(0.0)");

  y = gsl_expm1 (1e-10); y_expected = 1.000000000050000000002e-10;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_expm1(1e-10)");

  y = gsl_expm1 (-1e-10); y_expected = -9.999999999500000000017e-11;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_expm1(-1e-10)");

  y = gsl_expm1 (0.1); y_expected = 0.1051709180756476248117078264902;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_expm1(0.1)");

  y = gsl_expm1 (-0.1); y_expected = -0.09516258196404042683575094055356;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_expm1(-0.1)");

  y = gsl_expm1 (10.0); y_expected = 22025.465794806716516957900645284;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_expm1(10.0)");

  y = gsl_expm1 (-10.0); y_expected = -0.99995460007023751514846440848444;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_expm1(-10.0)");
   
  /* Test for log1p */

  y = gsl_log1p (0.0); y_expected = 0.0;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_log1p(0.0)");

  y = gsl_log1p (1e-10); y_expected = 9.9999999995000000000333333333308e-11;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_log1p(1e-10)");

  y = gsl_log1p (0.1); y_expected = 0.095310179804324860043952123280765;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_log1p(0.1)");

  y = gsl_log1p (10.0); y_expected = 2.3978952727983705440619435779651;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_log1p(10.0)");

  /* Test for gsl_hypot */

  y = gsl_hypot (0.0, 0.0) ; y_expected = 0.0;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(0.0, 0.0)");

  y = gsl_hypot (1e-10, 1e-10) ; y_expected = 1.414213562373095048801688e-10;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(1e-10, 1e-10)");

  y = gsl_hypot (1e-38, 1e-38) ; y_expected = 1.414213562373095048801688e-38;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(1e-38, 1e-38)");

  y = gsl_hypot (1e-10, -1.0) ; y_expected = 1.000000000000000000005;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(1e-10, -1)");

  y = gsl_hypot (-1.0, 1e-10) ; y_expected = 1.000000000000000000005;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(-1, 1e-10)");

  y = gsl_hypot (1e307, 1e301) ; y_expected = 1.000000000000499999999999e307;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(1e307, 1e301)");

  y = gsl_hypot (1e301, 1e307) ; y_expected = 1.000000000000499999999999e307;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(1e301, 1e307)");

  y = gsl_hypot (1e307, 1e307) ; y_expected = 1.414213562373095048801688e307;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_hypot(1e307, 1e307)");


  /* Test for acosh */

  y = gsl_acosh (1.0); y_expected = 0.0;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_acosh(1.0)");

  y = gsl_acosh (1.1); y_expected = 4.435682543851151891329110663525e-1;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_acosh(1.1)");

  y = gsl_acosh (10.0); y_expected = 2.9932228461263808979126677137742e0;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_acosh(10.0)");

  y = gsl_acosh (1e10); y_expected = 2.3718998110500402149594646668302e1;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_acosh(1e10)");

  /* Test for asinh */

  y = gsl_asinh (0.0); y_expected = 0.0;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_asinh(0.0)");

  y = gsl_asinh (1e-10); y_expected = 9.9999999999999999999833333333346e-11;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_asinh(1e-10)");

  y = gsl_asinh (-1e-10); y_expected = -9.9999999999999999999833333333346e-11;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_asinh(1e-10)");

  y = gsl_asinh (0.1); y_expected = 9.983407889920756332730312470477e-2;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_asinh(0.1)");

  y = gsl_asinh (-0.1); y_expected = -9.983407889920756332730312470477e-2;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_asinh(-0.1)");

  y = gsl_asinh (1.0); y_expected = 8.8137358701954302523260932497979e-1;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_asinh(1.0)");

  y = gsl_asinh (-1.0); y_expected = -8.8137358701954302523260932497979e-1;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_asinh(-1.0)");

  y = gsl_asinh (10.0); y_expected = 2.9982229502979697388465955375965e0;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_asinh(10)");

  y = gsl_asinh (-10.0); y_expected = -2.9982229502979697388465955375965e0;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_asinh(-10)");

  y = gsl_asinh (1e10); y_expected = 2.3718998110500402149599646668302e1;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_asinh(1e10)");

  y = gsl_asinh (-1e10); y_expected = -2.3718998110500402149599646668302e1;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_asinh(-1e10)");

  /* Test for atanh */

  y = gsl_atanh (0.0); y_expected = 0.0;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_atanh(0.0)");

  y = gsl_atanh (1e-20); y_expected = 1e-20;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_atanh(1e-20)");

  y = gsl_atanh (-1e-20); y_expected = -1e-20;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_atanh(-1e-20)");

  y = gsl_atanh (0.1); y_expected = 1.0033534773107558063572655206004e-1;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_atanh(0.1)");

  y = gsl_atanh (-0.1); y_expected = -1.0033534773107558063572655206004e-1;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_atanh(-0.1)");

  y = gsl_atanh (0.9); y_expected = 1.4722194895832202300045137159439e0;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_atanh(0.9)");

  y = gsl_atanh (-0.9); y_expected = -1.4722194895832202300045137159439e0;
  gsl_test_rel (y, y_expected, 1e-15, "gsl_atanh(0.9)");

  /* Test for pow_int */

  y = gsl_pow_2 (-3.14); y_expected = pow(-3.14, 2.0);
  gsl_test_rel (y, y_expected, 1e-15, "gsl_pow_2(-3.14)");

  y = gsl_pow_3 (-3.14); y_expected = pow(-3.14, 3.0);
  gsl_test_rel (y, y_expected, 1e-15, "gsl_pow_3(-3.14)");

  y = gsl_pow_4 (-3.14); y_expected = pow(-3.14, 4.0);
  gsl_test_rel (y, y_expected, 1e-15, "gsl_pow_4(-3.14)");

  y = gsl_pow_5 (-3.14); y_expected = pow(-3.14, 5.0);
  gsl_test_rel (y, y_expected, 1e-15, "gsl_pow_5(-3.14)");

  y = gsl_pow_6 (-3.14); y_expected = pow(-3.14, 6.0);
  gsl_test_rel (y, y_expected, 1e-15, "gsl_pow_6(-3.14)");

  y = gsl_pow_7 (-3.14); y_expected = pow(-3.14, 7.0);
  gsl_test_rel (y, y_expected, 1e-15, "gsl_pow_7(-3.14)");

  y = gsl_pow_8 (-3.14); y_expected = pow(-3.14, 8.0);
  gsl_test_rel (y, y_expected, 1e-15, "gsl_pow_8(-3.14)");

  y = gsl_pow_9 (-3.14); y_expected = pow(-3.14, 9.0);
  gsl_test_rel (y, y_expected, 1e-15, "gsl_pow_9(-3.14)");

  { 
    int n;
    for (n = -9; n < 10; n++) {
      y = gsl_pow_int (-3.14, n); y_expected = pow(-3.14, n);
      gsl_test_rel (y, y_expected, 1e-15, "gsl_pow_n(-3.14,%d)", n);
    }
  }

  /* Test for isinf, isnan, finite*/

  {
    double zero, one, inf, nan;
    int s;

    zero = 0.0;
    one = 1.0;
    inf = exp(1.0e10);
    nan = inf / inf;
    
    s = gsl_isinf(zero);
    gsl_test_int (s, 0, "gsl_isinf(0)");
    
    s = gsl_isinf(one);
    gsl_test_int (s, 0, "gsl_isinf(1)");
    
    s = gsl_isinf(inf);
    gsl_test_int (s, 1, "gsl_isinf(inf)");

    s = gsl_isinf(-inf);
    gsl_test_int (s, -1, "gsl_isinf(-inf)");
    
    s = gsl_isinf(nan);
    gsl_test_int (s, 0, "gsl_isinf(nan)");


    s = gsl_isnan(zero);
    gsl_test_int (s, 0, "gsl_isnan(0)");
    
    s = gsl_isnan(one);
    gsl_test_int (s, 0, "gsl_isnan(1)");
    
    s = gsl_isnan(inf);
    gsl_test_int (s, 0, "gsl_isnan(inf)");
    
    s = gsl_isnan(nan);
    gsl_test_int (s, 1, "gsl_isnan(nan)");


    s = gsl_finite(zero);
    gsl_test_int (s, 1, "gsl_finite(0)");
    
    s = gsl_finite(one);
    gsl_test_int (s, 1, "gsl_finite(1)");
    
    s = gsl_finite(inf);
    gsl_test_int (s, 0, "gsl_finite(inf)");
    
    s = gsl_finite(nan);
    gsl_test_int (s, 0, "gsl_finite(nan)");
  }

  {
    double x = gsl_fdiv (2.0, 3.0);
    gsl_test_rel (x, 2.0/3.0, 4*GSL_DBL_EPSILON, "gsl_fdiv(2,3)");
  }

  exit (gsl_test_summary ());
}


