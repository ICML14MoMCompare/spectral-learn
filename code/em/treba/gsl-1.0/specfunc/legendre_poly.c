/* specfunc/legendre_poly.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
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

/* Author:  G. Jungman */

#include <config.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include "gsl_sf_bessel.h"
#include "gsl_sf_exp.h"
#include "gsl_sf_gamma.h"
#include "gsl_sf_log.h"
#include "gsl_sf_pow_int.h"
#include "gsl_sf_legendre.h"

#include "error.h"

/*-*-*-*-*-*-*-*-*-*-*-* Functions with Error Codes *-*-*-*-*-*-*-*-*-*-*-*/

int
gsl_sf_legendre_P1_e(double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  {
    result->val = x;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
}

int
gsl_sf_legendre_P2_e(double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  {
    result->val = 0.5*(3.0*x*x - 1.0);
    result->err = GSL_DBL_EPSILON * (fabs(3.0*x*x) + 1.0);
    return GSL_SUCCESS;
  }
}

int
gsl_sf_legendre_P3_e(double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  {
    result->val = 0.5*x*(5.0*x*x - 3.0);
    result->err = GSL_DBL_EPSILON * (fabs(result->val) + 0.5 * fabs(x) * (fabs(5.0*x*x) + 3.0));
    return GSL_SUCCESS;
  }
}


int
gsl_sf_legendre_Pl_e(const int l, const double x, gsl_sf_result * result)
{ 
  /* CHECK_POINTER(result) */

  if(l < 0 || x < -1.0 || x > 1.0) {
    DOMAIN_ERROR(result);
  }
  else if(l == 0) {
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(l == 1) {
    result->val = x;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(l == 2) {
    result->val = 0.5 * (3.0*x*x - 1.0);
    result->err = 3.0 * GSL_DBL_EPSILON * fabs(result->val);
    return GSL_SUCCESS;
  }
  else if(x == 1.0) {
    result->val = 1.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(x == -1.0) {
    result->val = ( GSL_IS_ODD(l) ? -1.0 : 1.0 );
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else if(l < 100000) {
    /* Compute by upward recurrence on l.
     */
    double p_mm   = 1.0;     /* P_0(x) */
    double p_mmp1 = x;	    /* P_1(x) */
    double p_ell = p_mmp1;
    int ell;

    for(ell=2; ell <= l; ell++){
      p_ell = (x*(2*ell-1)*p_mmp1 - (ell-1)*p_mm) / ell;
      p_mm = p_mmp1;
      p_mmp1 = p_ell;
    }

    result->val = p_ell;
    result->err = (0.5 * ell + 1.0) * GSL_DBL_EPSILON * fabs(p_ell);
    return GSL_SUCCESS;
  }
  else {
    /* Asymptotic expansion.
     * [Olver, p. 473]
     */
    double u  = l + 0.5;
    double th = acos(x);
    gsl_sf_result J0;
    gsl_sf_result Jm1;
    int stat_J0  = gsl_sf_bessel_J0_e(u*th, &J0);
    int stat_Jm1 = gsl_sf_bessel_Jn_e(-1, u*th, &Jm1);
    double pre;
    double B00;
    double c1;

    /* B00 = 1/8 (1 - th cot(th) / th^2
     * pre = sqrt(th/sin(th))
     */
    if(th < GSL_ROOT4_DBL_EPSILON) {
      B00 = (1.0 + th*th/15.0)/24.0;
      pre = 1.0 + th*th/12.0;
    }
    else {
      double sin_th = sqrt(1.0 - x*x);
      double cot_th = x / sin_th;
      B00 = 1.0/8.0 * (1.0 - th * cot_th) / (th*th);
      pre = sqrt(th/sin_th);
    }

    c1 = th/u * B00;

    result->val  = pre * (J0.val + c1 * Jm1.val);
    result->err  = pre * (J0.err + fabs(c1) * Jm1.err);
    result->err += GSL_SQRT_DBL_EPSILON * fabs(result->val);

    return GSL_ERROR_SELECT_2(stat_J0, stat_Jm1);
  }
}


int
gsl_sf_legendre_Pl_array(const int lmax, const double x, double * result_array)
{
  /* CHECK_POINTER(result_array) */

  if(lmax < 0 || x < -1.0 || x > 1.0) {
    GSL_ERROR ("domain error", GSL_EDOM);
  }
  else if(lmax == 0) {
    result_array[0] = 1.0;
    return GSL_SUCCESS;
  }
  else if(lmax == 1) {
    result_array[0] = 1.0;
    result_array[1] = x;
    return GSL_SUCCESS;
  }
  else {
    double p_mm   = 1.0;    /* P_0(x) */
    double p_mmp1 = x;	    /* P_1(x) */
    double p_ell = p_mmp1;
    int ell;

    result_array[0] = 1.0;
    result_array[1] = x;

    for(ell=2; ell <= lmax; ell++){
      p_ell = (x*(2*ell-1)*p_mmp1 - (ell-1)*p_mm) / ell;
      p_mm = p_mmp1;
      p_mmp1 = p_ell;
      result_array[ell] = p_ell;
    }

    return GSL_SUCCESS;
  }
}


int
gsl_sf_legendre_Plm_e(const int l, const int m, const double x, gsl_sf_result * result)
{
  /* If l is large and m is large, then we have to worry
   * about overflow. Calculate an approximate exponent which
   * measures the normalization of this thing.
   */
  double dif = l-m;
  double sum = l+m;
  double exp_check = 0.5 * log(2.0*l+1.0) 
                   + 0.5 * dif * (log(dif)-1.0)
                   - 0.5 * sum * (log(sum)-1.0);

  /* CHECK_POINTER(result) */

  if(m < 0 || l < m || x < -1.0 || x > 1.0) {
    DOMAIN_ERROR(result);
  }
  else if(exp_check < GSL_LOG_DBL_MIN + 10.0){
    /* Bail out. */
    OVERFLOW_ERROR(result);
  }
  else {
    /* Account for the error due to the
     * representation of 1-x.
     */
    const double err_amp = 1.0 / (GSL_DBL_EPSILON + fabs(1.0-fabs(x)));

    double p_mm;     /* P_m^m(x) */
    double p_mmp1;   /* P_{m+1}^m(x) */

    /* Calculate P_m^m from the analytic result:
     *          P_m^m(x) = (-1)^m (2m-1)!! (1-x^2)^(m/2) , m > 0
     */
    p_mm = 1.0;
    if(m > 0){
      double root_factor = sqrt(1.0-x)*sqrt(1.0+x);
      double fact_coeff = 1.0;
      int i;
      for(i=1; i<=m; i++) {
        p_mm *= -fact_coeff * root_factor;
        fact_coeff += 2.0;
      }
    }

    /* Calculate P_{m+1}^m. */
    p_mmp1 = x * (2*m + 1) * p_mm;

    if(l == m){
      result->val = p_mm;
      result->err = err_amp * 2.0 * GSL_DBL_EPSILON * fabs(p_mm);
      return GSL_SUCCESS;
    }
    else if(l == m + 1) {
      result->val = p_mmp1;
      result->err = err_amp * 2.0 * GSL_DBL_EPSILON * fabs(p_mmp1);
      return GSL_SUCCESS;
    }
    else{
      double p_ell = 0.0;
      int ell;
    
      /* Compute P_l^m, l > m+1 by upward recurrence on l. */
      for(ell=m+2; ell <= l; ell++){
        p_ell = (x*(2*ell-1)*p_mmp1 - (ell+m-1)*p_mm) / (ell-m);
        p_mm = p_mmp1;
        p_mmp1 = p_ell;
      }

      result->val = p_ell;
      result->err = err_amp * (0.5*(l-m) + 1.0) * GSL_DBL_EPSILON * fabs(p_ell);

      return GSL_SUCCESS;
    }
  }
}


int
gsl_sf_legendre_Plm_array(const int lmax, const int m, const double x, double * result_array)
{
  /* If l is large and m is large, then we have to worry
   * about overflow. Calculate an approximate exponent which
   * measures the normalization of this thing.
   */
  double dif = lmax-m;
  double sum = lmax+m;
  double exp_check = 0.5 * log(2.0*lmax+1.0) 
                     + 0.5 * dif * (log(dif)-1.0)
                     - 0.5 * sum * (log(sum)-1.0);

  /* CHECK_POINTER(result_array) */

  if(m < 0 || lmax < m || x < -1.0 || x > 1.0) {
    GSL_ERROR ("error", GSL_EDOM);
  }
  else if(m > 0 && (x == 1.0 || x == -1.0)) {
    int ell;
    for(ell=m; ell<=lmax; ell++) result_array[ell-m] = 0.0;
    return GSL_SUCCESS;
  }
  else if(exp_check < GSL_LOG_DBL_MIN + 10.0){
    /* Bail out.
     */
    GSL_ERROR ("error", GSL_EOVRFLW);
  }
  else {
    double p_mm;                 /* P_m^m(x)     */
    double p_mmp1;               /* P_{m+1}^m(x) */

    /* Calculate P_m^m from the analytic result:
     *          P_m^m(x) = (-1)^m (2m-1)!! (1-x^2)^(m/2) , m > 0
     */
    p_mm = 1.0;
    if(m > 0){
      double root_factor = sqrt(1.0-x)*sqrt(1.0+x);
      double fact_coeff = 1.0;
      int i;
      for(i=1; i<=m; i++){
        p_mm *= -fact_coeff * root_factor;
        fact_coeff += 2.0;
      }
    }

    /* Calculate P_{m+1}^m. */
    p_mmp1 = x * (2*m + 1) * p_mm;

    if(lmax == m){
      result_array[0] = p_mm;
      return GSL_SUCCESS;
    }
    else if(lmax == m + 1) {
      result_array[0] = p_mm;
      result_array[1] = p_mmp1;
      return GSL_SUCCESS;
    }
    else{
      double p_ell;
      int ell;

      result_array[0] = p_mm;
      result_array[1] = p_mmp1;

      /* Compute P_l^m, l >= m+2, by upward recursion on l. */
      for(ell=m+2; ell <= lmax; ell++){
        p_ell = (x*(2*ell-1)*p_mmp1 - (ell+m-1)*p_mm) / (ell-m);
        p_mm = p_mmp1;
        p_mmp1 = p_ell;
        result_array[ell-m] = p_ell;
      }

      return GSL_SUCCESS;
    }
  }
}


int
gsl_sf_legendre_sphPlm_e(const int l, int m, const double x, gsl_sf_result * result)
{
  /* CHECK_POINTER(result) */

  if(m < 0 || l < m || x < -1.0 || x > 1.0) {
    DOMAIN_ERROR(result);
  }
  else if(m == 0) {
    gsl_sf_result P;
    int stat_P = gsl_sf_legendre_Pl_e(l, x, &P);
    double pre = sqrt((2.0*l + 1.0)/(4.0*M_PI));
    result->val  = pre * P.val;
    result->err  = pre * P.err;
    result->err += 2.0 * GSL_DBL_EPSILON * fabs(result->val);
    return stat_P;
  }
  else if(x == 1.0 || x == -1.0) {
    /* m > 0 here */
    result->val = 0.0;
    result->err = 0.0;
    return GSL_SUCCESS;
  }
  else {
    /* m > 0 and |x| < 1 here */

    /* Starting value for recursion.
     * Y_m^m(x) = sqrt( (2m+1)/(4pi m) gamma(m+1/2)/gamma(m) ) (-1)^m (1-x^2)^(m/2) / pi^(1/4)
     */
    gsl_sf_result lncirc;
    gsl_sf_result lnpoch;
    double lnpre_val;
    double lnpre_err;
    gsl_sf_result ex_pre;
    double sr;
    const double sgn = ( GSL_IS_ODD(m) ? -1.0 : 1.0);
    const double y_mmp1_factor = x * sqrt(2.0*m + 3.0);
    double y_mm, y_mm_err;
    double y_mmp1;
    gsl_sf_log_1plusx_e(-x*x, &lncirc);
    gsl_sf_lnpoch_e(m, 0.5, &lnpoch);  /* Gamma(m+1/2)/Gamma(m) */
    lnpre_val = -0.25*M_LNPI + 0.5 * (lnpoch.val + m*lncirc.val);
    lnpre_err = 0.25*M_LNPI*GSL_DBL_EPSILON + 0.5 * (lnpoch.err + fabs(m)*lncirc.err);
    gsl_sf_exp_err_e(lnpre_val, lnpre_err, &ex_pre);
    sr    = sqrt((2.0+1.0/m)/(4.0*M_PI));
    y_mm   = sgn * sr * ex_pre.val;
    y_mmp1 = y_mmp1_factor * y_mm;
    y_mm_err  = 2.0 * GSL_DBL_EPSILON * fabs(y_mm) + sr * ex_pre.err;
    y_mm_err *= 1.0 + 1.0/(GSL_DBL_EPSILON + fabs(1.0-x));

    if(l == m){
      result->val  = y_mm;
      result->err  = y_mm_err;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(y_mm);
      return GSL_SUCCESS;
    }
    else if(l == m + 1) {
      result->val  = y_mmp1;
      result->err  = fabs(y_mmp1_factor) * y_mm_err;
      result->err += 2.0 * GSL_DBL_EPSILON * fabs(y_mmp1);
      return GSL_SUCCESS;
    }
    else{
      double y_ell = 0.0;
      int ell;

      /* Compute Y_l^m, l > m+1, upward recursion on l. */
      for(ell=m+2; ell <= l; ell++){
        const double rat1 = (double)(ell-m)/(double)(ell+m);
	const double rat2 = (ell-m-1.0)/(ell+m-1.0);
        const double factor1 = sqrt(rat1*(2*ell+1)*(2*ell-1));
        const double factor2 = sqrt(rat1*rat2*(2*ell+1)/(2*ell-3));
        y_ell = (x*y_mmp1*factor1 - (ell+m-1)*y_mm*factor2) / (ell-m);
        y_mm   = y_mmp1;
        y_mmp1 = y_ell;
      }

      result->val  = y_ell;
      result->err  = (0.5*(l-m) + 1.0) * GSL_DBL_EPSILON * fabs(y_ell);
      result->err += fabs(y_mm_err/y_mm) * fabs(y_ell);

      return GSL_SUCCESS;
    }
  }
}


int
gsl_sf_legendre_sphPlm_array(const int lmax, int m, const double x, double * result_array)
{
  /* CHECK_POINTER(result_array) */

  if(m < 0 || lmax < m || x < -1.0 || x > 1.0) {
    GSL_ERROR ("error", GSL_EDOM);
  }
  else if(m > 0 && (x == 1.0 || x == -1.0)) {
    int ell;
    for(ell=m; ell<=lmax; ell++) result_array[ell-m] = 0.0;
    return GSL_SUCCESS;
  }
  else {
    double y_mm;
    double y_mmp1;

    if(m == 0) {
      y_mm   = 0.5/M_SQRTPI;          /* Y00 = 1/sqrt(4pi) */
      y_mmp1 = x * M_SQRT3 * y_mm;
    }
    else {
      /* |x| < 1 here */

      gsl_sf_result lncirc;
      gsl_sf_result lnpoch;
      double lnpre;
      const double sgn = ( GSL_IS_ODD(m) ? -1.0 : 1.0);
      gsl_sf_log_1plusx_e(-x*x, &lncirc);
      gsl_sf_lnpoch_e(m, 0.5, &lnpoch);  /* Gamma(m+1/2)/Gamma(m) */
      lnpre = -0.25*M_LNPI + 0.5 * (lnpoch.val + m*lncirc.val);
      y_mm   = sqrt((2.0+1.0/m)/(4.0*M_PI)) * sgn * exp(lnpre);
      y_mmp1 = x * sqrt(2.0*m + 3.0) * y_mm;
    }

    if(lmax == m){
      result_array[0] = y_mm;
      return GSL_SUCCESS;
    }
    else if(lmax == m + 1) {
      result_array[0] = y_mm;
      result_array[1] = y_mmp1;
      return GSL_SUCCESS;
    }
    else{
      double y_ell;
      int ell;

      result_array[0] = y_mm;
      result_array[1] = y_mmp1;

      /* Compute Y_l^m, l > m+1, upward recursion on l. */
      for(ell=m+2; ell <= lmax; ell++){
        const double rat1 = (double)(ell-m)/(double)(ell+m);
	const double rat2 = (ell-m-1.0)/(ell+m-1.0);
        const double factor1 = sqrt(rat1*(2*ell+1)*(2*ell-1));
        const double factor2 = sqrt(rat1*rat2*(2*ell+1)/(2*ell-3));
        y_ell = (x*y_mmp1*factor1 - (ell+m-1)*y_mm*factor2) / (ell-m);
        y_mm   = y_mmp1;
        y_mmp1 = y_ell;
	result_array[ell-m] = y_ell;
      }
    }

    return GSL_SUCCESS;
  }
}

#ifndef HIDE_INLINE_STATIC
int
gsl_sf_legendre_array_size(const int lmax, const int m)
{
  return lmax-m+1;
}
#endif

/*-*-*-*-*-*-*-*-*-* Functions w/ Natural Prototypes *-*-*-*-*-*-*-*-*-*-*/

#include "eval.h"

double gsl_sf_legendre_P1(const double x)
{
  EVAL_RESULT(gsl_sf_legendre_P1_e(x, &result));
}

double gsl_sf_legendre_P2(const double x)
{
  EVAL_RESULT(gsl_sf_legendre_P2_e(x, &result));
}

double gsl_sf_legendre_P3(const double x)
{
  EVAL_RESULT(gsl_sf_legendre_P3_e(x, &result));
}

double gsl_sf_legendre_Pl(const int l, const double x)
{
  EVAL_RESULT(gsl_sf_legendre_Pl_e(l, x, &result));
}

double gsl_sf_legendre_Plm(const int l, const int m, const double x)
{
  EVAL_RESULT(gsl_sf_legendre_Plm_e(l, m, x, &result));
}

double gsl_sf_legendre_sphPlm(const int l, const int m, const double x)
{
  EVAL_RESULT(gsl_sf_legendre_sphPlm_e(l, m, x, &result));
}

