/* specfunc/elljac.c
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
#include "gsl_sf_pow_int.h"
#include "gsl_sf_elljac.h"


/* See [Thompson, Atlas for Computing Mathematical Functions] */


int
gsl_sf_elljac_e(double u, double m, double * sn, double * cn, double * dn)
{
  if(fabs(m) > 1.0) {
    *sn = 0.0;
    *cn = 0.0;
    *dn = 0.0;
    GSL_ERROR ("|m| > 1.0", GSL_EDOM);
  }
  else if(fabs(m) < 2.0*GSL_DBL_EPSILON) {
    *sn = sin(u);
    *cn = cos(u);
    *dn = 1.0;
    return GSL_SUCCESS;
  }
  else if(fabs(m - 1.0) < 2.0*GSL_DBL_EPSILON) {
    *sn = tanh(u);
    *cn = 1.0/cosh(u);
    *dn = *cn;
    return GSL_SUCCESS;
  }
  else {
    int status = GSL_SUCCESS;
    const int N = 16;
    double   a[16];
    double   b[16];
    double   c[16];
    double phi[16];
    double psi[16]; /* psi[i] := phi[i] - Pi 2^{i-1} */
    double two_N;
    int n = 0;

    a[0] = 1.0;
    b[0] = sqrt(1.0 - m);
    c[0] = sqrt(m);

    while( fabs(c[n]) > 4.0 * GSL_DBL_EPSILON) {
      a[n+1] = 0.5 * (a[n] + b[n]);
      b[n+1] = sqrt(a[n] * b[n]);
      c[n+1] = 0.5 * (a[n] - b[n]);
      if(n >= N - 2) {
        status = GSL_EMAXITER;
	c[N-1] = 0.0;
	break;
      }
      ++n;
    }

    --n;
    two_N = (double)(1 << n ); /* 2^n */  /* gsl_sf_pow_int(2.0, n); */
    phi[n] = two_N * a[n] * u;
    psi[n] = two_N * (a[n]*u - 0.5*M_PI);

    while(n > 0) {
      const double psi_sgn = ( n == 1 ? -1.0 : 1.0 );
      const double phi_asin_arg = c[n] * sin(phi[n])/a[n];
      const double psi_asin_arg = c[n]/a[n] * psi_sgn * sin(psi[n]);
      const double phi_asin = asin(phi_asin_arg);
      const double psi_asin = asin(psi_asin_arg);
      phi[n-1] = 0.5 * (phi[n] + phi_asin);
      psi[n-1] = 0.5 * (psi[n] + psi_asin);
      --n;
    }

    *sn = sin(phi[0]);
    *cn = cos(phi[0]);
    {
      /* const double dn_method_1 = *cn / cos(phi[1] - phi[0]); */
      const double dn_method_2 = sin(psi[0])/sin(psi[1] - psi[0]);
      *dn = dn_method_2;
      /* printf("%18.16g  %18.16g\n", dn_method_1, dn_method_2); */
    }

    return status;
  }
}

