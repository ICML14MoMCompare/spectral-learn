/* eigen/jacobi.c
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

/* Author:  G. Jungman
 */
/* Simple linear algebra operations, operating directly
 * on the gsl_vector and gsl_matrix objects. These are
 * meant for "generic" and "small" systems. Anyone
 * interested in large systems will want to use more
 * sophisticated methods, presumably involving native
 * BLAS operations, specialized data representations,
 * or other optimizations.
 */
#include <config.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "gsl_eigen.h"

#define REAL double

inline static void
jac_rotate(gsl_matrix * a,
           unsigned int i, unsigned int j, unsigned int k, unsigned int l,
           double * g, double * h,
           double s, double tau)
{
  *g = gsl_matrix_get(a, i, j);
  *h = gsl_matrix_get(a, k, l);
  gsl_matrix_set(a, i, j, (*g) - s*((*h) + (*g)*tau));
  gsl_matrix_set(a, k, l, (*h) + s*((*g) - (*h)*tau));
}

int
gsl_eigen_jacobi(gsl_matrix * a,
                         gsl_vector * eval,
                         gsl_matrix * evec,
                         unsigned int max_rot, 
                         unsigned int * nrot)
{
  if(a->size1 != a->size2) {
     GSL_ERROR ("eigenproblem requires square matrix", GSL_ENOTSQR);
  }
  else if(a->size1 != evec->size1 || a->size1 != evec->size2) {
     GSL_ERROR ("eigenvector matrix must match input matrix", GSL_EBADLEN);
  }
  else if(a->size1 != eval->size) {
    GSL_ERROR ("eigenvalue vector must match input matrix", GSL_EBADLEN);
  }
  else {
    const unsigned int n = a->size1;
    unsigned int i, j, iq, ip;
    double t, s;

    REAL * b = (REAL *) malloc(n * sizeof(REAL));
    REAL * z = (REAL *) malloc(n * sizeof(REAL));
    if(b == 0 || z == 0) {
      if(b != 0) free(b);
      if(z != 0) free(z);
      GSL_ERROR ("could not allocate memory for workspace", GSL_ENOMEM);
    }

    /* Set eigenvectors to coordinate basis. */
    for(ip=0; ip<n; ip++) {
      for(iq=0; iq<n; iq++) {
        gsl_matrix_set(evec, ip, iq, 0.0);
      }
      gsl_matrix_set(evec, ip, ip, 1.0);
    }

    /* Initialize eigenvalues and workspace. */
    for(ip=0; ip<n; ip++) {
      REAL a_ipip = gsl_matrix_get(a, ip, ip);
      z[ip] = 0.0;
      b[ip] = a_ipip;
      gsl_vector_set(eval, ip, a_ipip);
    }

    *nrot = 0;

    for(i=1; i<=max_rot; i++) {
      REAL thresh;
      REAL tau;
      REAL g, h, c;
      REAL sm = 0.0;
      for(ip=0; ip<n-1; ip++) {
        for(iq=ip+1; iq<n; iq++) {
          sm += fabs(gsl_matrix_get(a, ip, iq));
        }
      }
      if(sm == 0.0) {
        free(z);
        free(b);
        return GSL_SUCCESS;
      }

      if(i < 4)
        thresh = 0.2*sm/(n*n);
      else
        thresh = 0.0;

      for(ip=0; ip<n-1; ip++) {
        for(iq=ip+1; iq<n; iq++) {
          const REAL d_ip = gsl_vector_get(eval, ip);
          const REAL d_iq = gsl_vector_get(eval, iq);
	  const REAL a_ipiq = gsl_matrix_get(a, ip, iq);
          g = 100.0 * fabs(a_ipiq);
          if(   i > 4
             && fabs(d_ip)+g == fabs(d_ip)
             && fabs(d_iq)+g == fabs(d_iq)
	     ) {
            gsl_matrix_set(a, ip, iq, 0.0);
          }
          else if(fabs(a_ipiq) > thresh) {
            h = d_iq - d_ip;
            if(fabs(h) + g == fabs(h)) {
              t = a_ipiq/h;
            }
            else {
              REAL theta = 0.5*h/a_ipiq;
              t = 1.0/(fabs(theta) + sqrt(1.0 + theta*theta));
              if(theta < 0.0) t = -t;
            }

            c   = 1.0/sqrt(1.0+t*t);
            s   = t*c;
            tau = s/(1.0+c);
            h   = t * a_ipiq;
            z[ip] -= h;
            z[iq] += h;
	    gsl_vector_set(eval, ip, d_ip - h);
	    gsl_vector_set(eval, iq, d_iq + h);
	    gsl_matrix_set(a, ip, iq, 0.0);

            for(j=0; j<ip; j++){
              jac_rotate(a, j, ip, j, iq, &g, &h, s, tau);
            }
            for(j=ip+1; j<iq; j++){
              jac_rotate(a, ip, j, j, iq, &g, &h, s, tau);
            }
            for(j=iq+1; j<n; j++){
              jac_rotate(a, ip, j, iq, j, &g, &h, s, tau);
            }
            for (j=0; j<n; j++){
              jac_rotate(evec, j, ip, j, iq, &g, &h, s, tau);
            }
            ++(*nrot);
          }
        }
      }
      for (ip=0; ip<n; ip++) {
        b[ip] += z[ip];
        z[ip]  = 0.0;
	gsl_vector_set(eval, ip, b[ip]);
      }

      /* continue iteration */
    }

    return GSL_EMAXITER;
  }
}

int
gsl_eigen_invert_jacobi(const gsl_matrix * a,
                             gsl_matrix * ainv,
                             unsigned int max_rot)
{
  if(a->size1 != a->size2 || ainv->size1 != ainv->size2) {
    return GSL_ENOTSQR;
  }
  else if(a->size1 != ainv->size2) {
    return GSL_EBADLEN;
  }
  else {
    const unsigned int n = a->size1;
    unsigned int nrot;
    unsigned int i,j,k,l;

    /* This is annoying because I do not want
     * the error handling in these functions.
     * But there are no "impl"-like versions
     * of these allocators... sigh.
     */
    gsl_vector * eval = gsl_vector_alloc(n);
    gsl_matrix * evec = gsl_matrix_alloc(n, n);
    gsl_matrix * inv_diag = gsl_matrix_alloc(n, n);

    if(eval == 0 || evec == 0 || inv_diag == 0) {
      if(eval != 0) gsl_vector_free(eval);
      if(evec != 0) gsl_matrix_free(evec);
      if(inv_diag != 0) gsl_matrix_free(inv_diag);
      return GSL_ENOMEM;
    }

    memcpy(ainv->data, a->data, n*n*sizeof(REAL));

    gsl_eigen_jacobi(ainv, eval, evec, max_rot, &nrot);

    for(i=0; i<n; i++) {
      if(fabs(gsl_vector_get(eval, i)) < 100.0 * GSL_DBL_EPSILON) {
        /* apparent singularity */
        gsl_vector_free(eval);
        gsl_matrix_free(evec);
        gsl_matrix_free(inv_diag);
        return GSL_ESING;
      }
    }

    /* Invert the diagonalized matrix. */
    for(i=0; i<n; i++) {
      for(j=0; j<n; j++) {
        gsl_matrix_set(inv_diag, i, j, 0.0);
      }
      gsl_matrix_set(inv_diag, i, i, 1.0/gsl_vector_get(eval, i));
    }

    for(i=0; i<n; i++) {
      for(j=0; j<n; j++) {
        gsl_matrix_set(ainv, i, j, 0.0);
        for(k=0; k<n; k++) {
          for(l=0; l<n; l++) {
	    REAL ainv_ij = gsl_matrix_get(ainv, i, j);
	    REAL evec_il = gsl_matrix_get(evec, i, l);
	    REAL evec_jk = gsl_matrix_get(evec, j, k);
	    REAL inv_diag_lk = gsl_matrix_get(inv_diag, l, k);
	    REAL delta = evec_il * inv_diag_lk * evec_jk;
	    gsl_matrix_set(ainv, i, j, ainv_ij + delta);
	  }
	}
      }
    }

    gsl_vector_free(eval);
    gsl_matrix_free(evec);
    gsl_matrix_free(inv_diag);
    return GSL_SUCCESS;
  }
}
