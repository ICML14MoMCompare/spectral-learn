/* linalg/householder.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman, Brian Gough
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
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

#include "gsl_linalg.h"

double
gsl_linalg_householder_transform (gsl_vector * v)
{
  /* replace v[0:n-1] with a householder vector (v[0:n-1]) and
     coefficient tau that annihilate v[1:n-1] */

  const size_t n = v->size ;

  if (n == 1)
    {
      return 0.0; /* tau = 0 */
    }
  else
    { 
      double alpha, beta, tau ;
      
      gsl_vector_view x = gsl_vector_subvector (v, 1, n - 1) ; 
      
      double xnorm = gsl_blas_dnrm2 (&x.vector);
      
      if (xnorm == 0) 
        {
          return 0.0; /* tau = 0 */
        }
      
      alpha = gsl_vector_get (v, 0) ;
      beta = - (alpha >= 0.0 ? +1.0 : -1.0) * hypot(alpha, xnorm) ;
      tau = (beta - alpha) / beta ;
      
      gsl_blas_dscal (1.0 / (alpha - beta), &x.vector);
      gsl_vector_set (v, 0, beta) ;
      
      return tau;
    }
}

int
gsl_linalg_householder_hm (double tau, const gsl_vector * v, gsl_matrix * A)
{
  /* applies a householder transformation v,tau to matrix m */

  size_t i, j;

  if (tau == 0.0)
    {
      return GSL_SUCCESS;
    }

  for (j = 0; j < A->size2; j++)
    {
      /* Compute wj = Akj vk */

      double wj = gsl_matrix_get(A,0,j);  

      for (i = 1; i < A->size1; i++)  /* note, computed for v(0) = 1 above */
        {
          wj += gsl_matrix_get(A,i,j) * gsl_vector_get(v,i);
        }

      /* Aij = Aij - tau vi wj */

      /* i = 0 */
      {
        double A0j = gsl_matrix_get (A, 0, j);
        gsl_matrix_set (A, 0, j, A0j - tau *  wj);
      }

      /* i = 1 .. M-1 */

      for (i = 1; i < A->size1; i++)
        {
          double Aij = gsl_matrix_get (A, i, j);
          double vi = gsl_vector_get (v, i);
          gsl_matrix_set (A, i, j, Aij - tau * vi * wj);
        }
    }

  return GSL_SUCCESS;
}

int
gsl_linalg_householder_mh (double tau, const gsl_vector * v, gsl_matrix * A)
{
  /* applies a householder transformation v,tau to matrix m from the
     right hand side in order to zero out rows */

  size_t i, j;

  if (tau == 0)
    return GSL_SUCCESS;

  /* A = A - tau w v' */

  for (i = 0; i < A->size1; i++)
    {
      double wi = gsl_matrix_get(A,i,0);  

      for (j = 1; j < A->size2; j++)  /* note, computed for v(0) = 1 above */
        {
          wi += gsl_matrix_get(A,i,j) * gsl_vector_get(v,j);
        }
      
      /* j = 0 */
      
      {
        double Ai0 = gsl_matrix_get (A, i, 0);
        gsl_matrix_set (A, i, 0, Ai0 - tau *  wi);
      }

      /* j = 1 .. N-1 */
      
      for (j = 1; j < A->size2; j++) 
        {
          double vj = gsl_vector_get (v, j);
          double Aij = gsl_matrix_get (A, i, j);
          gsl_matrix_set (A, i, j, Aij - tau * wi * vj);
        }
    }

  return GSL_SUCCESS;
}

int
gsl_linalg_householder_hv (double tau, const gsl_vector * v, gsl_vector * w)
{
  /* applies a householder transformation v to vector w */
  size_t i;
  double d = 0;
 
  if (tau == 0)
    return GSL_SUCCESS ;

  /* d = v'w */

  d = gsl_vector_get(w,0);

  for (i = 1 ; i < v->size ; i++)
    {
      d += gsl_vector_get(v,i) * gsl_vector_get(w,i);
    }

  /* w = w - tau (v) (v'w) */
  
  {
    double w0 = gsl_vector_get (w,0);
    gsl_vector_set (w, 0, w0 - tau * d);
  }

  for (i = 1; i < v->size ; i++)
    {
      double wi = gsl_vector_get (w,i);
      double vi = gsl_vector_get (v,i);
      gsl_vector_set (w, i, wi - tau * vi * d);
    }

  return GSL_SUCCESS;
}


int
gsl_linalg_householder_hm1 (double tau, gsl_matrix * A)
{
  /* applies a householder transformation v,tau to a matrix being
     build up from the identity matrix, using the first column of A as
     a householder vector */

  size_t i, j;

  if (tau == 0)
    {
      gsl_matrix_set (A, 0, 0, 1.0);
      
      for (j = 1; j < A->size2; j++)
        {
          gsl_matrix_set (A, 0, j, 0.0);
        }

      for (i = 1; i < A->size1; i++)
        {
          gsl_matrix_set (A, i, 0, 0.0);
        }

      return GSL_SUCCESS;
    }

  /* w = A' v */

  for (j = 1; j < A->size2; j++)
    {
      double wj = 0.0;   /* A0j * v0 */

      for (i = 1; i < A->size1; i++)
        {
          double vi = gsl_matrix_get(A, i, 0);
          wj += gsl_matrix_get(A,i,j) * vi;
        }

      /* A = A - tau v w' */

      gsl_matrix_set (A, 0, j, - tau *  wj);
      
      for (i = 1; i < A->size1; i++)
        {
          double vi = gsl_matrix_get (A, i, 0);
          double Aij = gsl_matrix_get (A, i, j);
          gsl_matrix_set (A, i, j, Aij - tau * vi * wj);
        }
    }

  for (i = 1; i < A->size1; i++)
    {
      double vi = gsl_matrix_get(A, i, 0);
      gsl_matrix_set(A, i, 0, -tau * vi);
    }

  gsl_matrix_set (A, 0, 0, 1.0 - tau);

  return GSL_SUCCESS;
}
