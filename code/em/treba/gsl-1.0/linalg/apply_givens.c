/* linalg/apply_givens.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001 Gerard Jungman, Brian Gough
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

inline static void
apply_givens_qr (size_t M, size_t N, gsl_matrix * Q, gsl_matrix * R,
		 size_t i, size_t j, double c, double s)
{
  size_t k;

  /* Apply rotation to matrix Q,  Q' = Q G */

  for (k = 0; k < M; k++)
    {
      double qki = gsl_matrix_get (Q, k, i);
      double qkj = gsl_matrix_get (Q, k, j);
      gsl_matrix_set (Q, k, i, qki * c - qkj * s);
      gsl_matrix_set (Q, k, j, qki * s + qkj * c);
    }

  /* Apply rotation to matrix R, R' = G^T R (note: upper triangular so
     zero for column < row) */

  for (k = GSL_MIN (i, j); k < N; k++)
    {
      double rik = gsl_matrix_get (R, i, k);
      double rjk = gsl_matrix_get (R, j, k);
      gsl_matrix_set (R, i, k, c * rik - s * rjk);
      gsl_matrix_set (R, j, k, s * rik + c * rjk);
    }
}

inline static void
apply_givens_vec (gsl_vector * v, size_t i, size_t j, double c, double s)
{
  /* Apply rotation to vector v' = G^T v */

  double vi = gsl_vector_get (v, i);
  double vj = gsl_vector_get (v, j);
  gsl_vector_set (v, i, c * vi - s * vj);
  gsl_vector_set (v, j, s * vi + c * vj);
}
