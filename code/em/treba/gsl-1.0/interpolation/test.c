/* interpolation/test.c
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
#include <config.h>
#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>


int
test_bsearch(void)
{
  double x_array[5] = { 0.0, 1.0, 2.0, 3.0, 4.0 };
  size_t index_result;
  int status = 0;
  int s;

  /* check an interior point */
  index_result = gsl_interp_bsearch(x_array, 1.5, 0, 4);
  s = (index_result != 1);
  status += s;
  gsl_test (s, "simple bsearch");

  /* check that we get the last interval if x == last value */
  index_result = gsl_interp_bsearch(x_array, 4.0, 0, 4);
  s = (index_result != 3);
  status += s;
  gsl_test (s, "upper endpoint bsearch");

  /* check that we get the first interval if x == first value */
  index_result = gsl_interp_bsearch(x_array, 0.0, 0, 4);
  s = (index_result != 0);
  status += s;
  gsl_test (s, "lower endpoint bsearch");

  /* check that we get correct interior boundary behaviour */
  index_result = gsl_interp_bsearch(x_array, 2.0, 0, 4);
  s = (index_result != 2);
  status += s;
  gsl_test (s, "degenerate bsearch");

  /* check out of bounds above */
  index_result = gsl_interp_bsearch(x_array, 10.0, 0, 4);
  s = (index_result != 3);
  status += s;
  gsl_test (s, "out of bounds bsearch +");

  /* check out of bounds below */
  index_result = gsl_interp_bsearch(x_array, -10.0, 0, 4);
  s = (index_result != 0);
  status += s;
  gsl_test (s, "out of bounds bsearch -");

  return status;
}



typedef double TEST_FUNC (double);
typedef struct _xy_table xy_table;

struct _xy_table
  {
    double * x;
    double * y;
    size_t n;
  };

xy_table make_xy_table (double x[], double y[], size_t n);

xy_table
make_xy_table (double x[], double y[], size_t n)
{
  xy_table t;
  t.x = x;
  t.y = y;
  t.n = n;
  return t;
}

static int
test_interp (
  const xy_table * data_table,
  const gsl_interp_type * T,
  xy_table * test_table,
  xy_table * test_d_table,
  xy_table * test_i_table
  )
{
  int status = 0;
  size_t i;

  gsl_interp_accel *a = gsl_interp_accel_alloc ();
  gsl_interp *interp = gsl_interp_alloc (T, data_table->n);

  gsl_interp_init (interp, data_table->x, data_table->y, data_table->n);

  for (i = 0; i < test_table->n; i++)
    {
      double x = test_table->x[i];
      double y;
      double deriv;
      double integ;
      double diff_y, diff_deriv, diff_integ;
      gsl_interp_eval_e (interp, data_table->x, data_table->y, x, a, &y);
      gsl_interp_eval_deriv_e (interp, data_table->x, data_table->y, x, a, &deriv);
      gsl_interp_eval_integ_e (interp, data_table->x, data_table->y, 0.0, x, a, &integ);
      diff_y = y - test_table->y[i];
      diff_deriv = deriv - test_d_table->y[i];
      diff_integ = integ - test_i_table->y[i];
      if (fabs (diff_y) > 1.e-10 || fabs(diff_deriv) > 1.0e-10 || fabs(diff_integ) > 1.0e-10) {
	status++;
      }
    }

  gsl_interp_accel_free (a);
  gsl_interp_free (interp);

  return status;
}

static int
test_linear (void)
{
  int s;

  double data_x[4] = { 0.0, 1.0, 2.0, 3.0 };
  double data_y[4] = { 0.0, 1.0, 2.0, 3.0 };
  double test_x[6] = { 0.0, 0.5, 1.0, 1.5, 2.5, 3.0 };
  double test_y[6] = { 0.0, 0.5, 1.0, 1.5, 2.5, 3.0 };
  double test_dy[6] = { 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
  double test_iy[6] = { 0.0, 0.125, 0.5, 9.0/8.0, 25.0/8.0, 9.0/2.0 };

  xy_table data_table = make_xy_table(data_x, data_y, 4);
  xy_table test_table = make_xy_table(test_x, test_y, 6);
  xy_table test_d_table = make_xy_table(test_x, test_dy, 6);
  xy_table test_i_table = make_xy_table(test_x, test_iy, 6);

  s = test_interp (&data_table, gsl_interp_linear, &test_table, &test_d_table, &test_i_table);
  gsl_test (s, "linear interpolation");
  return s;
}


static int
test_cspline (void)
{
  int s;

  double data_x[3] = { 0.0, 1.0, 2.0 };
  double data_y[3] = { 0.0, 1.0, 2.0 };
  double test_x[4] = { 0.0, 0.5, 1.0, 2.0 };
  double test_y[4] = { 0.0, 0.5, 1.0, 2.0 };
  double test_dy[4] = { 1.0, 1.0, 1.0, 1.0 };
  double test_iy[4] = { 0.0, 0.125, 0.5, 2.0 };

  xy_table data_table = make_xy_table(data_x, data_y, 3);
  xy_table test_table = make_xy_table(test_x, test_y, 4);
  xy_table test_d_table = make_xy_table(test_x, test_dy, 4);
  xy_table test_i_table = make_xy_table(test_x, test_iy, 4);

  s = test_interp (&data_table, gsl_interp_cspline, &test_table, &test_d_table, &test_i_table);
  gsl_test (s, "cspline interpolation");
  return s;
}


static int
test_akima (void)
{
  int s;

  double data_x[5] = { 0.0, 1.0, 2.0, 3.0, 4.0 };
  double data_y[5] = { 0.0, 1.0, 2.0, 3.0, 4.0 };
  double test_x[4] = { 0.0, 0.5, 1.0, 2.0 };
  double test_y[4] = { 0.0, 0.5, 1.0, 2.0 };
  double test_dy[4] = { 1.0, 1.0, 1.0, 1.0 };
  double test_iy[4] = { 0.0, 0.125, 0.5, 2.0 };

  xy_table data_table = make_xy_table(data_x, data_y, 5);
  xy_table test_table = make_xy_table(test_x, test_y, 4);
  xy_table test_d_table = make_xy_table(test_x, test_dy, 4);
  xy_table test_i_table = make_xy_table(test_x, test_iy, 4);

  s = test_interp (&data_table, gsl_interp_akima, &test_table, &test_d_table, &test_i_table);
  gsl_test (s, "akima interpolation");
  return s;
}


int 
main (int argc, char **argv)
{
  int status = 0;

  argc = 0;    /* prevent warnings about unused parameters */
  argv = 0;

  status += test_bsearch();
  status += test_linear();
  status += test_cspline();
  status += test_akima();

  exit (gsl_test_summary());
}
