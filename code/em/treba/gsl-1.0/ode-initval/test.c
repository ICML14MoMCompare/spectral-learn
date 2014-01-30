/* ode-initval/test_odeiv.c
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
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_ieee_utils.h>
#include "gsl_odeiv.h"

int rhs_linear (double t, const double y[], double f[], void *params);
int jac_linear (double t, const double y[], double *dfdy, double dfdt[],
                void *params);
int rhs_sin (double t, const double y[], double f[], void *params);
int jac_sin (double t, const double y[], double *dfdy, double dfdt[],
             void *params);
int rhs_exp (double t, const double y[], double f[], void *params);
int jac_exp (double t, const double y[], double *dfdy, double dfdt[],
             void *params);
int rhs_stiff (double t, const double y[], double f[], void *params);
int jac_stiff (double t, const double y[], double *dfdy, double dfdt[],
               void *params);
void test_stepper_linear (const gsl_odeiv_step_type * T, double h,
                         double base_prec);
void test_stepper_sin (const gsl_odeiv_step_type * T, double h, double base_prec);
void test_stepper_exp (const gsl_odeiv_step_type * T, double h, double base_prec);
void test_stepper_stiff (const gsl_odeiv_step_type * T, double h, double base_prec);
void test_evolve_system_flat (gsl_odeiv_step * step,
                              const gsl_odeiv_system * sys,
                              double t0, double t1, double hstart,
                              double y[], double yfin[],
                              double err_target, const char *desc);

void test_evolve_system (const gsl_odeiv_step_type * T,
                         const gsl_odeiv_system * sys,
                         double t0, double t1, double hstart,
                         double y[], double yfin[],
                         double err_target, const char *desc);

void test_evolve_exp (const gsl_odeiv_step_type * T, double h, double err);
void test_evolve_sin (const gsl_odeiv_step_type * T, double h, double err);
void test_evolve_stiff1 (const gsl_odeiv_step_type * T, double h, double err);
void test_evolve_stiff5 (const gsl_odeiv_step_type * T, double h, double err);

/* RHS for a + b t */

int
rhs_linear (double t, const double y[], double f[], void *params)
{
  f[0] = 0.0;
  f[1] = y[0];
  return GSL_SUCCESS;
}

int
jac_linear (double t, const double y[], double *dfdy, double dfdt[],
	    void *params)
{
  gsl_matrix dfdy_mat;
  dfdy_mat.size1 = 2;
  dfdy_mat.size2 = 2;
  dfdy_mat.tda = 2;
  dfdy_mat.data = dfdy;
  dfdy_mat.block = 0;
  gsl_matrix_set (&dfdy_mat, 0, 0, 0.0);
  gsl_matrix_set (&dfdy_mat, 0, 1, 0.0);
  gsl_matrix_set (&dfdy_mat, 1, 0, 1.0);
  gsl_matrix_set (&dfdy_mat, 1, 1, 0.0);
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  return GSL_SUCCESS;
}

gsl_odeiv_system rhs_func_lin = {
  rhs_linear,
  jac_linear,
  2,
  0
};


/* RHS for sin(t),cos(t) */

int
rhs_sin (double t, const double y[], double f[], void *params)
{
  f[0] = -y[1];
  f[1] = y[0];
  return GSL_SUCCESS;
}

int
jac_sin (double t, const double y[], double *dfdy, double dfdt[],
	 void *params)
{
  gsl_matrix dfdy_mat;
  dfdy_mat.data = dfdy;
  dfdy_mat.size1 = 2;
  dfdy_mat.size2 = 2;
  dfdy_mat.tda = 2;
  dfdy_mat.block = 0;
  gsl_matrix_set (&dfdy_mat, 0, 0, 0.0);
  gsl_matrix_set (&dfdy_mat, 0, 1, -1.0);
  gsl_matrix_set (&dfdy_mat, 1, 0, 1.0);
  gsl_matrix_set (&dfdy_mat, 1, 1, 0.0);
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  return GSL_SUCCESS;
}

gsl_odeiv_system rhs_func_sin = {
  rhs_sin,
  jac_sin,
  2,
  0
};


/* RHS for a exp(t)+ b exp(-t) */

int
rhs_exp (double t, const double y[], double f[], void *params)
{
  f[0] = y[1];
  f[1] = y[0];
  return GSL_SUCCESS;
}

int
jac_exp (double t, const double y[], double *dfdy, double dfdt[],
	 void *params)
{
  gsl_matrix dfdy_mat;
  dfdy_mat.data = dfdy;
  dfdy_mat.size1 = 2;
  dfdy_mat.size2 = 2;
  dfdy_mat.tda = 2;
  dfdy_mat.block = 0;
  gsl_matrix_set (&dfdy_mat, 0, 0, 0.0);
  gsl_matrix_set (&dfdy_mat, 0, 1, 1.0);
  gsl_matrix_set (&dfdy_mat, 1, 0, 1.0);
  gsl_matrix_set (&dfdy_mat, 1, 1, 0.0);
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  return GSL_SUCCESS;
}

gsl_odeiv_system rhs_func_exp = {
  rhs_exp,
  jac_exp,
  2,
  0
};


/* RHS for stiff example */

int
rhs_stiff (double t, const double y[], double f[], void *params)
{
  f[0] = 998.0 * y[0] + 1998.0 * y[1];
  f[1] = -999.0 * y[0] - 1999.0 * y[1];
  return GSL_SUCCESS;
}

int
jac_stiff (double t, const double y[], double *dfdy, double dfdt[],
	   void *params)
{
  gsl_matrix dfdy_mat;
  dfdy_mat.data = dfdy;
  dfdy_mat.size1 = 2;
  dfdy_mat.size2 = 2;
  dfdy_mat.tda = 2;
  dfdy_mat.block = 0;
  gsl_matrix_set (&dfdy_mat, 0, 0, 998.0);
  gsl_matrix_set (&dfdy_mat, 0, 1, 1998.0);
  gsl_matrix_set (&dfdy_mat, 1, 0, -999.0);
  gsl_matrix_set (&dfdy_mat, 1, 1, -1999.0);
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  return GSL_SUCCESS;
}

gsl_odeiv_system rhs_func_stiff = {
  rhs_stiff,
  jac_stiff,
  2,
  0
};


void
test_stepper_linear (const gsl_odeiv_step_type * T, double h,
		     double base_prec)
{
  int s = 0;
  double y[2];
  double yerr[2];
  double t;
  double del;
  double delmax = 0.0;
  int count = 0;

  gsl_odeiv_step *stepper = gsl_odeiv_step_alloc (T, 2);

  y[0] = 1.0;
  y[1] = 0.0;

  for (t = 0.0; t < 4.0; t += h)
    {
      gsl_odeiv_step_apply (stepper, t, h, y, yerr, 0, 0, &rhs_func_lin);
      del = fabs ((y[1] - (t + h)) / y[1]);
      delmax = GSL_MAX_DBL (del, delmax);
      if (del > (count + 1.0) * base_prec)
	{
	  printf ("  LINEAR(%20.17g)  %20.17g  %20.17g  %8.4g\n", t + h, y[1],
		  t + h, del);
	  s++;
	}
      count++;
    }

  gsl_test (s, "%s, linear [0,4], max relative error = %g",
	    gsl_odeiv_step_name (stepper), delmax);

  gsl_odeiv_step_free (stepper);
}


void
test_stepper_sin (const gsl_odeiv_step_type * T, double h, double base_prec)
{
  int s = 0;
  double y[2];
  double yerr[2];
  double t;
  double del;
  double delmax = 0.0;
  int count = 0;

  gsl_odeiv_step *stepper = gsl_odeiv_step_alloc (T, 2);

  y[0] = 1.0;
  y[1] = 0.0;

  for (t = 0.0; t < M_PI; t += h)
    {
      int stat;
      double sin_th = sin (t + h);
      gsl_odeiv_step_apply (stepper, t, h, y, yerr, 0, 0, &rhs_func_sin);
      del = fabs ((y[1] - sin_th) / sin_th);
      delmax = GSL_MAX_DBL (del, delmax);
      {
	if (t < 0.5 * M_PI)
	  {
	    stat = (del > (count + 1.0) * base_prec);
	  }
	else if (t < 0.7 * M_PI)
	  {
	    stat = (del > 1.0e+04 * base_prec);
	  }
	else if (t < 0.9 * M_PI)
	  {
	    stat = (del > 1.0e+06 * base_prec);
	  }
	else
	  {
	    stat = (del > 1.0e+09 * base_prec);
	  }
	if (stat != 0)
	  {
	    printf ("  SIN(%22.18g)  %22.18g  %22.18g  %10.6g\n", t + h, y[1],
		    sin_th, del);
	  }
	s += stat;
      }
      count++;
    }
  if (delmax > 1.0e+09 * base_prec)
    {
      s++;
      printf ("  SIN(0 .. M_PI)  delmax = %g\n", delmax);
    }

  gsl_test (s, "%s, sine [0,pi], max relative error = %g",
	    gsl_odeiv_step_name (stepper), delmax);

  delmax = 0.0;
  for (; t < 100.5 * M_PI; t += h)
    {
      gsl_odeiv_step_apply (stepper, t, h, y, yerr, 0, 0, &rhs_func_sin);
      del = fabs (y[1] - sin (t));
      delmax = GSL_MAX_DBL (del, delmax);
      count++;
    }

  if (del > count * 2.0 * base_prec)
    {
      s++;
      printf ("  SIN(%22.18g)  %22.18g  %22.18g  %10.6g\n", t + h, y[1],
	      sin (t), del);
    }

  gsl_test (s, "%s, sine [pi,100.5*pi], max absolute error = %g",
	    gsl_odeiv_step_name (stepper), delmax);

  gsl_odeiv_step_free (stepper);
}


void
test_stepper_exp (const gsl_odeiv_step_type * T, double h, double base_prec)
{
  int s = 0;
  double y[2];
  double yerr[2];
  double t;
  double del, delmax = 0.0;
  int count = 0;

  gsl_odeiv_step *stepper = gsl_odeiv_step_alloc (T, 2);

  y[0] = 1.0;
  y[1] = 1.0;

  for (t = 0.0; t < 20.0; t += h)
    {
      double ex = exp (t + h);
      gsl_odeiv_step_apply (stepper, t, h, y, yerr, 0, 0, &rhs_func_exp);
      del = fabs ((y[1] - ex) / y[1]);
      delmax = GSL_MAX_DBL (del, delmax);
      if (del > (count + 1.0) * 2.0 * base_prec)
	{
	  printf ("  EXP(%20.17g)  %20.17g  %20.17g  %8.4g\n", t + h, y[1],
		  ex, del);
	  s++;
	}
      count++;
    }

  gsl_test (s, "%s, exponential [0,20], max relative error = %g",
	    gsl_odeiv_step_name (stepper), delmax);

  gsl_odeiv_step_free (stepper);
}

void
test_stepper_stiff (const gsl_odeiv_step_type * T, double h, double base_prec)
{
  int s = 0;
  double y[2];
  double yerr[2];
  double t;
  double del, delmax = 0.0;
  int count = 0;

  gsl_odeiv_step *stepper = gsl_odeiv_step_alloc (T, 2);

  y[0] = 1.0;
  y[1] = 0.0;

  for (t = 0.0; t < 20.0; t += h)
    {
      gsl_odeiv_step_apply (stepper, t, h, y, yerr, NULL, NULL,
			    &rhs_func_stiff);

      if (t > 0.04)
	{
	  double arg = t + h;
	  double e1 = exp (-arg);
	  double e2 = exp (-1000.0 * arg);
	  double u = 2.0 * e1 - e2;
	  /* double v = -e1 + e2; */
	  del = fabs ((y[0] - u) / y[0]);
	  delmax = GSL_MAX_DBL (del, delmax);

	  if (del > (count + 1.0) * 100.0 * base_prec)
	    {
	      printf ("  STIFF(%20.17g)  %20.17g  %20.17g  %8.4g\n", arg,
		      y[0], u, del);
	      s++;
	    }
	}
      count++;
    }

  gsl_test (s, "%s, stiff [0,20], max relative error = %g",
	    gsl_odeiv_step_name (stepper), delmax);

  gsl_odeiv_step_free (stepper);
}

void
test_evolve_system_flat (gsl_odeiv_step * step,
			 const gsl_odeiv_system * sys,
			 double t0, double t1, double hstart,
			 double y[], double yfin[],
			 double err_target, const char *desc)
{
  int s = 0;
  double frac;

  double t = t0;
  double h = hstart;

  gsl_odeiv_evolve *e = gsl_odeiv_evolve_alloc (sys->dimension);

  while (t < t1)
    {
      gsl_odeiv_evolve_apply (e, NULL, step, sys, &t, t1, &h, y);
    }

  frac = fabs ((y[1] - yfin[1]) / yfin[1]) + fabs ((y[0] - yfin[0]) / yfin[0]);

  if (frac > 2.0 * e->count * err_target)
    {
      printf ("FLAT t = %.5e  y0 = %g y1= %g\n", t, y[0], y[1]);
      s++;
    }

  gsl_test (s, "%s, %s, evolve, no control, max relative error = %g",
	    gsl_odeiv_step_name (step), desc, frac);

  gsl_odeiv_evolve_free (e);
}


void
test_evolve_system (const gsl_odeiv_step_type * T,
		    const gsl_odeiv_system * sys,
		    double t0, double t1, double hstart,
		    double y[], double yfin[],
		    double err_target, const char *desc)
{
  int s = 0;
  double frac;

  double t = t0;
  double h = hstart;

  gsl_odeiv_step * step = gsl_odeiv_step_alloc (T, sys->dimension);

  gsl_odeiv_control *c =
    gsl_odeiv_control_standard_new (0.0, err_target, 1.0, 1.0);
  gsl_odeiv_evolve *e = gsl_odeiv_evolve_alloc (sys->dimension);

  while (t < t1)
    {
      gsl_odeiv_evolve_apply (e, c, step, sys, &t, t1, &h, y);
      /* printf ("SYS  t = %.18e h = %g y0 = %g y1= %g\n", t, h, y[0], y[1]); */
    }

  frac = fabs ((y[1] - yfin[1]) / yfin[1]) + fabs ((y[0] - yfin[0]) / yfin[0]);

  if (frac > 2.0 * e->count * err_target)
    {
      printf ("SYS  t = %.5e h = %g y0 = %g y1= %g\n", t, h, y[0], y[1]);
      s++;
    }

  gsl_test (s, "%s, %s, evolve, standard control, relative error = %g",
	    gsl_odeiv_step_name (step), desc, frac);

  gsl_odeiv_evolve_free (e);
  gsl_odeiv_control_free (c);
  gsl_odeiv_step_free (step);
}


void
test_evolve_exp (const gsl_odeiv_step_type * T, double h, double err)
{
  double y[2];
  double yfin[2];

  y[0] = 1.0;
  y[1] = 1.0;
  yfin[0] = exp (10.0);
  yfin[1] = yfin[0];
  test_evolve_system (T, &rhs_func_exp, 0.0, 10.0, h, y, yfin, err, "exp [0,10]");
}

void
test_evolve_sin (const gsl_odeiv_step_type * T, double h, double err)
{
  double y[2];
  double yfin[2];

  y[0] = 1.0;
  y[1] = 0.0;
  yfin[0] = cos (2.0);
  yfin[1] = sin (2.0);
  test_evolve_system (T, &rhs_func_sin, 0.0, 2.0, h, y, yfin, err, "sine [0,2]");
}

void
test_evolve_stiff1 (const gsl_odeiv_step_type * T, double h, double err)
{
  double y[2];
  double yfin[2];

  y[0] = 1.0;
  y[1] = 0.0;
  {
    double arg = 1.0;
    double e1 = exp (-arg);
    double e2 = exp (-1000.0 * arg);
    yfin[0] = 2.0 * e1 - e2;
    yfin[1] = -e1 + e2;
  }
  test_evolve_system (T, &rhs_func_stiff, 0.0, 1.0, h, y, yfin, err,
                      "stiff [0,1]");
}

void
test_evolve_stiff5 (const gsl_odeiv_step_type * T, double h, double err)
{
  double y[2];
  double yfin[2];

  y[0] = 1.0;
  y[1] = 0.0;
  {
    double arg = 5.0;
    double e1 = exp (-arg);
    double e2 = exp (-1000.0 * arg);
    yfin[0] = 2.0 * e1 - e2;
    yfin[1] = -e1 + e2;
  }
  test_evolve_system (T, &rhs_func_stiff, 0.0, 5.0, h, y, yfin, err,
                      "stiff [0,5]");
}



int
main (void)
{
  int i;

  struct ptype
  {
    const gsl_odeiv_step_type *type;
    double h;
  }
  p[20];

  p[0].type = gsl_odeiv_step_rk2;
  p[0].h = 1.0e-03;
  p[1].type = gsl_odeiv_step_rk2imp;
  p[1].h = 1.0e-03;
  p[2].type = gsl_odeiv_step_rk4;
  p[2].h = 1.0e-03;
  p[3].type = gsl_odeiv_step_rk4imp;
  p[3].h = 1.0e-03;
  p[4].type = gsl_odeiv_step_rkf45;
  p[4].h = 1.0e-03;
  p[5].type = gsl_odeiv_step_rk8pd;
  p[5].h = 1.0e-03;
  p[6].type = gsl_odeiv_step_rkck;
  p[6].h = 1.0e-03;
  p[7].type = gsl_odeiv_step_bsimp;
  p[7].h = 0.1;
  p[8].type = gsl_odeiv_step_gear1;
  p[8].h = 1.0e-03;
  p[9].type = gsl_odeiv_step_gear2;
  p[9].h = 1.0e-03;
  p[10].type = 0;

  gsl_ieee_env_setup ();

  for (i = 0; p[i].type != 0; i++)
    {
      test_stepper_linear (p[i].type, p[i].h, GSL_SQRT_DBL_EPSILON);
      test_stepper_sin (p[i].type, p[i].h / 10.0, GSL_SQRT_DBL_EPSILON);
      test_stepper_exp (p[i].type, p[i].h / 10.0, GSL_SQRT_DBL_EPSILON);
      test_stepper_stiff (p[i].type, p[i].h / 10.0, GSL_SQRT_DBL_EPSILON);
    }

  for (i = 0; p[i].type != 0; i++)
    {
      test_evolve_exp (p[i].type, p[i].h, GSL_SQRT_DBL_EPSILON);
      test_evolve_sin (p[i].type, p[i].h, GSL_SQRT_DBL_EPSILON);
      test_evolve_stiff1 (p[i].type, p[i].h, 1e-5);
      test_evolve_stiff5 (p[i].type, p[i].h, 1e-5);
    }

  exit (gsl_test_summary ());
}
