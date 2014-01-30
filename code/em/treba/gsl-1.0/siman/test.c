/* siman/test.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Mark Galassi
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
#include <string.h>
#include <math.h>
#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_siman.h>
#include <gsl/gsl_ieee_utils.h>
#include <stdio.h>

/* set up parameters for this simulated annealing run */
#define N_TRIES 200		/* how many points do we try before stepping */
#define ITERS_FIXED_T 1000	/* how many iterations for each T? */
#define STEP_SIZE 1.0		/* max step size in random walk */
#define K 1.0			/* Boltzmann constant */
#define T_INITIAL 0.008		/* initial temperature */
#define MU_T 1.003		/* damping factor for temperature */
#define T_MIN 2.0e-6

gsl_siman_params_t params = {N_TRIES, ITERS_FIXED_T, STEP_SIZE,
			     K, T_INITIAL, MU_T, T_MIN};

inline double square (double x) ;
inline double square (double x) { return x * x ; } 

double test_E_1D(Element x);
void test_step_1D(const gsl_rng * r, Element *x_p, gsl_siman_params_t p);
void print_pos_1D(Element x);
double distance_1D(Element x, Element y);

double test_E_1D(Element x)
{
  double val = x.D1;
  double u = val - 1 ;
  return exp(-u * u)*sin(8*val);
/* return sin(sin(val*val) - cos(val)) + cos(sin(val) + sin(val)*sin(val));*/
/* return 1.0/(square(x-1.0) + 1.0)*sin(8.0*x); */
/* return 1.0/(square(x-1.0) + 1.0)*sin(8.0*x) + x*x/100.0; */
}

/* takes a step for the test function; max distance: step_size.
 * the new point is put in x_p and returned.
 */
void test_step_1D(const gsl_rng * r, Element *x_p, gsl_siman_params_t p)
{
  double old_x = x_p->D1;
  double new_x;

  new_x = gsl_rng_uniform(r);
  new_x = new_x*2*p.step_size;
  new_x = new_x - p.step_size + old_x;

  x_p->D1 = new_x;
}


/* simple routine to print out a position value */
void print_pos_1D(Element x)
{
  printf("%12g", x.D1);
}

void debug_pos_1D(void * x)
{
  printf("%12g", ((Element *)x)->D1);
}


/* a metric for the 2D space */
double distance_1D(Element x, Element y)
{
  return fabs(y.D1 - x.D1);
}


double test_E_2D(Element x);
void test_step_2D(const gsl_rng * r, Element *x_p, gsl_siman_params_t p);
void print_pos_2D(Element x);
double distance_2D(Element x, Element y);

/* a 2-D function to be minimized */
double test_E_2D(Element x)
{
  double old_x = x.D2[0], old_y = x.D2[1];
  double u = old_x-1 ;
  double v = old_y - 0.8 ;
  return exp(-u * u - v * v)*sin(8*old_x + 8 * old_y);
}

/* takes a step for the test function; max distance: step_size.  the
   new point is put in x_p and returned. */
void test_step_2D(const gsl_rng * r, Element *x_p, gsl_siman_params_t p)
{
  double old_x =  x_p->D2[0], old_y = x_p->D2[1], new_x, new_y;

  new_x = gsl_rng_uniform(r);
  new_x = new_x*2*p.step_size;
  new_x = new_x - p.step_size + old_x;

  new_y = gsl_rng_uniform(r);
  new_y = new_y*2*p.step_size;
  new_y = new_y - p.step_size + old_y;

  x_p->D2[0] = new_x;
  x_p->D2[1] = new_y;
}

/* simple routine to print out a position value */
void print_pos_2D(Element x)
{
  printf("%g::%g", x.D2[0], x.D2[1]);
}

/* a metric for the 2D space */
double distance_2D(Element x, Element y)
{
  double u = y.D2[0]-x.D2[0] ;
  double v = y.D2[1]-x.D2[1] ;
  return sqrt(u*u + v*v);
}


/**********************************************/
/************ 3-dimensional search ************/
/**********************************************/

double test_E_3D(Element x);
void test_step_3D(const gsl_rng * r, Element *x_p, gsl_siman_params_t p);
void print_pos_3D(Element x);
double distance_3D(Element x, Element y);

/* a 3-D function to be minimized */
double test_E_3D(Element x)
{
  return exp(-square(x.D3[0]-1) - square(x.D3[1] - 0.8)
    - square(x.D3[2] - 0.8)) * sin(8*x.D3[0] + 8*x.D3[1] + 8*x.D3[2])
      + (square(x.D3[0]) + square(x.D3[1]) + square(x.D3[2]))/10000.0;
}

/* takes a step for the test function; max distance: step_size.
 * the new point is put in x_p and returned.
 */
void test_step_3D(const gsl_rng * r, Element *x_p, gsl_siman_params_t p)
{
  double old_x = x_p->D3[0], old_y = x_p->D3[1], old_z = x_p->D3[2];
  double new_x, new_y, new_z;

  new_x = gsl_rng_uniform(r);
  new_x = new_x*2*p.step_size;
  new_x = new_x - p.step_size + old_x;

  new_y = gsl_rng_uniform(r);
  new_y = new_y*2*p.step_size;
  new_y = new_y - p.step_size + old_y;

  new_z = gsl_rng_uniform(r);
  new_z = new_z*2*p.step_size;
  new_z = new_z - p.step_size + old_z;

  x_p->D3[0] = new_x;
  x_p->D3[1] = new_y;
  x_p->D3[2] = new_z;
}

/* simple routine to print out a position value */
void print_pos_3D(Element x)
{
  printf("%g::%g::%g", x.D3[0], x.D3[1], x.D3[2]);
}

/* a metric for the 2D space */
double distance_3D(Element x, Element y)
{
  return sqrt(square(y.D3[0]-x.D3[0]) + square(y.D3[1]-x.D3[1])
	      + square(y.D3[2]-x.D3[2]));
}

double E1(void *xp);
double M1(void *xp, void *yp);
void S1(const gsl_rng * r, void *xp, double step_size);
void P1(void *xp);

/* now some functions to test in one dimension */
double E1(void *xp)
{
  double x = * ((double *) xp);

/*   return exp(-square(x-1))*sin(8*x); */
  return exp(-square(x-1))*sin(8*x) - exp(-square(x-1000))*0.89;
}

double M1(void *xp, void *yp)
{
  double x = *((double *) xp);
  double y = *((double *) yp);

  return fabs(x - y);
}

void S1(const gsl_rng * r, void *xp, double step_size)
{
  double old_x = *((double *) xp);
  double new_x;

  new_x = gsl_rng_uniform(r)*2*step_size - step_size + old_x;
/*   new_x = new_x*2*step_size; */
/*   new_x = new_x - step_size + old_x; */

  /* printf("test step from old_x = %g to new_x = %g\n", old_x, new_x); */

  memcpy(xp, &new_x, sizeof(new_x));
}

void P1(void *xp)
{
  printf(" %12g ", *((double *) xp));
}

int main(void)
{
  double x_min = 1.36312999455315182 ;
  double x ;

  gsl_rng * r = gsl_rng_alloc (gsl_rng_env_setup()) ;

  gsl_ieee_env_setup ();

  /* The function tested here has multiple mimima. 
     The global minimum is at    x = 1.36312999, (f = -0.87287)
     There is a local minimum at x = 0.60146196, (f = -0.84893) */

  gsl_ieee_env_setup ();

  x = -10.0 ;
  gsl_siman_solve(r, &x, E1, S1, M1, NULL, NULL, NULL, NULL,
		  sizeof(double), params);
  gsl_test_rel(x, x_min, 1e-3, "f(x)= exp(-(x-1)^2) sin(8x), x0=-10") ;

  x = +10.0 ;
  gsl_siman_solve(r, &x, E1, S1, M1, NULL, NULL, NULL, NULL,
		  sizeof(double), params);
  gsl_test_rel(x, x_min, 1e-3, "f(x)= exp(-(x-1)^2) sin(8x), x0=10") ;

  /* Start at the false minimum */

  x = +0.6 ; 
  gsl_siman_solve(r, &x, E1, S1, M1, NULL, NULL, NULL, NULL,
		  sizeof(double), params);
  gsl_test_rel(x, x_min, 1e-3, "f(x)= exp(-(x-1)^2) sin(8x), x0=0.6") ;

  x = +0.5 ; 
  gsl_siman_solve(r, &x, E1, S1, M1, NULL, NULL, NULL, NULL,
		  sizeof(double), params);
  gsl_test_rel(x, x_min, 1e-3, "f(x)= exp(-(x-1)^2) sin(8x), x0=0.5") ;

  x = +0.4 ; 
  gsl_siman_solve(r, &x, E1, S1, M1, NULL, NULL, NULL, NULL,
		  sizeof(double), params);
  gsl_test_rel(x, x_min, 1e-3, "f(x)= exp(-(x-1)^2) sin(8x), x0=0.4") ;

  exit (gsl_test_summary ());

#ifdef JUNK 
  x0.D1 = 12.0;
  printf("#one dimensional problem, x0 = %f\n", x0.D1);
  gsl_siman_Usolve(r, &x0, test_E_1D, test_step_1D, distance_1D,
		   print_pos_1D, params);


  x0.D2[0] = 12.0;
  x0.D2[1] = 5.5;
  printf("#two dimensional problem, (x0,y0) = (%f,%f)\n",
	 x0.D2[0], x0.D2[1]);
  gsl_siman_Usolve(r, &x0, test_E_2D, test_step_2D, distance_2D,
		   print_pos_2D, params); 

  x0.D3[0] = 12.2;
  x0.D3[1] = 5.5;
  x0.D3[2] = -15.5;
  printf("#three dimensional problem, (x0,y0,z0) = (%f,%f,%f)\n",
	 x0.D3[0], x0.D3[1], x0.D3[2]);
  gsl_siman_Usolve(r, &x0, test_E_3D, test_step_3D, distance_3D, 
		   print_pos_3D, params); 

  x0.D2[0] = 12.2;
  x0.D2[1] = 5.5;

  gsl_siman_solve(r, &x0, test_E_2D, test_step_2D, distance_2D, print_pos_2D, params);
  
  x0.D3[0] = 12.2;
  x0.D3[1] = 5.5;
  x0.D3[2] = -15.5;

  gsl_siman_solve(r, &x0, test_E_3D, test_step_3D, distance_3D, print_pos_3D, params);

  return 0;
#endif
}
