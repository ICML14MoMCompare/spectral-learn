/* randist/exppow.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 James Theiler, Brian Gough
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
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/* The exponential power probability distribution is  

   p(x) dx = (1/(2 a Gamma(1+1/b))) * exp(-|x/a|^b) dx

   for -infty < x < infty. For b = 1 it reduces to the Laplace
   distribution. 

   The exponential power distribution is related to the gamma
   distribution by E = a * pow(G(1/b),1/b), where E is an exponential
   power variate and G is a gamma variate.

   We use this relation for b < 1. For b >=1 we use rejection methods
   based on the laplace and gaussian distributions which should be
   faster.

   See P. R. Tadikamalla, "Random Sampling from the Exponential Power
   Distribution", Journal of the American Statistical Association,
   September 1980, Volume 75, Number 371, pages 683-686.
   
*/

double
gsl_ran_exppow (const gsl_rng * r, const double a, const double b)
{
  if (b < 1) 
    {
      double u = gsl_rng_uniform (r) ;
      double v = gsl_ran_gamma (r, 1/b, 1.0) ;
      double z = a * pow(v, 1/b) ;

      if (u > 0.5) 
	{
	  return z ;
	} 
      else 
	{
	  return -z ;
	}
    }
  else if (b == 1) 
    {
      /* Laplace distribution */
      return gsl_ran_laplace (r, a) ;
    }
  else if (b < 2) 
    {
      /* Use laplace distribution for rejection method */

      double x, y, h, ratio, u ;

      /* Scale factor chosen by upper bound on ratio at b = 2 */

      double s = 1.4489 ; 
      do 
	{
	  x = gsl_ran_laplace (r, a) ;
	  y = gsl_ran_laplace_pdf (x,a) ;
	  h = gsl_ran_exppow_pdf (x,a,b) ;
	  ratio = h/(s * y) ;
	  u = gsl_rng_uniform (r) ;
	} 
      while (u > ratio) ;
      
      return x ;
    }
  else if (b == 2)
    {
      /* Gaussian distribution */
      return gsl_ran_gaussian (r, a/sqrt(2.0)) ;
    }
  else
    {
      /* Use gaussian for rejection method */

      double x, y, h, ratio, u ;
      const double sigma = a / sqrt(2.0) ;

      /* Scale factor chosen by upper bound on ratio at b = infinity.
	 This could be improved by using a rational function
	 approximation to the bounding curve. */

      double s = 2.4091 ;  /* this is sqrt(pi) e / 2 */

      do 
	{
	  x = gsl_ran_gaussian (r, sigma) ;
	  y = gsl_ran_gaussian_pdf (x, sigma) ;
	  h = gsl_ran_exppow_pdf (x, a, b) ;
	  ratio = h/(s * y) ;
	  u = gsl_rng_uniform (r) ;
	} 
      while (u > ratio) ;

      return x;
    }
}

double
gsl_ran_exppow_pdf (const double x, const double a, const double b)
{
  double p ;
  double lngamma = gsl_sf_lngamma (1+1/b) ;
  p = (1/(2*a)) * exp(-pow(fabs(x/a),b) - lngamma);
  return p;
}

