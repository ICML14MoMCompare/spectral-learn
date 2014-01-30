/* interpolation/integ_eval_macro.h
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

/* function for doing the spline integral evaluation
   which is common to both the cspline and akima methods
 */

static inline double
integ_eval (double ai, double bi, double ci, double di, double xi, double a,
	    double b)
{
  const double t0 = b + a;
  const double t1 = a * a + a * b + b * b;
  const double t2 = a * a * a + a * a * b + b * b * a + b * b * b;
  const double bterm = 0.5 * bi * (t0 - 2.0 * xi);
  const double cterm = ci / 3.0 * (t1 - 3.0 * xi * (t0 - xi));
  const double dterm =
    di / 4.0 * (t2 - 2.0 * xi * (2.0 * t1 - xi * (3.0 * t0 - 2.0 * xi)));
  return (b - a) * (ai + bterm + cterm + dterm);
}
