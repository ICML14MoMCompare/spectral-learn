/* specfunc/test_dilog.c
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
#include <gsl/gsl_test.h>
#include <gsl/gsl_sf.h>
#include "test_sf.h"


int test_dilog(void)
{
  gsl_sf_result r;
  gsl_sf_result r1, r2;
  int s = 0;

  /* real dilog */

  TEST_SF(s, gsl_sf_dilog_e, (-3.0, &r),   -1.9393754207667089531,     TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_dilog_e, (-0.5, &r),   -0.4484142069236462024,     TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_dilog_e, (-0.001, &r), -0.0009997501110486510834,  TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_dilog_e, (0.1, &r),     0.1026177910993911,        TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_dilog_e, (0.7, &r),     0.8893776242860387386,     TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_dilog_e, (1.0, &r),     1.6449340668482260,        TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_dilog_e, (1.5, &r),     2.3743952702724802007,     TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_dilog_e, (2.0, &r),     2.4674011002723397,        TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_dilog_e, ( 5.0, &r),    1.7837191612666306277,     TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_dilog_e, ( 11.0, &r),   0.3218540439999117111,     TEST_TOL1, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_dilog_e, (12.59, &r),   0.0010060918167266208634,  TEST_TOL3, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_dilog_e, (12.595, &r),  0.00003314826006436236810, TEST_TOL5, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_dilog_e, (13.0, &r),   -0.07806971248458575855,    TEST_TOL2, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_dilog_e, (20.0, &r),   -1.2479770861745251168,     TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_dilog_e, (150.0, &r),  -9.270042702348657270,      TEST_TOL0, GSL_SUCCESS);
  TEST_SF(s, gsl_sf_dilog_e, (1100.0, &r), -21.232504073931749553,     TEST_TOL0, GSL_SUCCESS);


  /* complex dilog */
  /* FIXME: probably need more tests here... 
   * also need to work on accuracy for r->1; need to
   * adjust the switch-over point I suppose.
   */

  TEST_SF_2(s, gsl_sf_complex_dilog_e, (1.00001, M_PI/2.0, &r1, &r2),
            -0.20562022409960237363, TEST_TOL1,
             0.91597344814458309320, TEST_TOL1,
	     GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_complex_dilog_e, (0.99999, M_PI/2.0, &r1, &r2),
            -0.20561329262779687646, TEST_TOL0,
             0.91595774018131512060, TEST_TOL0,
	     GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_complex_dilog_e, (0.991, M_PI/2.0, &r1, &r2),
            -0.20250384721077806127, TEST_TOL0,
             0.90888544355846447810, TEST_TOL0,
	     GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_complex_dilog_e, (0.98, M_PI/2.0, &r1, &r2),
            -0.19871638377785918403, TEST_TOL2,
             0.90020045882981847610, TEST_TOL2,
	     GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_complex_dilog_e, (0.95, M_PI/2.0, &r1, &r2),
            -0.18848636456893572091, TEST_TOL1,
             0.87633754133420277830, TEST_TOL1,
	     GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_complex_dilog_e, (0.8, M_PI/2.0, &r1, &r2),
            -0.13980800855429037810, TEST_TOL0,
             0.75310609092419884460, TEST_TOL0,
	     GSL_SUCCESS);


  TEST_SF_2(s, gsl_sf_complex_dilog_e, (0.5, M_PI/2.0, &r1, &r2),
            -0.05897507442156586346, TEST_TOL0,
             0.48722235829452235710, TEST_TOL0,
	     GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_complex_dilog_e, (0.01, M_PI/2.0, &r1, &r2),
            -0.000024999375027776215378, TEST_TOL3,
             0.009999888892888684820, TEST_TOL3,
	     GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_complex_dilog_e, (10.0, M_PI/2.0, &r1, &r2),
            -3.0596887943287347304, TEST_TOL0,
             3.7167814930680685900, TEST_TOL0,
	     GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_complex_dilog_e, (100.0, M_PI/2.0, &r1, &r2),
            -11.015004738293824854, TEST_TOL0,
             7.2437843013083534970, TEST_TOL0,
	     GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_complex_dilog_e, (0.99, M_PI/8.0, &r1, &r2),
            1.0571539648820244720, TEST_TOL0,
            0.7469145254610851318, TEST_TOL0,
	    GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_complex_dilog_e, (0.99, M_PI/64.0, &r1, &r2),
            1.5381800285902999666, TEST_TOL0,
            0.1825271634987756651, TEST_TOL0,
	    GSL_SUCCESS);

  TEST_SF_2(s, gsl_sf_complex_dilog_e, (0.9, 3.0*M_PI/4.0, &r1, &r2),
            -0.6062840301356530985, TEST_TOL1,
             0.4836632833122775721, TEST_TOL1,
	    GSL_SUCCESS);

  return s;
}
