/* rng/benchmark.c
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
#include <time.h>
#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>

void benchmark (const gsl_rng_type * T);

#define N  1000000
int isum;
double dsum;

int
main (void)
{
  benchmark(gsl_rng_ranlux);
  benchmark(gsl_rng_ranlux389);
  benchmark(gsl_rng_ranlxs0);
  benchmark(gsl_rng_ranlxs1);
  benchmark(gsl_rng_ranlxs2);
  benchmark(gsl_rng_ranlxd1);
  benchmark(gsl_rng_ranlxd2);

  benchmark(gsl_rng_slatec);
  benchmark(gsl_rng_gfsr4);
  benchmark(gsl_rng_cmrg);
  benchmark(gsl_rng_minstd);
  benchmark(gsl_rng_mrg);
  benchmark(gsl_rng_mt19937);
  benchmark(gsl_rng_r250);
  benchmark(gsl_rng_ran0);
  benchmark(gsl_rng_ran1);
  benchmark(gsl_rng_ran2);
  benchmark(gsl_rng_ran3);
  benchmark(gsl_rng_rand48);
  benchmark(gsl_rng_rand);
  benchmark(gsl_rng_random8_bsd);
  benchmark(gsl_rng_random8_glibc2);
  benchmark(gsl_rng_random8_libc5);
  benchmark(gsl_rng_random128_bsd);
  benchmark(gsl_rng_random128_glibc2);
  benchmark(gsl_rng_random128_libc5);
  benchmark(gsl_rng_random256_bsd);
  benchmark(gsl_rng_random256_glibc2);
  benchmark(gsl_rng_random256_libc5);
  benchmark(gsl_rng_random32_bsd);
  benchmark(gsl_rng_random32_glibc2);
  benchmark(gsl_rng_random32_libc5);
  benchmark(gsl_rng_random64_bsd);
  benchmark(gsl_rng_random64_glibc2);
  benchmark(gsl_rng_random64_libc5);
  benchmark(gsl_rng_random_bsd);
  benchmark(gsl_rng_random_glibc2);
  benchmark(gsl_rng_random_libc5);
  benchmark(gsl_rng_randu);
  benchmark(gsl_rng_ranf);

  benchmark(gsl_rng_ranmar);
  benchmark(gsl_rng_taus);
  benchmark(gsl_rng_transputer);
  benchmark(gsl_rng_tt800);
  benchmark(gsl_rng_uni32);
  benchmark(gsl_rng_uni);
  benchmark(gsl_rng_vax);
  benchmark(gsl_rng_zuf);
  return 0;
}

void
benchmark (const gsl_rng_type * T)
{
  int start, end;
  int i = 0, d = 0 ;
  double t1, t2;

  gsl_rng *r = gsl_rng_alloc (T);

  start = clock ();
  do
    {
      int j;
      for (j = 0; j < N; j++)
	isum += gsl_rng_get (r);

      i += N;
      end = clock ();
    }
  while (end < start + CLOCKS_PER_SEC/10);

  t1 = (end - start) / (double) CLOCKS_PER_SEC;

  start = clock ();
  do
    {
      int j;
      for (j = 0; j < N; j++)
	dsum += gsl_rng_uniform (r);

      d += N;
      end = clock ();
    }
  while (end < start + CLOCKS_PER_SEC/10);

  t2 = (end - start) / (double) CLOCKS_PER_SEC;


  printf ("%6.0f k ints/sec, %6.0f k doubles/sec, %s\n",
	  i / t1 / 1000.0, d / t2 / 1000.0, gsl_rng_name (r));

  gsl_rng_free (r);
}
