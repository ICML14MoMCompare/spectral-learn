/* vector/test_complex_source.c
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

void FUNCTION (test, func) (void);
void FUNCTION (test, binary) (void);
void FUNCTION (test, trap) (void);

void
FUNCTION (test, func) (void)
{
  size_t i;

  TYPE (gsl_vector) * v = FUNCTION (gsl_vector, calloc) (N);

  gsl_test (v->data == 0, NAME (gsl_vector) "_alloc returns valid pointer");
  gsl_test (v->size != N, NAME (gsl_vector) "_alloc returns valid size");
  gsl_test (v->stride != 1, NAME (gsl_vector) "_alloc returns unit stride");

  for (i = 0; i < N; i++)
    {
      BASE x = ZERO;
      GSL_REAL (x) = (ATOMIC)i;
      GSL_IMAG (x) = (ATOMIC)(i + 1234);
      FUNCTION (gsl_vector, set) (v, i, x);
    };

  status = 0;

  for (i = 0; i < N; i++)
    {
      if (v->data[2 * i] != (ATOMIC) i || v->data[2 * i + 1] != (ATOMIC) (i + 1234))
	status = 1;
    };
  
  gsl_test (status, NAME (gsl_vector) "_set writes into array correctly");

  status = 0;

  for (i = 0; i < N; i++)
    {
      BASE x, y;
      GSL_REAL (x) = (ATOMIC)i;
      GSL_IMAG (x) = (ATOMIC)(i + 1234);
      y = FUNCTION (gsl_vector, get) (v, i);
      if (!GSL_COMPLEX_EQ (x, y))
	status = 1;
    };
  gsl_test (status, NAME (gsl_vector) "_get reads from array correctly");

  /* Now set stride to 2 */

  v->stride = 2 ;

  status = 0;

  for (i = 0; i < N / 2; i++)
    {
      BASE x, y;
      GSL_REAL (x) = (ATOMIC)(2 * i);
      GSL_IMAG (x) = (ATOMIC)(2 * i + 1234);
      y = FUNCTION (gsl_vector, get) (v, i);
      if (!GSL_COMPLEX_EQ (x, y))
	status = 1;
    };
  gsl_test (status, NAME (gsl_vector) "_get reads from array correctly with stride");

  for (i = 0; i < N / 2; i++)
    {
      BASE x;
      GSL_REAL (x) = (ATOMIC)i ;
      GSL_IMAG (x) = (ATOMIC)(i + 1234);
      FUNCTION (gsl_vector, set) (v, i, x);
    };

  status = 0;

  for (i = 0; i < N / 2; i++)
    {
      if (v->data[2 * 2 * i] != (ATOMIC) i || v->data[2 * 2 * i + 1] != (ATOMIC) (i + 1234))
	status = 1;
    };
  
  gsl_test (status, NAME (gsl_vector) "_set writes into array correctly with stride");

  /* Reset stride to 1 */

  v->stride = 1 ;

  for (i = 0; i < N; i++)
    {
      BASE x = ZERO;
      GSL_REAL (x) = (ATOMIC)i;
      GSL_IMAG (x) = (ATOMIC)(i + 1234);
      FUNCTION (gsl_vector, set) (v, i, x);
    };


  {
    BASE x, y, r, s ;
    GSL_REAL(x) = 2 ;
    GSL_IMAG(x) = 2 + 1234;
    GSL_REAL(y) = 5 ;
    GSL_IMAG(y) = 5 + 1234;

    FUNCTION (gsl_vector,swap_elements) (v, 2, 5) ;
    
    r = FUNCTION(gsl_vector,get)(v,2);
    s = FUNCTION(gsl_vector,get)(v,5);

    status = ! GSL_COMPLEX_EQ(r,y) ;
    status |= ! GSL_COMPLEX_EQ(s,x) ;
    
    FUNCTION (gsl_vector,swap_elements) (v, 2, 5) ;
    
    r = FUNCTION(gsl_vector,get)(v,2);
    s = FUNCTION(gsl_vector,get)(v,5);

    status |= ! GSL_COMPLEX_EQ(r,x) ;
    status |= ! GSL_COMPLEX_EQ(s,y) ;
  }

  gsl_test (status, NAME(gsl_vector) "_swap_elements" DESC " exchanges elements correctly") ;

  status = 0;

  FUNCTION (gsl_vector,reverse) (v) ;
  
  for (i = 0; i < N; i++)
    {
      BASE x,r ;
      GSL_REAL(x) = (ATOMIC)(N - i - 1) ;
      GSL_IMAG(x) = (ATOMIC)(N - i - 1 + 1234);

      r = FUNCTION (gsl_vector, get) (v, i);

      status |= !GSL_COMPLEX_EQ(r,x);
    }

  gsl_test (status, NAME(gsl_vector) "_reverse" DESC " reverses elements correctly") ;
  
  FUNCTION (gsl_vector, free) (v);	/* free whatever is in v */
}

void
FUNCTION (test, binary) (void)
{
  TYPE (gsl_vector) * v = FUNCTION (gsl_vector, calloc) (N);
  TYPE (gsl_vector) * w = FUNCTION (gsl_vector, calloc) (N);

  size_t i;

  {
    FILE *f = fopen ("test.dat", "wb");

    for (i = 0; i < N; i++)
      {
	BASE x = ZERO;
	GSL_REAL (x) = (ATOMIC)(N - i);
	GSL_IMAG (x) = (ATOMIC)(N - i + 1);
	FUNCTION (gsl_vector, set) (v, i, x);
      };

    FUNCTION (gsl_vector, fwrite) (f, v);

    fclose (f);
  }

  {
    FILE *f = fopen ("test.dat", "rb");

    FUNCTION (gsl_vector, fread) (f, w);

    status = 0;
    for (i = 0; i < N; i++)
      {
	if (w->data[2 * i] != (ATOMIC) (N - i) || w->data[2 * i + 1] != (ATOMIC) (N - i + 1))
	  status = 1;
      };
    fclose (f);
  }

  FUNCTION (gsl_vector, free) (v);
  FUNCTION (gsl_vector, free) (w);

  gsl_test (status, NAME (gsl_vector) "_write and read work correctly");

}


void
FUNCTION (test, trap) (void)
{
  TYPE (gsl_vector) * vc = FUNCTION (gsl_vector, alloc) (N);

  BASE z = {{(ATOMIC)1.2, (ATOMIC)3.4}};
  BASE z1 = {{(ATOMIC)4.5, (ATOMIC)6.7}};

  size_t j = 0;

  status = 0;
  FUNCTION (gsl_vector, set) (vc, j - 1, z);
  gsl_test (!status,
	    NAME (gsl_vector) "_set traps index below lower bound");

  status = 0;
  FUNCTION (gsl_vector, set) (vc, N + 1, z);
  gsl_test (!status,
	    NAME (gsl_vector) "_set traps index above upper bound");

  status = 0;
  FUNCTION (gsl_vector, set) (vc, N, z);
  gsl_test (!status, NAME (gsl_vector) "_set traps index at upper bound");

  status = 0;
  z1 = FUNCTION (gsl_vector, get) (vc, j - 1);
  gsl_test (!status,
	    NAME (gsl_vector) "_get traps index below lower bound");

  gsl_test (GSL_REAL (z1) != 0,
	    NAME (gsl_vector) "_get returns zero real below lower bound");
  gsl_test (GSL_IMAG (z1) != 0,
	    NAME (gsl_vector) "_get returns zero imag below lower bound");

  status = 0;
  z1 = FUNCTION (gsl_vector, get) (vc, N + 1);
  gsl_test (!status,
	    NAME (gsl_vector) "_get traps index above upper bound");
  gsl_test (GSL_REAL (z1) != 0,
	    NAME (gsl_vector) "_get returns zero real above upper bound");
  gsl_test (GSL_IMAG (z1) != 0,
	    NAME (gsl_vector) "_get returns zero imag above upper bound");

  status = 0;
  z1 = FUNCTION (gsl_vector, get) (vc, N);
  gsl_test (!status, NAME (gsl_vector) "_get traps index at upper bound");
  gsl_test (GSL_REAL (z1) != 0,
	    NAME (gsl_vector) "_get returns zero real at upper bound");
  gsl_test (GSL_IMAG (z1) != 0,
	    NAME (gsl_vector) "_get returns zero imag at upper bound");

  FUNCTION (gsl_vector, free) (vc);
}




