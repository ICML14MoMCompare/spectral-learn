/* vector/test_source.c
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
  TYPE (gsl_vector) * v;
  size_t i;

  v = FUNCTION (gsl_vector, calloc) (N);

  gsl_test (v->data == 0, NAME (gsl_vector) "_alloc returns valid pointer");
  gsl_test (v->size != N, NAME (gsl_vector) "_alloc returns valid size");
  gsl_test (v->stride != 1, NAME (gsl_vector) "_alloc returns unit stride");

  for (i = 0; i < N; i++)
    {
      FUNCTION (gsl_vector, set) (v, i, (ATOMIC) i);
    }


  status = 0;

  for (i = 0; i < N; i++)
    {
      if (v->data[i] != (ATOMIC) i)
	status = 1;
    };
  
  gsl_test (status,
	    NAME (gsl_vector) "_set" DESC " writes into array correctly");


  status = 0;
  for (i = 0; i < N; i++)
    {
      if (FUNCTION (gsl_vector, get) (v, i) != (ATOMIC) i)
	status = 1;
    };
  gsl_test (status,
	    NAME (gsl_vector) "_get" DESC " reads from array correctly");


  /* Now set stride to 2 */

  v->stride = 2;
  
  status = 0;
  for (i = 0; i < N / 2; i++)
    {
      if (FUNCTION (gsl_vector, get) (v, i) != (ATOMIC) (2 * i))
	status = 1;
    };
  gsl_test (status, NAME (gsl_vector) "_get" DESC " reads correctly with stride");

  for (i = 0; i < N / 2; i++)
    {
      FUNCTION (gsl_vector, set) (v, i, (ATOMIC) (i + 1000));
    };
  
  status = 0;

  for (i = 0; i < N / 2; i++)
    {
      if (v->data[2 * i] != (ATOMIC) (i + 1000))
	status = 1;
    };
  
  gsl_test (status, NAME (gsl_vector) "_set" DESC " writes correctly with stride");

  /* Reset stride to 1 */

  v->stride = 1 ;

  for (i = 0; i < N; i++)
    {
      FUNCTION (gsl_vector, set) (v, i, (ATOMIC) i);
    }

  FUNCTION (gsl_vector,swap_elements) (v, 2, 5) ;
  
  status = (FUNCTION(gsl_vector,get)(v,2) != 5) ;
  status |= (FUNCTION(gsl_vector,get)(v,5) != 2) ;

  FUNCTION (gsl_vector,swap_elements) (v, 2, 5) ;

  status |= (FUNCTION(gsl_vector,get)(v,2) != 2) ;
  status |= (FUNCTION(gsl_vector,get)(v,5) != 5) ;

  gsl_test (status, NAME(gsl_vector) "_swap_elements" DESC " exchanges elements correctly") ;

  status = 0;

  FUNCTION (gsl_vector,reverse) (v) ;
  
  for (i = 0; i < N; i++)
    {
      status |= (FUNCTION (gsl_vector, get) (v, i) !=  (ATOMIC) (N - i - 1));
    }

  gsl_test (status, NAME(gsl_vector) "_reverse" DESC " reverses elements correctly") ;

  {
    BASE exp_max = FUNCTION(gsl_vector,get)(v, 0);
    BASE exp_min = FUNCTION(gsl_vector,get)(v, 0);
    size_t exp_imax = 0, exp_imin = 0;

    for (i = 0; i < N; i++)
      {
        BASE k = FUNCTION(gsl_vector, get) (v, i) ;
        if (k < exp_min) {
          exp_min = FUNCTION(gsl_vector, get) (v, i);
          exp_imin = i;
        }
      }

    for (i = 0; i < N; i++)
      {
        BASE k = FUNCTION(gsl_vector, get) (v, i) ;
        if (k > exp_max) {
          exp_max = FUNCTION(gsl_vector, get) (v, i) ;
          exp_imax = i;
        } 
      }

    {
      BASE max = FUNCTION(gsl_vector, max) (v) ;

      gsl_test (max != exp_max, NAME(gsl_vector) "_max returns correct maximum value");
    }

    {
      BASE min = FUNCTION(gsl_vector, min) (v) ;
      
      gsl_test (min != exp_min, NAME(gsl_vector) "_min returns correct minimum value");
    }

    {
      BASE min, max;
      FUNCTION(gsl_vector, minmax) (v, &min, &max);

      gsl_test (max != exp_max, NAME(gsl_vector) "_minmax returns correct maximum value");
      gsl_test (min != exp_min, NAME(gsl_vector) "_minmax returns correct minimum value");
    }


    {
      size_t imax =  FUNCTION(gsl_vector, max_index) (v) ;

      gsl_test (imax != exp_imax, NAME(gsl_vector) "_max_index returns correct maximum i");
    }

    {
      size_t imin = FUNCTION(gsl_vector, min_index) (v) ;

      gsl_test (imin != exp_imin, NAME(gsl_vector) "_min_index returns correct minimum i");
    }

    {
      size_t imin, imax;

      FUNCTION(gsl_vector, minmax_index) (v,  &imin, &imax);

      gsl_test (imax != exp_imax, NAME(gsl_vector) "_minmax_index returns correct maximum i");
      gsl_test (imin != exp_imin, NAME(gsl_vector) "_minmax_index returns correct minimum i");
    }
  }


  {
    TYPE (gsl_vector) * a = FUNCTION (gsl_vector, calloc) (N);
    TYPE (gsl_vector) * b = FUNCTION (gsl_vector, calloc) (N);
    
    for (i = 0; i < N; i++)
      {
        FUNCTION (gsl_vector, set) (a, i, (BASE)(3 + i));
        FUNCTION (gsl_vector, set) (b, i, (BASE)(3 + 2 * i));
      }
    
    FUNCTION(gsl_vector, memcpy) (v, a);
    FUNCTION(gsl_vector, add) (v, b);
    
    {
      int status = 0;
      
      for (i = 0; i < N; i++)
        {
          BASE r = FUNCTION(gsl_vector,get) (v,i);
          BASE x = FUNCTION(gsl_vector,get) (a,i);
          BASE y = FUNCTION(gsl_vector,get) (b,i);
          BASE z = x + y;
          if (r != z)
            status = 1;
        }
      gsl_test (status, NAME (gsl_vector) "_add adds correctly");
    }


    FUNCTION(gsl_vector, memcpy) (v, a);
    FUNCTION(gsl_vector, sub) (v, b);
    
    {
      int status = 0;
      
      for (i = 0; i < N; i++)
        {
          BASE r = FUNCTION(gsl_vector,get) (v,i);
          BASE x = FUNCTION(gsl_vector,get) (a,i);
          BASE y = FUNCTION(gsl_vector,get) (b,i);
          BASE z = x - y;
          if (r != z)
            status = 1;
        }
      gsl_test (status, NAME (gsl_vector) "_sub subtracts correctly");
    }

    FUNCTION(gsl_vector, memcpy) (v, a);
    FUNCTION(gsl_vector, mul) (v, b);
    
    {
      int status = 0;
      
      for (i = 0; i < N; i++)
        {
          BASE r = FUNCTION(gsl_vector,get) (v,i);
          BASE x = FUNCTION(gsl_vector,get) (a,i);
          BASE y = FUNCTION(gsl_vector,get) (b,i);
          BASE z = x * y;
          if (r != z)
            status = 1;
        }
      gsl_test (status, NAME (gsl_vector) "_mul multiplies correctly");
    }

    FUNCTION(gsl_vector, memcpy) (v, a);
    FUNCTION(gsl_vector, div) (v, b);
    
    {
      int status = 0;
      
      for (i = 0; i < N; i++)
        {
          BASE r = FUNCTION(gsl_vector,get) (v,i);
          BASE x = FUNCTION(gsl_vector,get) (a,i);
          BASE y = FUNCTION(gsl_vector,get) (b,i);
          BASE z = x / y;
          if (fabs(r - z) > 2 * GSL_FLT_EPSILON * fabs(z))
            status = 1;
        }
      gsl_test (status, NAME (gsl_vector) "_div divides correctly");
    }


    FUNCTION(gsl_vector, free) (a);
    FUNCTION(gsl_vector, free) (b);
  }




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
	FUNCTION (gsl_vector, set) (v, i, (ATOMIC) (N - i));
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
	if (w->data[i] != (ATOMIC) (N - i))
	  status = 1;
      };

    gsl_test (status, NAME (gsl_vector) "_write and read work correctly");

    fclose (f);
  }

  FUNCTION (gsl_vector, free) (v);	/* free whatever is in v */
  FUNCTION (gsl_vector, free) (w);	/* free whatever is in w */
}

void
FUNCTION (test, trap) (void)
{
  TYPE (gsl_vector) * v = FUNCTION (gsl_vector, alloc) (N);

  size_t j = 0;
  double x;

  status = 0;
  FUNCTION (gsl_vector, set) (v, j - 1, (ATOMIC)0);
  gsl_test (!status, NAME (gsl_vector) "_set traps index below lower bound");

  status = 0;
  FUNCTION (gsl_vector, set) (v, N + 1, (ATOMIC)0);
  gsl_test (!status, NAME (gsl_vector) "_set traps index above upper bound");

  status = 0;
  FUNCTION (gsl_vector, set) (v, N, (ATOMIC)0);
  gsl_test (!status, NAME (gsl_vector) "_set traps index at upper bound");

  status = 0;
  x = FUNCTION (gsl_vector, get) (v, j - 1);
  gsl_test (!status, NAME (gsl_vector) "_get traps index below lower bound");
  gsl_test (x != 0,
	 NAME (gsl_vector) "_get returns zero for index below lower bound");

  status = 0;
  x = FUNCTION (gsl_vector, get) (v, N + 1);
  gsl_test (!status, NAME (gsl_vector) "_get traps index above upper bound");
  gsl_test (x != 0,
	 NAME (gsl_vector) "_get returns zero for index above upper bound");

  status = 0;
  x = FUNCTION (gsl_vector, get) (v, N);
  gsl_test (!status, NAME (gsl_vector) "_get traps index at upper bound");
  gsl_test (x != 0,
	    NAME (gsl_vector) "_get returns zero for index at upper bound");

  FUNCTION (gsl_vector, free) (v);	/* free whatever is in v */
}





