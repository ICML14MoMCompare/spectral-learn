/* Author: G. Jungman
 */
#include <config.h>
#include <gsl/gsl_ieee_utils.h>

#include "gsl_qrng.h"


int test_sobol(void)
{
  int status = 0;
  double v[3];
  /* int i; */

  /* test in dimension 2 */
  gsl_qrng * g = gsl_qrng_alloc(gsl_qrng_sobol, 2);
  gsl_qrng_get(g, v);
  gsl_qrng_get(g, v);
  gsl_qrng_get(g, v);
  status += ( v[0] != 0.25 || v[1] != 0.75 );
  gsl_qrng_get(g, v);
  status += ( v[0] != 0.375 || v[1] != 0.375 );
  gsl_qrng_free(g);

  /* test in dimension 3 */
  g = gsl_qrng_alloc(gsl_qrng_sobol, 3);
  gsl_qrng_get(g, v);
  gsl_qrng_get(g, v);
  gsl_qrng_get(g, v);
  status += ( v[0] != 0.25 || v[1] != 0.75 || v[2] != 0.25 );
  gsl_qrng_get(g, v);
  status += ( v[0] != 0.375 || v[1] != 0.375 || v[2] != 0.625 );
  gsl_qrng_free(g);

  return status;
}


int test_nied2(void)
{
  int status = 0;
  double v[3];
  /* int i; */

  /* test in dimension 2 */
  gsl_qrng * g = gsl_qrng_alloc(gsl_qrng_niederreiter_2, 2);
  gsl_qrng_get(g, v);
  gsl_qrng_get(g, v);
  gsl_qrng_get(g, v);
  status += ( v[0] != 0.75 || v[1] != 0.25 );
  gsl_qrng_get(g, v);
  status += ( v[0] != 0.25 || v[1] != 0.75 );
  gsl_qrng_get(g, v);
  gsl_qrng_get(g, v);
  gsl_qrng_get(g, v);
  status += ( v[0] != 0.625 || v[1] != 0.125 );
  gsl_qrng_free(g);

  /* test in dimension 3 */
  g = gsl_qrng_alloc(gsl_qrng_niederreiter_2, 3);
  gsl_qrng_get(g, v);
  gsl_qrng_get(g, v);
  gsl_qrng_get(g, v);
  status += ( v[0] != 0.75 || v[1] != 0.25 || v[2] != 0.3125 );
  gsl_qrng_get(g, v);
  status += ( v[0] != 0.25 || v[1] != 0.75 || v[2] != 0.5625 );
  gsl_qrng_get(g, v);
  gsl_qrng_get(g, v);
  gsl_qrng_get(g, v);
  status += ( v[0] != 0.625 || v[1] != 0.125 || v[2] != 0.6875 );
  gsl_qrng_free(g);

  return status;
}


int main()
{
  int status = 0;

  gsl_ieee_env_setup ();

  status += test_sobol();
  status += test_nied2();

  return status;
}
