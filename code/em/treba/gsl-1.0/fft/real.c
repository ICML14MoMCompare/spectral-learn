#include <config.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_complex.h>

#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_real_float.h>

#define BASE_DOUBLE
#include "templates_on.h"
#include "real_pass.h"
#include "real_init.c"
#include "real_main.c"
#include "real_pass_2.c"
#include "real_pass_3.c"
#include "real_pass_4.c"
#include "real_pass_5.c"
#include "real_pass_n.c"
#include "real_radix2.c"
#include "real_unpack.c"
#include "templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "templates_on.h"
#include "real_pass.h"
#include "real_init.c"
#include "real_main.c"
#include "real_pass_2.c"
#include "real_pass_3.c"
#include "real_pass_4.c"
#include "real_pass_5.c"
#include "real_pass_n.c"
#include "real_radix2.c"
#include "real_unpack.c"
#include "templates_off.h"
#undef  BASE_FLOAT
