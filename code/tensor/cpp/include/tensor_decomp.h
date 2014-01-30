#ifndef TENSOR_DECOMP_H
#define TENSOR_DECOMP_H

#include "ltensor_defs.h"
#include "Eigen/Dense"
#include <iostream>
#include "util.h"

using namespace std;
using namespace Eigen;

#define DEFAULT_POWER_IT 100

double tensor_decomp(Ltensor & T, Lmatrix & evecs, Lvector & evals, int K, Ltensor & R, int POWER_IT = DEFAULT_POWER_IT);
double power_method(Ltensor & T, Lvector & evec, int M);

#endif

