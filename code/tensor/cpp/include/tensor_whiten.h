#ifndef TENSOR_WHITEN_H
#define TENSOR_WHITEN_H

#include "util.h"
#include "ltensor_defs.h"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include <iostream>

using namespace std;
using namespace Eigen;

double whiten_tensor(MatrixXd & X, Ltensor & Y, int n, Ltensor & Z, MatrixXd & pinv_trans_W);

#endif

