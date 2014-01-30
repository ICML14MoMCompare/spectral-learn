#ifndef LEARN_OTILDE_H
#define LEARN_OTILDE_H

#include "Eigen/Dense"
#include "ltensor_defs.h"
#include "tensor_decomp.h"
#include "tensor_whiten.h"
#include <iostream>

using namespace std;
using namespace Eigen;

double learn_Otilde(MatrixXd & X, Ltensor & Y, int n, MatrixXd & O, VectorXd & gamma, bool symmetrizeY = true);

#endif

