#ifndef FISTA_H
#define FISTA_H

#include <iostream>
#include "Eigen/Dense"
#include "projections.h"

using namespace std;
using namespace Eigen;

void opt_icml12_fista_wa(MatrixXd & H, MatrixXd & H_Sigma, MatrixXd & A_Sigma, double tau, int maxK);

#endif
