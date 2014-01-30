#ifndef TENSOR_SYMMETRIZE_H
#define TENSOR_SYMMETRIZE_H

#include "Eigen/Dense"
#include "Eigen/SVD"

using namespace Eigen;

void computeNfromLowRankInverse(MatrixXd & H, int n, MatrixXd & N);

#endif

