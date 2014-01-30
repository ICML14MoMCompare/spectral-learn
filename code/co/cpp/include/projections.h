#ifndef PROJECTIONS_H
#define PROJECTIONS_H

#include <iostream>
#include <vector>
#include <algorithm>
#include "Eigen/Dense"
#include "Eigen/SVD"
#include "redsvd/redsvd.hpp"

using namespace std;
using namespace Eigen;
using namespace REDSVD;

bool simplex_projection(VectorXd & v, double norm);
bool pnfa_projection(MatrixXd & X, VectorXd & a);
bool nn_projection(MatrixXd & X, double thr);
bool fast_nn_projection(MatrixXd & X, double thr, int R);

#endif
