#ifndef UTIL_H
#define UTIL_H

#include "Eigen/Dense"
#include "projections.h"
#include "redsvd/redsvd.hpp"

using namespace Eigen;
using namespace REDSVD;

double nuclear_norm(MatrixXd & m);
double nuclear_norm(MatrixXf & m);
double fast_nuclear_norm(MatrixXd & m, int R);
void generate_low_rank_pnfa(MatrixXd & A, VectorXd & ainf, int rank);

#endif

