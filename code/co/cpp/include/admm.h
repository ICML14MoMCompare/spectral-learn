#ifndef ADMM_H
#define ADMM_H

#include <iostream>
#include "Eigen/Dense"
#include "projections.h"
#include "util.h"

using namespace std;
using namespace Eigen;

void opt_icml12_admm(MatrixXd & H, MatrixXd & H_Sigma, VectorXd & ainf, MatrixXd & A_Sigma, double tau = 0.1, int maxK = 50);
void opt_icml12_admm_wa(MatrixXd & H, MatrixXd & H_Sigma, MatrixXd & A_Sigma, double tau = 0.1, int maxK = 50);
void accel_admm(MatrixXd & H, MatrixXd & H_Sigma, VectorXd & ainf, MatrixXd & A_Sigma, double tau = 0.1, int maxK = 50, bool Ainitial = false);
void accel_trunc_admm(MatrixXd & H, MatrixXd & H_Sigma, VectorXd & ainf, MatrixXd & A_Sigma, double tau, int maxK, int maxRank);
void accel_trunc_admm_wa(MatrixXd & H, MatrixXd & H_Sigma, MatrixXd & A_Sigma, double tau, int maxK, int maxRank);
void accel_admm_wa(MatrixXd & H, MatrixXd & H_Sigma, MatrixXd & A_Sigma, double tau, int maxK, bool Ainitial = false);

#endif

