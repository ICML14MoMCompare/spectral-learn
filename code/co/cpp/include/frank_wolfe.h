#ifndef FW_H
#define FW_H

#include <iostream>
#include "Eigen/Dense"
#include "redsvd/redsvd.hpp"

using namespace std;
using namespace Eigen;
using namespace REDSVD;

void opt_icml12_fw_wa(MatrixXd & H, MatrixXd & H_Sigma, MatrixXd & A_Sigma, double alpha, int maxK);
void opt_icml12_fw_wa(MatrixXf & H, MatrixXf & H_Sigma, MatrixXf & A_Sigma, double alpha, int maxK);
void fw_reg_wa(MatrixXd & H, MatrixXd & H_Sigma, MatrixXd & A_Sigma, double tau, int maxK);
void fw_reg_wa(MatrixXf & H, MatrixXf & H_Sigma, MatrixXf & A_Sigma, double tau, int maxK);

#endif

