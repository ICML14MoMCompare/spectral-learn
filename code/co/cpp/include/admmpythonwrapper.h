#ifndef ADMM_WRAPPER_H
#define ADMM_WRAPPER_H

#include <iostream>
#include "Eigen/Dense"
#include "admm.h"

using namespace std;
using namespace Eigen;

extern "C"
{
void admm_pnfa_learn(double* H, double* H_Sigma, double* ainf, double* A_Sigma, int basisSize, int numSymbols, double tau, int maxK, int maxDim);
void admm_wfa_learn(double* H, double* H_Sigma, double* A_Sigma, int basisSize, int numSymbols, double tau = 0.1, int maxK = 50, int maxDim = 50);
}

#endif