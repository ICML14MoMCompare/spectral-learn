#ifndef TENSOR_WRAPPER_H
#define TENSOR_WRAPPER_H

#include <iostream>
#include "ltensor_defs.h"
#include "tensor_decomp.h"
#include "Eigen/Eigenvalues"
#include "learn_Otilde.h"
#include "projections.h"
#include "tensor_symmetrize.h"

using namespace std;
using namespace Eigen;

extern "C"
{
void learn_tensor(double* Xdata, double* Ydata, int n, double * Odata, double* gammadata, int num_symbols, bool symmetrizeY);
bool test_transfer_tensor(double* tensorData);
bool test_Y_symmetry(double *Ydata, int num_symbols);
void get_real_HMM_hankels(double* Hpsdata, double* Hpsigdata, double* Hsigsdata, double* tensorData, double *HpsnontensorData, double* HppsData);
void do_simplex_projection(double* vectorData, double norm, int length);
void do_probmat_projection(double* matrixData, double* normData, int m, int n);
}

#endif
