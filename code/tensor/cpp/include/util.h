#ifndef MYUTIL_H
#define MYUTIL_H

#include "Eigen/Dense"
#include "ltensor_defs.h"

#include "boost/random/mersenne_twister.hpp"
#include "boost/random/normal_distribution.hpp"

#include <iostream>
#include <vector>

using namespace std;
using namespace Eigen;

void random_spherical_vector(Lvector & v, int d);
/*
void eigen2lvector(VectorXd & e, Lvector & l);
*/
void eigen2lmatrix(const MatrixXd & e, Lmatrix & l);
void lmatrix2eigen(Lmatrix & l, MatrixXd & e);
void lvector2eigen(Lvector & l, VectorXd & e);
void ltensor2eigenmatrix(Ltensor & t, MatrixXd & m);
double frob_norm_ltensor(Ltensor & t);
double sum_ltensor(Ltensor & t);

#define DEBUG_TEN_SYM true
bool is_tensor_symmetric(Ltensor & t, bool debug = DEBUG_TEN_SYM);
#define DEBUG_MAT_SYM true
bool is_matrix_symmetric(MatrixXd & m, bool debug = DEBUG_MAT_SYM);

Ltensor symmetrize_tensor(Ltensor & t);
MatrixXd symmetrize_matrix(MatrixXd & m);

double evalwa(const VectorXd & a0, const VectorXd & ainf, const vector<MatrixXd> & wa, const vector<int> s);

#endif
