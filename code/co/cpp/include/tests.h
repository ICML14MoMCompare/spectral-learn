#ifndef TESTS_H
#define TESTS_H

#include <iostream>
#include "Eigen/Dense"
#include "redsvd/redsvd.hpp"
#include "projections.h"
#include "admm.h"
#include "fista.h"
#include "frank_wolfe.h"
#include "util.h"

using namespace std;
using namespace Eigen;
using namespace REDSVD;

void run_tests();
void test_simplex();
void test_pnfa();
void test_nn();
void test_tau();
void test_maxK();
void test_wa();
void test_fista();
void test_accel();
void test_trunc_admm();
void test_scaling();
void test_fw();
void test_fw_reg();
void test_fw_tau();
void test_redsvd();
void test_glrp();
void test_accel_wa();
void test_refined_admm();

#endif

