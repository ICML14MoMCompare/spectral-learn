#include "include/tests.h"

void run_tests() {
//	test_simplex();
//	test_pnfa();
//	test_nn();
//	test_tau();
//	test_maxK();
//	test_wa();
//	test_fista();
//	test_accel();
//	test_scaling();
	//test_fw();
	//test_fw_reg();
	//test_redsvd();
	//test_fw_tau();
	//test_glrp();
	//test_trunc_admm();
	//test_accel_wa();
	test_refined_admm();
}

void test_simplex() {
	cout << "========= TEST_SIMPLEX =========" << endl;
	VectorXd v(20);
	v << -9,-3,-10,-4,-2,-3,7,5,4,9,6,6,10,3,8,8,4,6,1,1;
	cout << "Original:" << endl << v.transpose() << endl;
	VectorXd w(20);
	double norm;
	// Larger numbers must get higher probability
	// Larger norms should yield more non-zeros
	w = v;
	norm = 1.0;
	simplex_projection(w,norm);
	cout << "Projection to norm (" << norm << "):" << endl << w.transpose() << endl;
	cout << "Sum: " << w.sum() << endl;
	w = v;
	norm = 5.0;
	simplex_projection(w,norm);
	cout << "Projection to norm (" << norm << "):" << endl << w.transpose() << endl;
	cout << "Sum: " << w.sum() << endl;
	w = v;
	norm = 10.0;
	simplex_projection(w,norm);
	cout << "Projection to norm (" << norm << "):" << endl << w.transpose() << endl;
	cout << "Sum: " << w.sum() << endl;
}

void test_pnfa() {
	cout << "========= TEST_PNFA =========" << endl;
	// A is already stochastic, projection should have no effect
	MatrixXd A(3,3);
	A << 0.4375,0.3750,0.1875,0.4500,0.1000,0.4500,0.6667,0.1333,0.2000;
	VectorXd a(3);
        a = VectorXd::Zero(3);
	MatrixXd B(3,3);
	B = A;
	pnfa_projection(A,a);
	cout << "Original:" << endl << B << endl;
	cout << "Projected:" << endl << A << endl;
	cout << "Distance: " << (A-B).squaredNorm() << endl;
	// Projections should keep proportionality
	A << 1,2,3,0,5,5,3,1,3;
	cout << "Original:" << endl << A << endl;
	pnfa_projection(A,a);
	cout << "Projected:" << endl << A << endl;
	// When there are two operators stacked together, projection must group
	// rows corresponding to same state
	MatrixXd C(6,3);
	C << 1,2,3,0,5,5,3,1,3,1,2,3,0,5,5,3,1,3;
	cout << "Original:" << endl << C << endl;
	pnfa_projection(C,a);
	cout << "Projected:" << endl << C << endl;
}

void test_nn() {
	cout << "========= TEST_NN =========" << endl;
	VectorXd u1(3), u2(3), v1(3), v2(3);
	u1 << 1,0,0;
	u2 << 0,1,0;
	v1 << 1,0,1;
	v2 << 1,0,-1;
	MatrixXd A(3,3);
	A = u1*v1.transpose() + 0.5*u2*v2.transpose();
	cout << "Original:" << endl << A << endl;
	// Since sigma_min = 0.5, this should shrink entries but not alter
	// zero/non-zero pattern
	nn_projection(A,0.1);
	cout << "Projected 0.1:" << endl << A << endl;
	// This should zero entries corresponding to u2*v2.transpose
	nn_projection(A,0.7);
	cout << "Projected 0.7:" << endl << A << endl;
}

void test_tau() {
	cout << "========= TEST_TAU =========" << endl;
	// Setup some known data (of full rank)
	MatrixXd H(3,3);
	H << 2,5,10,4,6,3,8,3,6;
	MatrixXd TrueA(3,3);
	TrueA << 0.4375,0.3750,0.1875,0.4500,0.1000,0.4500,0.6667,0.1333,0.2000;
	MatrixXd H_Sigma(3,3);
	H_Sigma = TrueA * H;
	MatrixXd A(3,3);
	VectorXd ainf(3);
        ainf = VectorXd::Zero(3);
	// For testing how close to a stochastic matrix we get
	VectorXd ones(3);
	ones = VectorXd::Ones(3);
	// As tau -> 0, must have A -> TrueA
	opt_icml12_admm(H,H_Sigma,ainf,A,100.0);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	cout << "Found matrix:" << endl << A << endl;
	cout << "Row sums: " << (A*ones).transpose() << endl;
	opt_icml12_admm(H,H_Sigma,ainf,A,1.0);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	cout << "Found matrix:" << endl << A << endl;
	cout << "Row sums: " << (A*ones).transpose() << endl;
	opt_icml12_admm(H,H_Sigma,ainf,A,0.1);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	cout << "Found matrix:" << endl << A << endl;
	cout << "Row sums: " << (A*ones).transpose() << endl;
	opt_icml12_admm(H,H_Sigma,ainf,A,0.01);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	cout << "Found matrix:" << endl << A << endl;
	cout << "Row sums: " << (A*ones).transpose() << endl;
	opt_icml12_admm(H,H_Sigma,ainf,A,0.001);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	cout << "Found matrix:" << endl << A << endl;
	cout << "Row sums: " << (A*ones).transpose() << endl;
	opt_icml12_admm(H,H_Sigma,ainf,A,0.0);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	cout << "Found matrix:" << endl << A << endl;
	cout << "Row sums: " << (A*ones).transpose() << endl;
}

void test_maxK() {
	cout << "========= TEST_MAXK =========" << endl;
	// Setup some known data (of full rank)
	MatrixXd H(3,3);
	H << 2,5,10,4,6,3,8,3,6;
	MatrixXd TrueA(3,3);
	TrueA << 0.4375,0.3750,0.1875,0.4500,0.1000,0.4500,0.6667,0.1333,0.2000;
	MatrixXd H_Sigma(3,3);
	H_Sigma = TrueA * H;
	MatrixXd A(3,3);
	VectorXd ainf(3);
        ainf = VectorXd::Zero(3);
	// For testing how close to a stochastic matrix we get
	VectorXd ones(3);
	ones = VectorXd::Ones(3);
	// As MaxK -> INF, must have A -> TrueA
	opt_icml12_admm(H,H_Sigma,ainf,A,0.0,10);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	cout << "Found matrix:" << endl << A << endl;
	cout << "Row sums: " << (A*ones).transpose() << endl;
	opt_icml12_admm(H,H_Sigma,ainf,A,0.0,50);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	cout << "Found matrix:" << endl << A << endl;
	cout << "Row sums: " << (A*ones).transpose() << endl;
	opt_icml12_admm(H,H_Sigma,ainf,A,0.0,100);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	cout << "Found matrix:" << endl << A << endl;
	cout << "Row sums: " << (A*ones).transpose() << endl;
	opt_icml12_admm(H,H_Sigma,ainf,A,0.0,1000);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	cout << "Found matrix:" << endl << A << endl;
	cout << "Row sums: " << (A*ones).transpose() << endl;
	opt_icml12_admm(H,H_Sigma,ainf,A,0.0,10000);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	cout << "Found matrix:" << endl << A << endl;
	cout << "Row sums: " << (A*ones).transpose() << endl;
}

void test_wa() {
	cout << "========= TEST_WA =========" << endl;
	// Setup some known data (of full rank)
	MatrixXd H(3,3);
	H << 2,5,10,4,6,3,8,3,6;
	MatrixXd TrueA(3,3);
	TrueA << 0.4375,0.3750,0.1875,0.4500,0.1000,0.4500,0.6667,0.1333,0.2000;
	MatrixXd H_Sigma(3,3);
	H_Sigma = TrueA * H;
	MatrixXd A(3,3);
	VectorXd ainf(3);
        ainf = VectorXd::Zero(3);
	// For testing how close to a stochastic matrix we get
	VectorXd ones(3);
	ones = VectorXd::Ones(3);
	// As MaxK -> INF, must have A -> TrueA
	opt_icml12_admm_wa(H,H_Sigma,A,0.0,10);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	cout << "Found matrix:" << endl << A << endl;
	opt_icml12_admm_wa(H,H_Sigma,A,0.0,50);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	cout << "Found matrix:" << endl << A << endl;
	opt_icml12_admm_wa(H,H_Sigma,A,0.0,100);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	cout << "Found matrix:" << endl << A << endl;
	opt_icml12_admm_wa(H,H_Sigma,A,0.1,100);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	cout << "Found matrix:" << endl << A << endl;
	opt_icml12_admm_wa(H,H_Sigma,A,0.5,100);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	cout << "Found matrix:" << endl << A << endl;
}

void test_fista() {
	cout << "========= TEST_FISTA =========" << endl;
	// Setup some known data (of full rank)
	MatrixXd H(3,3);
	H << 2,5,10,4,6,3,8,3,6;
	MatrixXd TrueA(3,3);
	TrueA << 0.4375,0.3750,0.1875,0.4500,0.1000,0.4500,0.6667,0.1333,0.2000;
	MatrixXd H_Sigma(3,3);
	H_Sigma = TrueA * H;
	MatrixXd A(3,3);
	VectorXd ainf(3);
        ainf = VectorXd::Zero(3);
	// As MaxK -> INF, must have A -> TrueA
	opt_icml12_fista_wa(H,H_Sigma,A,0.0,10);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	cout << "Found matrix:" << endl << A << endl;
	opt_icml12_fista_wa(H,H_Sigma,A,0.0,50);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	cout << "Found matrix:" << endl << A << endl;
	opt_icml12_fista_wa(H,H_Sigma,A,0.0,100);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	cout << "Found matrix:" << endl << A << endl;
	opt_icml12_fista_wa(H,H_Sigma,A,0.0,1000);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	cout << "Found matrix:" << endl << A << endl;
	opt_icml12_fista_wa(H,H_Sigma,A,0.1,100);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	cout << "Found matrix:" << endl << A << endl;
	opt_icml12_fista_wa(H,H_Sigma,A,0.1,1000);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	cout << "Found matrix:" << endl << A << endl;
	opt_icml12_fista_wa(H,H_Sigma,A,0.1,10000);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	cout << "Found matrix:" << endl << A << endl;
}

void test_accel() {
	cout << "========= TEST_ACCEL =========" << endl;
	// Setup some known data (of full rank)
	MatrixXd H(3,3);
	H << 2,5,10,4,6,3,8,3,6;
	MatrixXd TrueA(3,3);
	TrueA << 0.4375,0.3750,0.1875,0.4500,0.1000,0.4500,0.6667,0.1333,0.2000;
	MatrixXd H_Sigma(3,3);
	H_Sigma = TrueA * H;
	MatrixXd A(3,3);
	VectorXd ainf(3);
        ainf = VectorXd::Zero(3);
	// For testing how close to a stochastic matrix we get
	VectorXd ones(3);
	ones = VectorXd::Ones(3);
	// As MaxK -> INF, must have A -> TrueA
	accel_admm(H,H_Sigma,ainf,A,0.0,50);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	cout << "Found matrix:" << endl << A << endl;
	cout << "Row sums: " << (A*ones).transpose() << endl;
	accel_admm(H,H_Sigma,ainf,A,0.1);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	cout << "Found matrix:" << endl << A << endl;
	cout << "Row sums: " << (A*ones).transpose() << endl;
	accel_admm(H,H_Sigma,ainf,A,0.01);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	cout << "Found matrix:" << endl << A << endl;
	cout << "Row sums: " << (A*ones).transpose() << endl;
	accel_admm(H,H_Sigma,ainf,A,0.001);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	cout << "Found matrix:" << endl << A << endl;
	cout << "Row sums: " << (A*ones).transpose() << endl;
}

void test_scaling() {
	cout << "========= TEST_SCALING =========" << endl;
	MatrixXd H, A, TrueA, H_Sigma;
	VectorXd ainf;
	const int sizes[5] = { 10, 50, 100, 200, 500 };
	const int maxK1 = 500;
	const int maxK2 = 5000;
	const double tau = 0.001;
	for (int i = 1; i < 4; i++) {
		int d = sizes[i];
		ainf = VectorXd::Zero(d);
		H = MatrixXd::Random(d,d);
		TrueA = MatrixXd::Random(d,d).cwiseAbs();
		double alpha = 2*nuclear_norm(TrueA);
		pnfa_projection(TrueA,ainf);
		H_Sigma = TrueA * H;
		A = MatrixXd::Zero(d,d);
		accel_admm(H,H_Sigma,ainf,A,tau,maxK1);
		cout << "Error: " << (A - TrueA).squaredNorm() << endl;
		//opt_icml12_admm(H,H_Sigma,ainf,A,tau,maxK);
		//cout << "Error: " << (A - TrueA).squaredNorm() << endl;
		//opt_icml12_admm_wa(H,H_Sigma,A,tau,maxK);
		//cout << "Error: " << (A - TrueA).squaredNorm() << endl;
		//opt_icml12_fista_wa(H,H_Sigma,A,tau,maxK);
		//cout << "Error: " << (A - TrueA).squaredNorm() << endl;
		//opt_icml12_fw_wa(H,H_Sigma,A,alpha,maxK1);
		//cout << "Error: " << (A - TrueA).squaredNorm() << endl;
		//opt_icml12_fw_wa(H,H_Sigma,A,alpha,maxK2);
		//cout << "Error: " << (A - TrueA).squaredNorm() << endl;
		//fw_reg_wa(H,H_Sigma,A,tau,maxK1);
		//cout << "Error: " << (A - TrueA).squaredNorm() << endl;
		//fw_reg_wa(H,H_Sigma,A,tau,maxK2);
		//cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	}
	/*
	MatrixXf Hf, Af, TrueAf, H_Sigmaf;
	VectorXf ainff;
	for (int i = 0; i < 5; i++) {
		int d = sizes[i];
		ainff = VectorXf::Zero(d);
		Hf = MatrixXf::Random(d,d);
		TrueAf = MatrixXf::Random(d,d).cwiseAbs();
		double alpha = 1.5*nuclear_norm(TrueAf);
		pnfa_projection((MatrixXd)TrueAf,(VectorXd)ainff);
		H_Sigmaf = TrueAf * Hf;
		Af = MatrixXf::Zero(d,d);
		opt_icml12_fw_wa(Hf,H_Sigmaf,Af,alpha,maxK);
		cout << "Error: " << (Af - TrueAf).squaredNorm() << endl;
	}
	*/
	/*
	int d = 10;
	ainf = VectorXd::Zero(d);
	ones = VectorXd::Ones(d);
	H = MatrixXd::Random(d,d);
	TrueA = MatrixXd::Random(d,d).cwiseAbs();
	pnfa_projection(TrueA,ainf);
	H_Sigma = TrueA * H;
	A = MatrixXd::Zero(d,d);
	accel_admm(H,H_Sigma,ainf,A,tau,maxK);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	*/
}

void test_fw() {
	cout << "========= TEST_FW =========" << endl;
	// Setup some known data (of full rank)
	MatrixXf H(3,3);
	H << 2,5,10,4,6,3,8,3,6;
	MatrixXf TrueA(3,3);
	TrueA << 0.4375,0.3750,0.1875,0.4500,0.1000,0.4500,0.6667,0.1333,0.2000;
	double nnorm = nuclear_norm(TrueA);
	MatrixXf H_Sigma(3,3);
	H_Sigma = TrueA * H;
	MatrixXf A(3,3);
	VectorXf ainf(3);
        ainf = VectorXf::Zero(3);
	// As MaxK -> INF, must have A -> TrueA
	opt_icml12_fw_wa(H,H_Sigma,A,nnorm,10);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	//cout << "Found matrix:" << endl << A << endl;
	opt_icml12_fw_wa(H,H_Sigma,A,nnorm,50);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	//cout << "Found matrix:" << endl << A << endl;
	opt_icml12_fw_wa(H,H_Sigma,A,nnorm,100);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	//cout << "Found matrix:" << endl << A << endl;
	opt_icml12_fw_wa(H,H_Sigma,A,nnorm,1000);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	//cout << "Found matrix:" << endl << A << endl;
	opt_icml12_fw_wa(H,H_Sigma,A,nnorm,50000);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	//cout << "Found matrix:" << endl << A << endl;
	opt_icml12_fw_wa(H,H_Sigma,A,1.5*nnorm,50000);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	//cout << "Found matrix:" << endl << A << endl;
	opt_icml12_fw_wa(H,H_Sigma,A,1.01*nnorm,50000);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	//cout << "Found matrix:" << endl << A << endl;
	opt_icml12_fw_wa(H,H_Sigma,A,0.98*nnorm,50000);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	//cout << "Found matrix:" << endl << A << endl;
	opt_icml12_fw_wa(H,H_Sigma,A,0.6*nnorm,50000);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	//cout << "Found matrix:" << endl << A << endl;
}

void test_redsvd() {
	cout << "========= TEST_REDSVD =========" << endl;
	int d = 500;
	MatrixXf H = MatrixXf::Random(d,d);
	RedSVD S(H,5);
	cout << "First 5 singular values:" << endl;
	cout << S.singularValues() << endl;
}

void test_fw_reg() {
	cout << "========= TEST_FW_REG =========" << endl;
	// Setup some known data (of full rank)
	MatrixXf H(3,3);
	H << 2,5,10,4,6,3,8,3,6;
	MatrixXf TrueA(3,3);
	TrueA << 0.4375,0.3750,0.1875,0.4500,0.1000,0.4500,0.6667,0.1333,0.2000;
	MatrixXf H_Sigma(3,3);
	H_Sigma = TrueA * H;
	MatrixXf A(3,3);
	VectorXf ainf(3);
	double tau = 0.01;
        ainf = VectorXf::Zero(3);
	// As MaxK -> INF, must have A -> TrueA
	fw_reg_wa(H,H_Sigma,A,tau,10);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	//cout << "Found matrix:" << endl << A << endl;
	fw_reg_wa(H,H_Sigma,A,tau,50);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	//cout << "Found matrix:" << endl << A << endl;
	fw_reg_wa(H,H_Sigma,A,tau,100);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	//cout << "Found matrix:" << endl << A << endl;
	fw_reg_wa(H,H_Sigma,A,tau,1000);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	//cout << "Found matrix:" << endl << A << endl;
	fw_reg_wa(H,H_Sigma,A,tau,50000);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	//cout << "Found matrix:" << endl << A << endl;
	fw_reg_wa(H,H_Sigma,A,1.5*tau,50000);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	//cout << "Found matrix:" << endl << A << endl;
	fw_reg_wa(H,H_Sigma,A,1.01*tau,50000);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	//cout << "Found matrix:" << endl << A << endl;
	fw_reg_wa(H,H_Sigma,A,0.98*tau,50000);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	//cout << "Found matrix:" << endl << A << endl;
	fw_reg_wa(H,H_Sigma,A,0.6*tau,50000);
	cout << "Error: " << (A - TrueA).squaredNorm() << endl;
	//cout << "Found matrix:" << endl << A << endl;
}

void test_fw_tau() {
	cout << "========= TEST_FW_TAU =========" << endl;
	MatrixXd H, A, TrueA, H_Sigma;
	VectorXd ainf;
	const int d = 100;
	const int maxK1 = 500;
	const int maxK2 = 5000;
	//double taus[12] = {20, 10, 5, 2, 1, 0.5, 0.1, 0.05, 0.01, 0.005, 0.001, 0.0005};
	//double taus[12] = {1000, 800, 500, 200, 100, 80, 60, 20, 10, 5, 2, 1};
	//double taus[10] = {100, 90, 80, 70, 60, 50, 40, 30, 20, 10};
	double taus[10] = {30, 28, 26, 24, 22, 20, 18, 16, 14, 12};
	const int Ntaus = 10;
	for (int i = 0; i < Ntaus; i++) {
		double tau = taus[i];
		ainf = VectorXd::Zero(d);
		H = MatrixXd::Random(d,d);
		TrueA = MatrixXd::Random(d,d).cwiseAbs();
		pnfa_projection(TrueA,ainf);
		H_Sigma = TrueA * H;
		A = MatrixXd::Zero(d,d);
		fw_reg_wa(H,H_Sigma,A,tau,maxK1);
		cout << "Parameter Error: " << (A - TrueA).squaredNorm() << endl;
		cout << "Prediction Error: " << (A*H-H_Sigma).squaredNorm() << endl;
		fw_reg_wa(H,H_Sigma,A,tau,maxK2);
		cout << "Parameter Error: " << (A - TrueA).squaredNorm() << endl;
		cout << "Prediction Error: " << (A*H-H_Sigma).squaredNorm() << endl;
	}
}

void test_glrp() {
	cout << "========= TEST_GLRP =========" << endl;
	int d = 20;
	int r = 5;
	MatrixXd A;
	A = MatrixXd::Zero(d,d);
	VectorXd ainf;
	ainf = VectorXd::Zero(d);
	generate_low_rank_pnfa(A,ainf,r);
	//cout << "Generated PNFA:" << endl << A << endl;
	//cout << "Row sums: " << A.rowwise().sum().transpose() << endl;
}

void test_trunc_admm() {
	cout << "========= TEST_TRUNC_ADMM =========" << endl;
	MatrixXd H, A, TrueA, H_Sigma;
	VectorXd ainf;
	/*
	const int sizes[5] = { 10, 50, 100, 200, 500 };
	const int maxK = 500;
	const double tau = 0.001;
	for (int i = 0; i < 5; i++) {
		int d = sizes[i];
		int rank = (int) ceil(sqrt((float) d));
		ainf = VectorXd::Zero(d);
		H = MatrixXd::Random(d,d);
		TrueA = MatrixXd::Zero(d,d);
		generate_low_rank_pnfa(TrueA,ainf,rank);
		int maxRank = 2*rank;
		H_Sigma = TrueA * H;
		A = MatrixXd::Zero(d,d);
		accel_trunc_admm(H,H_Sigma,ainf,A,tau,maxK,maxRank);
		cout << "Parameter Error: " << (A - TrueA).squaredNorm() << endl;
		cout << "Prediction Error: " << (A*H - H_Sigma).squaredNorm() << endl;
	}
	*/
	int d = 200;
	int rank = 20;
	ainf = VectorXd::Zero(d);
	H = MatrixXd::Random(d,d);
	TrueA = MatrixXd::Zero(d,d);
	generate_low_rank_pnfa(TrueA,ainf,rank);
	//int maxRank = 2*rank;
	int maxRank = 5+rank;
	H_Sigma = TrueA * H;
	const int Ks[4] = {50, 100, 200, 400};
	const double Ts[3] = {0.01, 0.001, 0.0001};
	for (int i = 0; i < 3; i++) {
		double tau = Ts[i];
		for (int j = 0; j < 4; j++) {
			int maxK = Ks[j];
			A = MatrixXd::Zero(d,d);
			accel_trunc_admm(H,H_Sigma,ainf,A,tau,maxK,maxRank);
			cout << "Parameter Error: " << (A - TrueA).squaredNorm() << endl;
			cout << "Prediction Error: " << (A*H - H_Sigma).squaredNorm() << endl;
		}
	}
}

void test_accel_wa() {
	cout << "========= TEST_ACCEL_WA =========" << endl;
	MatrixXd H, A, TrueA, H_Sigma;
	VectorXd ainf;
	const int sizes[5] = { 10, 50, 100, 200, 500 };
	const int maxK1 = 100;
	const int maxK2 = 500;
	const double tau = 0.001;
	for (int i = 0; i < 5; i++) {
		int d = sizes[i];
		int rank = 3 + (int) ceil((double) d / 10.0);
		int maxRank = 2*rank;
		ainf = VectorXd::Zero(d);
		H = MatrixXd::Random(d,d);
		TrueA = MatrixXd::Zero(d,d);
		generate_low_rank_pnfa(TrueA,ainf,rank);
		H_Sigma = TrueA * H;
		A = MatrixXd::Zero(d,d);
		accel_trunc_admm_wa(H,H_Sigma,A,tau,maxK1,maxRank);
		cout << "Parameter Error: " << (A - TrueA).squaredNorm() << endl;
		cout << "Prediction Error: " << (A*H - H_Sigma).squaredNorm() << endl;
		accel_trunc_admm_wa(H,H_Sigma,A,tau,maxK2,maxRank);
		cout << "Parameter Error: " << (A - TrueA).squaredNorm() << endl;
		cout << "Prediction Error: " << (A*H - H_Sigma).squaredNorm() << endl;
	}
}

void test_refined_admm() {
	MatrixXd H, A, TrueA, H_Sigma;
	VectorXd ainf;
	int d = 400;
	int rank = 20;
	const int Ks[4] = {100, 200};
	const double Ts[3] = {0.01, 0.001, 0.0001};
	ainf = VectorXd::Zero(d);
	H = MatrixXd::Random(d,d);
	TrueA = MatrixXd::Zero(d,d);
	generate_low_rank_pnfa(TrueA,ainf,rank);
	H_Sigma = TrueA * H;
	cout << "========= TEST_REFINED_ADMM (FOR PNFA) =========" << endl;
	int maxRank = 2*rank;
	for (int i = 0; i < 3; i++) {
		double tau = Ts[i];
		for (int j = 0; j < 2; j++) {
			int maxK = Ks[j];
			A = MatrixXd::Zero(d,d);
			accel_trunc_admm(H,H_Sigma,ainf,A,tau,maxK,maxRank);
			cout << "Parameter Error: " << (A - TrueA).squaredNorm() << endl;
			cout << "Prediction Error: " << (A*H - H_Sigma).squaredNorm() << endl;
			// Run 2 iterations of ADMM with exact SVD to refine
			// results
			accel_admm(H,H_Sigma,ainf,A,tau,4,true);
			cout << "(Refined) Parameter Error: " << (A - TrueA).squaredNorm() << endl;
			cout << "(Refined) Prediction Error: " << (A*H - H_Sigma).squaredNorm() << endl;
		}
	}
	cout << "========= TEST_REFINED_ADMM (FOR WA) =========" << endl;
	for (int i = 0; i < 3; i++) {
		double tau = Ts[i];
		for (int j = 0; j < 2; j++) {
			int maxK = Ks[j];
			A = MatrixXd::Zero(d,d);
			accel_trunc_admm_wa(H,H_Sigma,A,tau,maxK,maxRank);
			cout << "Parameter Error: " << (A - TrueA).squaredNorm() << endl;
			cout << "Prediction Error: " << (A*H - H_Sigma).squaredNorm() << endl;
			// Run 2 iterations of ADMM with exact SVD to refine
			// results
			accel_admm_wa(H,H_Sigma,A,tau,4,true);
			cout << "(Refined) Parameter Error: " << (A - TrueA).squaredNorm() << endl;
			cout << "(Refined) Prediction Error: " << (A*H - H_Sigma).squaredNorm() << endl;
		}
	}
}

