#include "include/fista.h"

//#define DEBUG_IT
#define DEBUG_CONV
#define DEBUG_INPUT

#define EPSCONV 1E-6

#define FNAME "(opt_icml12_fista_wa)"
void opt_icml12_fista_wa(MatrixXd & H, MatrixXd & H_Sigma, MatrixXd & A_Sigma, double tau, int maxK) {

	// Check dimensions (A_Sigma * H == H_Sigma)
	int pH = H.rows();
	int sH = H.cols();
	int pHs = H_Sigma.rows();
	int sHs = H_Sigma.cols();
	int rA = A_Sigma.rows();
	int cA = A_Sigma.cols();
	if (cA != pH or rA != pHs or sH != sHs) {
		cerr << FNAME << " ERROR: Dimensions mismatch!" << endl;
		return;
	}

#ifdef DEBUG_INPUT
	cerr << FNAME << " Called with input:" << endl;
	cerr << "::::: H\t(" << pH << "x" << sH << ")" << endl;
	cerr << "::::: H_Sigma\t(" << pHs << "x" << sHs << ")" << endl;
	cerr << "::::: A_Sigma\t(" << rA << "x" << cA << ")" << endl;
	cerr << "::::: tau = " << tau << endl;
	cerr << "::::: maxK = " << maxK << endl;
#endif
	
	// Placeholders for variables
	MatrixXd X = MatrixXd::Zero(rA, cA);
	MatrixXd Xprev = MatrixXd::Zero(rA, cA);
	MatrixXd Y = MatrixXd::Zero(rA, cA);
	double gamma = 1.0;
	double gammaPrev = 0.0;
	double diff = 0.0;
	
	// Precomputations
	double L = 2.0*(H*H.transpose()).norm();
	MatrixXd Aux1 = 2.0*H_Sigma*H.transpose() / L;
	MatrixXd Aux2 = MatrixXd::Identity(pH,pH) - (2.0*H*H.transpose() / L);

	int k = 0;
	do {
		Xprev = X;
		gammaPrev = gamma;

		// Proximal operator
		X = Y*Aux2 + Aux1;
		nn_projection(X,tau/L);
		// Step update
		gamma = (1.0 + sqrt(1.0 + 4.0*gamma*gamma)) / 2.0;
		// Dual variable
		Y = X + (gammaPrev - 1.0)*(X - Xprev)/gamma;

		diff = (Xprev - X).norm();
		k++;
	} while (k < maxK && diff > EPSCONV);

#ifdef DEBUG_CONV
	cerr << FNAME << " Converged after " << k << " iterations ";
	if (k >= maxK) cerr << "by MAXIMUM ITERATIONS";
	else cerr << "by CONVERGENCE CRITERION";
	cerr << endl;
#endif
		
	A_Sigma = X;
}

