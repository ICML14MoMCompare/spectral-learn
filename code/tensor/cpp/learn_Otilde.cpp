#include "include/learn_Otilde.h"

#define DEBUG_STATS

// Input: X (d*d), Y (d*d*d), n = #number of states
// Output: O (d*n) (scaled singular vectors), gamma (n) (singular values)
// Returns: residual from tensor decomposition
double learn_Otilde(MatrixXd & X, Ltensor & Y, int n, MatrixXd & O, VectorXd & gamma, bool symmetrizeY) {
	// TODO Check dimensions. But all other functions do, so not really
	// needed
	int d = X.rows();
	// We need X to be symmetric, and Y is optional
	// (though it may be a good idea)
	MatrixXd sX;
	sX = symmetrize_matrix(X);
	Ltensor sY;
	if (symmetrizeY) {
		sY = symmetrize_tensor(Y);
	} else {
		sY = Y;
	}
	// Whitening
	Ltensor Z(n,n,n);
	MatrixXd Wtpinv;
	Wtpinv = MatrixXd::Zero(d,n);
	double w_residual = whiten_tensor(sX,sY,n,Z,Wtpinv);
	// Decompose
	Lvector lgamma(n);
	Lmatrix evecs(d,n);
	Ltensor deflatedZ(n,n,n);
	double d_residual = tensor_decomp(Z, evecs, lgamma, n, deflatedZ);
	// Go to eigen :(
	MatrixXd E;
	E = MatrixXd::Zero(n,n);
	lmatrix2eigen(evecs,E);
	gamma = VectorXd::Zero(n);
	lvector2eigen(lgamma,gamma);
	// Build matrix O
	//O = MatrixXd::Zero(d,n);
	//for (int i = 0; i < n; i++) {
	//	O.col(i) = gamma(i)*Wtpinv*E.col(i);
	//}
	O = Wtpinv*E*gamma.asDiagonal();
	return d_residual;
}

