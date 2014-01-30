#include "include/tensor_whiten.h"

// TODO Agree on input types
// Take as input X (DxD), Y (DxDxD), n <= rank(X)
// (note: X and Y should be symmetric)
// Return whitened tensor Z, whitening matrix W (*), and residual error of X
// (*) actually, the pseudo-inverse of its transpose, which will be needed later
#define EPS_SVAL 1E-6
#define FNAME "(whiten_tensor)"
double whiten_tensor(MatrixXd & X, Ltensor & Y, int n, Ltensor & Z, MatrixXd & pinv_trans_W) {
	// Get dimensions
	int d1x = X.rows();
	int d2x = X.cols();
	int d1y = Y.get_dim1();
	int d2y = Y.get_dim2();
	int d3y = Y.get_dim3();
	// Check dimensions
	if (d1x != d2x or d2x != d1y or d1y != d2y or d2y != d3y or n > d1x) {
		cerr << FNAME << " Dimensions mismatch!" << endl;
		throw;
	}
	int d = d1x;
	// Check symmetry
	if (!is_matrix_symmetric(X) or !is_tensor_symmetric(Y)) {
		cerr << FNAME << " Inputs not symmetric!" << endl;
		throw;
	}
	// Compute eigenvalue decomposition for X
	SelfAdjointEigenSolver<MatrixXd> saes(X);
	// Check if X had rank at least n by looking at n-th eigenvalue
	if (saes.eigenvalues().tail(n).minCoeff() < EPS_SVAL) {
		cerr << FNAME << " X has rank less than " << n << "!" << endl;
		throw;
	}
	// Obtain whitening W for X
	MatrixXd V(d,n);
	MatrixXd D(n,n);
	MatrixXd W(d,n);
	V = saes.eigenvectors().rightCols(n);
	D = saes.eigenvalues().tail(n).asDiagonal();
	W = V*D.inverse().cwiseSqrt();
	// Compute residual
	double residual = (MatrixXd::Identity(n,n) - W.transpose()*X*W).norm();
	// Convert to Lmatrix
	Lmatrix Wl(d,n);
	eigen2lmatrix(W,Wl);
	// Contract with Y to obtain Z = Y(W,W,W)
	Index<'i'> i;
	Index<'j'> j;
	Index<'k'> k;
	Index<'a'> a;
	Index<'b'> b;
	Index<'c'> c;
	Z = Ltensor(n,n,n);
	Z(i,j,k) = Y(a,b,c)*Wl(a,i)*Wl(b,j)*Wl(c,k);
	// Compute (W^T)^+
	pinv_trans_W = V*D.cwiseSqrt();
}

