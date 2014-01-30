#include "include/tensor_symmetrize.h"

void computeNfromLowRankInverse(MatrixXd & H, int n, MatrixXd & N) {
	int r = H.rows();
	int c = H.cols();
	JacobiSVD<MatrixXd> svd(H, ComputeThinU | ComputeThinV);
	MatrixXd Dinv(n,n);
	Dinv = svd.singularValues().head(n).asDiagonal().inverse();
	N = svd.matrixV().leftCols(n)*Dinv*svd.matrixU().leftCols(n).transpose();
}

