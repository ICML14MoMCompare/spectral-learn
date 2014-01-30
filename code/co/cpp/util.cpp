#include "include/util.h"

// TODO Use a f****ng template!
double nuclear_norm(MatrixXd & m) {
	double z = m.jacobiSvd(ComputeThinU | ComputeThinV).singularValues().sum();
	return z;
}

double nuclear_norm(MatrixXf & m) {
	double z = m.jacobiSvd(ComputeThinU | ComputeThinV).singularValues().sum();
	return z;
}

double fast_nuclear_norm(MatrixXd & m, int R) {
	float z = RedSVD(m.cast<float>(),R).singularValues().sum();
	return (double) z;
}

//#define DEBUG_RANK
#define GLRP_ITERATIONS 100

void generate_low_rank_pnfa(MatrixXd & A, VectorXd & ainf, int rank) {
	A = MatrixXd::Random(A.rows(),A.cols());
	for (int i = 0; i < GLRP_ITERATIONS; i++) {
		fast_nn_projection(A,0.0,rank);
		pnfa_projection(A,ainf);
	}
#ifdef DEBUG_RANK
	cerr.precision(3);
	cerr << "(generate_low_rank_pnfa) Singular values of output:" << endl;
	cerr << A.jacobiSvd(ComputeThinU | ComputeThinV).singularValues().transpose() << endl;
#endif
}

