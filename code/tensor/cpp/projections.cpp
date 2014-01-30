#include "include/projections.h"

//#define DEBUG_SP
//#define DEBUG_PP

// Tolerance when checking non-zero in different projections
#define TOL_ZERO_SIMPLEX 1E-8
#define TOL_ZERO_NUCLEAR 1E-10

bool sort_func(double a, double b) { return (a > b); }

// Cf. Duchi et al., ICML '08
bool simplex_projection(VectorXd & v, double norm) {
#ifdef DEBUG_SP
	cerr << "(simplex_projection) Projecting to norm: " << norm << endl;
#endif
	// w = sort(v), decreasing
	// TODO Alternative, define a RandomIterator for RowVectorXd
	int i;
	vector<double> w;
	w.clear();
	for (i = 0; i < v.rows(); i++) w.push_back(v[i]);
	sort(w.begin(), w.end(), sort_func);
	// Find number of positive coordinates
	double cumw = 0.0;
	for (i = 0; i < w.size(); i++) {
		cumw += w[i];
#ifdef DEBUG_SP
		cerr << "(simplex_projection) cumw: " << cumw << "\tw: " << w[i] << "\tthr: ";
		cerr << (cumw - norm)/(double) (i+1) << endl;
#endif
		if (w[i] - ((cumw - norm)/(double) (i+1)) <= -TOL_ZERO_SIMPLEX) {
			// Adjust cumw for theta computation
			cumw -= w[i];
			break;
		}
	}
	// Threshold coordinates
#ifdef DEBUG_SP
	//int rho = min(i,(int) w.size());
	int rho = i;
	int nonzero = 0;
#endif
	double theta = max(0.0,(cumw - norm)/(double) i);
	for (i = 0; i < v.rows(); i++) {
	       v[i] = max(0.0,v[i] - theta);
#ifdef DEBUG_SP
	       if (v[i] > 0.0) nonzero++;
#endif
	}
	// If sum is less than intended norm, move along the (1,1,...,1) vector until touching simplex
	if (v.sum() < norm) {
		// Compute slack
		double s = norm - v.sum();
		// Spread it across components
		s = s / (double) v.rows();
		for (i = 0; i < v.rows(); i++) v[i] += s;
	}
#ifdef DEBUG_SP
	cerr << "(simplex_projection) Rho: " << rho << ", Non-zero: " << nonzero << ", Sum: " << v.sum() << endl;
#endif
	return true;
}

// Cf. definition in Balle, UPC PhD thesis '13
bool pnfa_projection(MatrixXd & X, const VectorXd & v) {
	// Check dimensions, X is a x b and v is c x 1
	// With a = p*k, b = p, c = p
	// Where k = # symbols, p = # prefixes
	int pk = X.rows();
	int p = X.cols();
	if (p != v.rows() or (pk % p) != 0) {
		cerr << "(pnfa_projection) ERROR: Dimensions mismatch!" << endl;
		throw;
	}
	// Compute number of operators
	int k = pk/p;
	// Projections of all rows corresponding to same prefix together
	// pi == prefix index
	// si == symbol index
	// ri == row index in X
	// ci == columns index
	// i == index inside w
	for (int pi = 0; pi < p; pi++) {
		double z = 1.0 - v[pi];
		VectorXd w(pk);
		int i = 0;
		for (int si = 0; si < k; si++) {
			int ri = si*p + pi;
			for (int ci = 0; ci < p; ci++) {
				w(i) = X(ri,ci);
				i++;
			}
		}
#ifdef DEBUG_PP
		cerr << "(pnfa_projection) Row " << i << " BEFORE projection:";
		cerr << endl << w << endl;
#endif
		simplex_projection(w, z);
#ifdef DEBUG_PP
		cerr << "(pnfa_projection) Row " << i << " AFTER projection:";
		cerr << endl << w << endl;
#endif
		// Rewrite everything to its place
		i = 0;
		for (int si = 0; si < k; si++) {
			int ri = si*p + pi;
			for (int ci = 0; ci < p; ci++) {
				X(ri,ci) = w(i);
				i++;
			}
		}
	}
	return true;
}


// Cf. definition in Balle, UPC PhD thesis '13
bool probmat_projection(MatrixXd & X, VectorXd & a) {
	int r = X.rows();
	int c = X.cols();

	if(r != a.rows())	{
		cerr << "(probmat_projection) ERROR: Dimensions mismatch!" << endl;
		throw;	
	}

	for (int ri = 0; ri < r; ri++) {
		VectorXd w(c);
		for(int ci = 0; ci < c; ci++){
			w(ci) = X(ri,ci);
		}
		simplex_projection(w,1.0-a(ri));
		for(int ci = 0; ci < c; ci++){
			X(ri,ci) = w(ci);
		}
	}
	return true;
}

/*
// Cf. Cai et al., SIAM J. Opt. '08
bool nn_projection(MatrixXd & X, double thr) {
	int r = X.rows();
	int c = X.cols();
	int n = min(r, c);

	JacobiSVD<MatrixXd> svd(X, ComputeThinU | ComputeThinV);
	const VectorXd & svv = svd.singularValues();
	if (svv.lpNorm<1>() < TOL_ZERO_NUCLEAR) return false;

	VectorXd S = VectorXd::Zero(n);
	for (int i = 0; i < n; i++) S(i) = svv(i) > thr ? svv(i) - thr : 0;
	X = svd.matrixU()*S.asDiagonal()*svd.matrixV().transpose();

	return true;
}

// Cf. Cai et al., SIAM J. Opt. '08
// With truncated SVD by Halko
bool fast_nn_projection(MatrixXd & X, double thr, int R) {
	// Figure out the dimensions of things
	// Will project to 2*R (if we can) to copmute SVD
	// Then take the first R
	int mrc = min(X.rows(), X.cols());
	int outrank = min(R,mrc);
	int svdrank = min(2*outrank,mrc);

	// Cast because RedSVD works with floats
	MatrixXf Xf = X.cast<float>();

	// Compute reduced SVD
	RedSVD svd(Xf,svdrank);

	// Check if matrix is already "0"
	VectorXf svv = svd.singularValues().head(outrank);
	if (svv.lpNorm<1>() < TOL_ZERO_NUCLEAR) return false;

	// Truncated necessary singular values
	VectorXf S = VectorXf::Zero(outrank);
	for (int i = 0; i < outrank; i++) S(i) = svv(i) > thr ? svv(i) - thr : 0;

	// Compute output (and cast back to double)
	Xf = svd.matrixU().leftCols(outrank)*S.asDiagonal()*svd.matrixV().leftCols(outrank).transpose();
	X = Xf.cast<double>();

	return true;
}
*/
