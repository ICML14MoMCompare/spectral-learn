#include "include/frank_wolfe.h"

//#define DEBUG_IT
#define DEBUG_CONV
#define DEBUG_INPUT

// Parameters that measure change in the objective to detect convergence
#define EPSREL 1E-4
#define EPSABS 1E-10

#define LINESEARCH

// Optimizing f(X) subject to ||X||_* <= alpha (with Frank-Wolfe type optimization)
// We like to work with doubles (but RedSVD uses floats!)
void opt_icml12_fw_wa(MatrixXd & H, MatrixXd & H_Sigma, MatrixXd & A_Sigma, double alpha, int maxK) {
	MatrixXf Hf = H.cast<float>();
	MatrixXf H_Sigmaf = H_Sigma.cast<float>();
	MatrixXf A_Sigmaf = A_Sigma.cast<float>();
	opt_icml12_fw_wa(Hf, H_Sigmaf, A_Sigmaf, alpha, maxK);
	A_Sigma = A_Sigmaf.cast<double>();
}

#define FNAME "(opt_icml12_fw_wa)"
void opt_icml12_fw_wa(MatrixXf & H, MatrixXf & H_Sigma, MatrixXf & A_Sigma, double alpha, int maxK) {
	
	// For output of debugging
	cerr.precision(15);

	// Check dimensions (A_Sigma * H == H_Sigma) and ainf.rows == cA
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
	cerr << "::::: alpha = " << alpha << endl;
	cerr << "::::: maxK = " << maxK << endl;
#endif
	
	// Placeholders for variables
	MatrixXf X = MatrixXf::Zero(rA, cA);
	MatrixXf Xprev = MatrixXf::Zero(rA, cA);
	MatrixXf Y = MatrixXf::Zero(rA, cA);
	MatrixXf G = MatrixXf::Zero(rA, cA);
	RedSVD D;

	// Size of updates (using "standard" from Jaggi '12)
	double gamma = 0.0;
	double gammak = 0.0;

	// Parameters controlling convergence
	double r2 = 0.0;
	double epspri = 0.0;

	// Precomputations for gradient of
	// loss function ||X*H - H_Sigma||_F
	MatrixXf Aux1 = 2.0*H_Sigma*H.transpose();
	MatrixXf Aux2 = 2.0*H*H.transpose();

	int k = 0;
	do {
		Xprev = X;
		
		// Compute the gradient at the current point
		G = X*Aux2 - Aux1;
		// Obtain the best rank 1 approximation
		D.run(G,1);
		// Compute direction of update
		Y = -alpha*D.matrixU()*(D.matrixV().transpose());
		// Compute ammount of update
		gammak = 2.0/(2.0 + (double) k);
#ifdef LINESEARCH
		MatrixXf E = Y - X;
		gamma = -(X*H-H_Sigma).cwiseProduct(E*H).sum();
		gamma = gamma / (E*H).squaredNorm();
		if (gamma > gammak) gammak = gamma;
		if (gammak > 1.0) gammak = 1.0;
#endif
		// Update current point
		X = (1.0-gammak)*X + gammak*Y;

		// Update stopping criterion
		r2 = (Xprev - X).norm();
		epspri = sqrt(rA*cA)*EPSABS + EPSREL*X.norm();
		
#ifdef DEBUG_IT
		cerr << FNAME << " Iteration: " << k;
		cerr << " (gammak: " << gammak << ")\t";
		cerr << "Change: " << r2 << endl;
#endif
		
		k++;

	} while (r2 > epspri and k < maxK);

#ifdef DEBUG_CONV
	cerr << FNAME << " Converged after " << k << " iterations ";
	if (k >= maxK) cerr << "by MAXIMUM ITERATIONS";
	else cerr << "by CONVERGENCE CRITERION";
	cerr << endl;
#endif

	// Assign to output
	A_Sigma = X;
}

// Optimizing f(X) + tau ||X||_* (with Frank-Wolfe type optimization)
// (Cf. Harchaoui et al. '13)
// We like to work with doubles (but RedSVD uses floats!)
void fw_reg_wa(MatrixXd & H, MatrixXd & H_Sigma, MatrixXd & A_Sigma, double tau, int maxK) {
	MatrixXf Hf = H.cast<float>();
	MatrixXf H_Sigmaf = H_Sigma.cast<float>();
	MatrixXf A_Sigmaf = A_Sigma.cast<float>();
	opt_icml12_fw_wa(Hf, H_Sigmaf, A_Sigmaf, tau, maxK);
	A_Sigma = A_Sigmaf.cast<double>();
}

#define FNAME "(fw_reg_wa)"
void fw_reg_wa(MatrixXf & H, MatrixXf & H_Sigma, MatrixXf & A_Sigma, double tau, int maxK) {
	
	// For output of debugging
	cerr.precision(6);

	// Check dimensions (A_Sigma * H == H_Sigma) and ainf.rows == cA
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
	cerr << FNAME << " Called with input: ";
	cerr << "H (" << pH << "x" << sH << ") ";
	cerr << "H_Sigma (" << pHs << "x" << sHs << ") ";
	cerr << "A_Sigma (" << rA << "x" << cA << ") ";
	cerr << "tau = " << tau;
	cerr << " maxK = " << maxK;
	cerr << endl;
#endif
	
	// Placeholders for variables
	MatrixXf X = MatrixXf::Zero(rA, cA);
	MatrixXf Xprev = MatrixXf::Zero(rA, cA);
	MatrixXf Y = MatrixXf::Zero(rA, cA);
	MatrixXf G = MatrixXf::Zero(rA, cA);
	double nuc = 0.0;
	double nucPrev = 0.0;
	RedSVD D;

	// Size of updates
	double gamma = 0.0;

	// Upper bound on nuclear norm
	double nucBound = H_Sigma.squaredNorm() / tau;

	// Parameters controlling convergence
	double r2 = 0.0;
	double epspri = 0.0;

	// Precomputations for gradient of
	// loss function ||X*H - H_Sigma||_F
	MatrixXf Aux1 = 2.0*H_Sigma*H.transpose();
	MatrixXf Aux2 = 2.0*H*H.transpose();

	// TODO Use interior-point method to optimize update
	// over a fixed number of previous points (e.g. 5)
	// This should give a big improvement on convergence

	int k = 0;
	do {
		Xprev = X;
		nucPrev = nuc;
		
		// Compute the gradient at the current point
		G = X*Aux2 - Aux1;
		// Obtain the best rank 1 approximation
		D.run(G,1);
		// Compute direction of update
		Y = -nucBound*D.matrixU()*(D.matrixV().transpose());
		// Compute ammount of update
		MatrixXf E = Y - X;
		gamma = -(tau*(nucBound-nuc) + 2.0*(X*H-H_Sigma).cwiseProduct(E*H).sum());
		gamma = gamma / (2.0*(E*H).squaredNorm());
		if (gamma < 0.0  ) gamma = 0.0;
		if (gamma > 1.0) gamma = 1.0;
		// Update current point
		X = (1.0-gamma)*X + gamma*Y;
		nuc = (1.0-gamma)*nuc + gamma*nucBound;

		// Update stopping criterion
		r2 = (Xprev - X).norm() + abs(nucPrev - nuc);
		epspri = sqrt(rA*cA)*EPSABS + EPSREL*(X.norm() + nuc);
		
#ifdef DEBUG_IT
		cerr << FNAME << " Iteration: " << k;
		cerr << " (gamma: " << gamma << ")\t";
		cerr << "Change: " << r2;
		cerr << "\tNuclear: " << nuc << endl;
#endif
		
		k++;

	} while (r2 > epspri and k < maxK);

#ifdef DEBUG_CONV
	cerr << FNAME << " Converged after " << k << " iterations ";
	if (k >= maxK) cerr << "by MAXIMUM ITERATIONS";
	else cerr << "by CONVERGENCE CRITERION";
	cerr << endl;
#endif

	// Assign to output
	A_Sigma = X;
}

