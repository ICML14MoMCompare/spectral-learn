#include "include/admm.h"

//#define DEBUG_IT
#define DEBUG_CONV
//#define DEBUG_NORM
#define DEBUG_INPUT

// Parameters that measure change in the objective to detect convergence
#define EPSREL 1E-4
#define EPSABS 1E-10

// ADMM cf. Boyd et al., Foundations and Trends in Machine Learning '11
// Convex algorithm for PNFA cf. Balle et al., ICML'12
// In this version the output is guaranteed to be a PNFA
#define FNAME "(opt_icml12_admm)"
void opt_icml12_admm(MatrixXd & H, MatrixXd & H_Sigma, VectorXd & ainf, MatrixXd & A_Sigma, double tau, int maxK) {
	
	// For output of debugging
	cerr.precision(15);

	// Check dimensions (A_Sigma * H == H_Sigma) and ainf.rows == cA
	int pH = H.rows();
	int sH = H.cols();
	int pHs = H_Sigma.rows();
	int sHs = H_Sigma.cols();
	int rA = A_Sigma.rows();
	int cA = A_Sigma.cols();
	int lainf = ainf.rows();
	if (cA != pH or rA != pHs or sH != sHs or lainf != cA) {
		cerr << FNAME << " ERROR: Dimensions mismatch!" << endl;
		return;
	}
	
#ifdef DEBUG_INPUT
	cerr << FNAME << " Called with input:" << endl;
	cerr << "::::: H\t(" << pH << "x" << sH << ")" << endl;
	cerr << "::::: H_Sigma\t(" << pHs << "x" << sHs << ")" << endl;
	cerr << "::::: A_Sigma\t(" << rA << "x" << cA << ")" << endl;
	cerr << "::::: ainf\t(" << lainf << ")" << endl;
	cerr << "::::: tau = " << tau << endl;
	cerr << "::::: maxK = " << maxK << endl;
#endif

	// Number of terms in the ADMM objective split
	const unsigned int NT = 3;

	// Placeholders for primal and dual variables in ADMM
	vector<MatrixXd> X = vector<MatrixXd>(NT, MatrixXd::Zero(rA, cA));
	MatrixXd Z = MatrixXd::Zero(rA, cA);
	MatrixXd Zprev = MatrixXd::Zero(rA, cA);
	vector<MatrixXd> Y = vector<MatrixXd>(NT, MatrixXd::Zero(rA, cA));

	// For learning PNFA, a good initial point is a low rank stochastic
	// matrix
	Z = MatrixXd::Ones(rA, cA);
	pnfa_projection(Z,ainf);

	// Parameters controlling convergence
	double rho = 1.0;
	double r2 = 0.0;
	double s2 = 0.0;
	double epspri = 0.0;
	double epsdual = 0.0;

	// Precomputations for proximal operator of
	// loss function ||X*H - H_Sigma||_F
	MatrixXd Aux1 = 2*H_Sigma*H.transpose();
	MatrixXd Aux2 = 2*H*H.transpose();
	MatrixXd Aux3 = MatrixXd::Identity(pH,pH);

	int k = 0;
	do {
		// Proximal of loss
		// TODO Compute without explicit inverse
		X[0] = (Aux1 + rho*(Z - Y[0]))*((Aux2 + rho*Aux3).inverse());
		// Proximal of trace norm
		X[1] = Z - Y[1];
		nn_projection(X[1], tau/rho);
		// Proximal of PNFA constraints
		X[2] = Z - Y[2];
		pnfa_projection(X[2], ainf);

		// Update aggregate and residuals
		Z = MatrixXd::Zero(rA, cA);
		for (int i = 0; i < NT; i++) Z += X[i] + Y[i];
		Z /= NT;
		for (int i = 0; i < NT; i++) Y[i] += X[i] - Z;
		
		// Update r2 and s2
		r2 = 0.0;
		for (int i = 0; i < NT; i++) r2 += (X[i] - Z).norm();
		r2 = sqrt(r2);
		s2 = 4*rho*rho*(Z - Zprev).norm();
		s2 = sqrt(s2);
		
		// Update rho
		if (r2 > 10*s2) rho = 2*rho;
		else if (s2 > 10*r2) rho = rho/2;

		// Save value of aggregate
		Zprev = Z;
		
		// Update stopping criterion
		double maxXZ = 0.0;
		for (int i = 0; i < NT; i++) if (X[i].squaredNorm() > maxXZ) maxXZ = X[i].squaredNorm();
		if (Z.squaredNorm() > maxXZ) maxXZ = Z.squaredNorm();
		double maxY = 0.0;
		for (int i = 0; i < NT; i++) if (Y[i].squaredNorm() > maxY) maxY = Y[i].squaredNorm();
		epspri = sqrt(rA*cA)*EPSABS + EPSREL*maxXZ;
		epsdual = sqrt(rA*cA)*EPSABS + EPSREL*maxY;
		
#ifdef DEBUG_IT
		cerr << FNAME << " Iteration: " << k  << "\tObjectives: " << r2 << ' ' << s2 << endl;
#endif
		
		k++;

	} while (r2 > epspri and s2 > epsdual and k < maxK);

#ifdef DEBUG_CONV
	cerr << FNAME << " Converged after " << k << " iterations ";
	if (k >= maxK) cerr << "by MAXIMUM ITERATIONS";
	else if (r2 <= epspri) cerr << "by PRIMAL CRITERION";
	else if (s2 <= epsdual) cerr << "by DUAL CRITERION";
	cerr << endl;
#endif

	// TODO Recompute ainf such that the probabilities add up to 1
	// (which is not guaranteed because Z is not projected onto simplex)
	// Question: project A_Sigma, OR adjust ainf (could be A_Sigma is not
	// all positive...)
	// Here we project Z again onto the simplex
	// Note that this might increase the rank :(
#ifdef DEBUG_NORM
	cerr << FNAME << " Output of optimization, before normalization (ainf = " << ainf.transpose() << ")" << endl;
	cerr << Z << endl;
#endif
	pnfa_projection(Z, ainf);
#ifdef DEBUG_NORM
	cerr << FNAME << " After normalization" << endl;
	cerr << Z << endl;
#endif

	// Assign to output
	A_Sigma = Z;
}

// ADMM cf. Boyd et al., Foundations and Trends in Machine Learning '11
// Convex algorithm for WA cf. Balle et al., ICML'12
// This version outputs a WA (not necessarily probabilistic)
#define FNAME "(opt_icml12_admm_wa)"
void opt_icml12_admm_wa(MatrixXd & H, MatrixXd & H_Sigma, MatrixXd & A_Sigma, double tau, int maxK) {
	
	// For output of debugging
	cerr.precision(15);

	// Check dimensions (A_Sigma * H == H_Sigma)
	int pH = H.rows();
	int sH = H.cols();
	int pHs = H_Sigma.rows();
	int sHs = H_Sigma.cols();
	int rA = A_Sigma.rows();
	int cA = A_Sigma.cols();
	if (cA != pH or rA != pHs or sH != sHs) {
		cerr << "(opt_icml12_admm) ERROR: Dimensions mismatch!" << endl;
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
	
	// Number of terms in the ADMM objective split
	const unsigned int NT = 2;

	// Placeholders for primal and dual variables in ADMM
	vector<MatrixXd> X = vector<MatrixXd>(NT, MatrixXd::Zero(rA, cA));
	MatrixXd Z = MatrixXd::Zero(rA, cA);
	MatrixXd Zprev = MatrixXd::Zero(rA, cA);
	vector<MatrixXd> Y = vector<MatrixXd>(NT, MatrixXd::Zero(rA, cA));
	
	// For learning WA, a good initial point is the unregularized
	// least-squares solution
	// BUT: doesn't seem to improve convergence
	JacobiSVD<MatrixXd> Hsvd(H, ComputeThinU | ComputeThinV);
	Z = H_Sigma*Hsvd.matrixU()*Hsvd.singularValues().asDiagonal().inverse()*Hsvd.matrixV().transpose();

	// Parameters controlling convergence
	double rho = 1.0;
	double r2 = 0.0;
	double s2 = 0.0;
	double epspri = 0.0;
	double epsdual = 0.0;

	// Precomputations for proximal operator of
	// loss function ||X*H - H_Sigma||_F
	MatrixXd Aux1 = 2*H_Sigma*H.transpose();
	MatrixXd Aux2 = 2*H*H.transpose();
	MatrixXd Aux3 = MatrixXd::Identity(pH,pH);

	int k = 0;
	do {
		// Proximal of loss
		// TODO Compute without explicit inverse
		X[0] = (Aux1 + rho*(Z - Y[0]))*((Aux2 + rho*Aux3).inverse());
		// Proximal of trace norm
		X[1] = Z - Y[1];
		nn_projection(X[1], tau/rho);

		// Update aggregate and residuals
		Z = MatrixXd::Zero(rA, cA);
		for (int i = 0; i < NT; i++) Z += X[i] + Y[i];
		Z /= NT;
		for (int i = 0; i < NT; i++) Y[i] += X[i] - Z;
		
		// Update r2 and s2
		r2 = 0.0;
		for (int i = 0; i < NT; i++) r2 += (X[i] - Z).norm();
		r2 = sqrt(r2);
		s2 = 4*rho*rho*(Z - Zprev).norm();
		s2 = sqrt(s2);
		
		// Update rho
		if (r2 > 10*s2) rho = 2*rho;
		else if (s2 > 10*r2) rho = rho/2;

		// Save value of aggregate
		Zprev = Z;
		
		// Update stopping criterion
		double maxXZ = 0.0;
		for (int i = 0; i < NT; i++) if (X[i].squaredNorm() > maxXZ) maxXZ = X[i].squaredNorm();
		if (Z.squaredNorm() > maxXZ) maxXZ = Z.squaredNorm();
		double maxY = 0.0;
		for (int i = 0; i < NT; i++) if (Y[i].squaredNorm() > maxY) maxY = Y[i].squaredNorm();
		epspri = sqrt(rA*cA)*EPSABS + EPSREL*maxXZ;
		epsdual = sqrt(rA*cA)*EPSABS + EPSREL*maxY;
		
#ifdef DEBUG_IT
		cerr << FNAME << " Iteration: " << k  << "\tObjectives: " << r2 << ' ' << s2 << endl;
#endif
		
		k++;

	} while (r2 > epspri and s2 > epsdual and k < maxK);

#ifdef DEBUG_CONV
	cerr << FNAME << " Converged after " << k << " iterations ";
	if (k >= maxK) cerr << "by MAXIMUM ITERATIONS";
	else if (r2 <= epspri) cerr << "by PRIMAL CRITERION";
	else if (s2 <= epsdual) cerr << "by DUAL CRITERION";
	cerr << endl;
#endif

	// Assign to output
	A_Sigma = Z;
}

