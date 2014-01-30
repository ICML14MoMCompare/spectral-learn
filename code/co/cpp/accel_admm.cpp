#include "include/admm.h"

//#define DEBUG_IT
#define DEBUG_CONV
//#define DEBUG_NORM
#define DEBUG_INPUT

// Parameters that measure change in the objective to detect convergence
#define EPSREL 1E-4
#define EPSABS 1E-10
// Parameters controlling amount and complexity of expensive iterations
// in ADMM using truncated randomized SVD
#define EXPENSIVE_RATIO 0.1
#define RANK_EXPENSE_FACTOR 3

// A (hopefully) accelerated version of ADMM for 3 terms objectives
// Without step correction
#define FNAME "(accel_admm)"
void accel_admm(MatrixXd & H, MatrixXd & H_Sigma, VectorXd & ainf, MatrixXd & A_Sigma, double tau, int maxK, bool Ainitial) {
	
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
	cerr << FNAME << " Called with input: ";
	cerr << "H (" << pH << "x" << sH << ") ";
	cerr << "H_Sigma (" << pHs << "x" << sHs << ") ";
	cerr << "A_Sigma (" << rA << "x" << cA << ") ";
	if (Ainitial) cerr << "[initial] ";
	cerr << "ainf (" << lainf << ") ";
	cerr << "tau = " << tau;
	cerr << " maxK = " << maxK;
	cerr << endl;
#endif
	
	// Number of terms in the ADMM objective split
	const unsigned int NT = 3;

	// Placeholders for primal and dual variables in ADMM
	vector<MatrixXd> X = vector<MatrixXd>(NT, MatrixXd::Zero(rA, cA));
	vector<MatrixXd> U = vector<MatrixXd>(NT, MatrixXd::Zero(rA, cA));
	MatrixXd Y = MatrixXd::Zero(rA, cA);
	MatrixXd Yprev = MatrixXd::Zero(rA, cA);

	// For learning PNFA, a good initial point is a low rank stochastic
	// matrix
	MatrixXd Z;
	if (!Ainitial) {
		Z = MatrixXd::Ones(rA, cA);
		pnfa_projection(Z,ainf);
		Y = Z;
	} else {
		Y = A_Sigma;
	}

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

		// Update primal variables
		Yprev = Y;
		// Proximal of loss
		X[0] = (Aux1 + rho*(Y - U[0]))*((Aux2 + rho*Aux3).inverse());
		// Proximal of trace norm
		X[1] = X[0] - U[1];
		nn_projection(X[1], tau/rho);
		//fast_nn_projection(X[1], tau/rho, maxRank);
		// Proximal of PNFA constraints
		X[2] = X[1] - U[2];
		pnfa_projection(X[2], ainf);

		// Obtain average
		Y = (X[0] + X[1] + X[2]) / 3.0;
		Y += (U[0] + U[1] + U[2]) / 3.0;

		// Update dual variables
		U[0] = U[0] + (X[0] - Y);
		U[1] = U[1] + (X[1] - Y);
		U[2] = U[2] + (X[2] - Y);

		// Compute residuals
		r2 = 0.0;
		r2 += (X[0]-Y).squaredNorm();
		r2 += (X[1]-Y).squaredNorm();
		r2 += (X[2]-Y).squaredNorm();
		r2 = sqrt(r2);
		s2 = rho*(Y-Yprev).norm();
		
		// Update rho
		if (r2 > 10*s2) rho = 2*rho;
		else if (s2 > 10*r2) rho = rho/2;

		// Update stopping criterion
		double maxU = 0.0;
		for (int i = 0; i < NT; i++) if (U[i].squaredNorm() > maxU) maxU = U[i].squaredNorm();
		epspri = sqrt(rA*cA)*EPSABS + EPSREL*Y.squaredNorm();
		epsdual = sqrt(rA*cA)*EPSABS + EPSREL*maxU;
		
#ifdef DEBUG_IT
		cerr << FNAME << " Iteration: " << k;
		cerr << "(rho: " << rho << ")\t";
		cerr << "Objectives: " << r2 << ' ' << s2 << endl;
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


	// Here we project Z again onto the simplex
	// Note that this might increase the rank :(
	Z = X[2];
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

// A (hopefully) accelerated version of ADMM for 3 terms objectives
// With truncated SVD
// With a budget of "expensive steps"
// And backtracking (at the end)
#define FNAME "(accel_trunc_admm)"
void accel_trunc_admm(MatrixXd & H, MatrixXd & H_Sigma, VectorXd & ainf, MatrixXd & A_Sigma, double tau, int maxK, int maxRank) {
	
	// For output of debugging
	cerr.precision(8);

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
	cerr << FNAME << " Called with input: ";
	cerr << "H (" << pH << "x" << sH << ") ";
	cerr << "H_Sigma (" << pHs << "x" << sHs << ") ";
	cerr << "A_Sigma (" << rA << "x" << cA << ") ";
	cerr << "ainf (" << lainf << ") ";
	cerr << "tau = " << tau;
	cerr << " maxK = " << maxK;
	cerr << " maxRank = " << maxRank << endl;
#endif
	
	// Number of terms in the ADMM objective split
	const unsigned int NT = 3;

	// Placeholders for primal and dual variables in ADMM
	vector<MatrixXd> X = vector<MatrixXd>(NT, MatrixXd::Zero(rA, cA));
	vector<MatrixXd> U = vector<MatrixXd>(NT, MatrixXd::Zero(rA, cA));
	MatrixXd Y = MatrixXd::Zero(rA, cA);
	MatrixXd Yprev = MatrixXd::Zero(rA, cA);

	// For learning PNFA, a good initial point is a low rank stochastic
	// matrix
	MatrixXd Z;
       	Z = MatrixXd::Ones(rA, cA);
	pnfa_projection(Z,ainf);
	Y = Z;

	// Parameters controlling convergence
	double rho = 1.0;
	double r2 = 0.0;
	double s2 = 0.0;
	double epspri = 0.0;
	double epsdual = 0.0;

	// Parameters controlling SVD truncation
	int svdTrunc = (int) ceil(1.5 * (double) maxRank);
	//int svdExpensive = 10*svdTrunc;
	int svdExpensive = RANK_EXPENSE_FACTOR*svdTrunc;
	int mraca = min(rA,cA);
	int remainingExpensive = (int) ceil(EXPENSIVE_RATIO * (double) maxK);
	bool expensive = false;

#ifdef DEBUG_INPUT
	cerr << FNAME << " SVD SETUP --- svdTrunc: " << svdTrunc;
	cerr << " svdExpensive: " << svdExpensive;
	cerr << " (full? " << (svdExpensive >= mraca) << ")";
	cerr << " expensiveIterations: " << remainingExpensive << endl;
#endif

	// Parameters controlling "backtracking"
	double objective = H_Sigma.squaredNorm();
	double objectiveBest = objective;
	MatrixXd Ybest = MatrixXd::Zero(rA, cA);

	// Precomputations for proximal operator of
	// loss function ||X*H - H_Sigma||_F
	MatrixXd Aux1 = 2*H_Sigma*H.transpose();
	MatrixXd Aux2 = 2*H*H.transpose();
	MatrixXd Aux3 = MatrixXd::Identity(pH,pH);

	int k = 0;
	do {
		// Save variables from previous iteration
		Yprev = Y;

		// Update primal variables
		// Proximal of loss
		X[0] = (Aux1 + rho*(Y - U[0]))*((Aux2 + rho*Aux3).inverse());
		// Proximal of trace norm
		X[1] = X[0] - U[1];
		// Check whether we want expensive projection of not
		if (expensive) {
			expensive = false;
			remainingExpensive--;
			if (svdExpensive < mraca) {
				fast_nn_projection(X[1], tau/rho, svdExpensive);
			} else {
				nn_projection(X[1], tau/rho);
			}
		} else {
			fast_nn_projection(X[1], tau/rho, svdTrunc);
		}

		// Proximal of PNFA constraints
		X[2] = X[1] - U[2];
		pnfa_projection(X[2], ainf);

		// Obtain new point as average
		Y = (X[0] + X[1] + X[2]) / 3.0;
		Y += (U[0] + U[1] + U[2]) / 3.0;

		// Update dual variables
		U[0] = U[0] + (X[0] - Y);
		U[1] = U[1] + (X[1] - Y);
		U[2] = U[2] + (X[2] - Y);

		// If objective is not improving (surely because of approximated SVD)
		// Then schedule an expensive iteration
		objective = (Y*H-H_Sigma).squaredNorm();
		if (objective > objectiveBest) {
			// Backtracking??
			//objective = objectiveBest;
			//Y = Ybest;
			//for (int i = 0; i < NT; i++) U[i] = Ubest[i];
			expensive = true;
		} else {
			// Save current point for backtracking
			objectiveBest = objective;
			Ybest = Y;
			//for (int i = 0; i < NT; i++) Ubest[i] = U[i];
		}

		// Compute residuals (to measure convergence)
		r2 = 0.0;
		r2 += (X[0]-Y).squaredNorm();
		r2 += (X[1]-Y).squaredNorm();
		r2 += (X[2]-Y).squaredNorm();
		r2 = sqrt(r2);
		s2 = rho*(Y-Yprev).norm();
		
		// Update rho
		if (r2 > 10*s2) rho = 2*rho;
		else if (s2 > 10*r2) rho = rho/2;

		// Update stopping criterion
		double maxU = 0.0;
		for (int i = 0; i < NT; i++) if (U[i].squaredNorm() > maxU) maxU = U[i].squaredNorm();
		epspri = sqrt(rA*cA)*EPSABS + EPSREL*Y.squaredNorm();
		epsdual = sqrt(rA*cA)*EPSABS + EPSREL*maxU;

#ifdef DEBUG_IT
		cerr << FNAME << " Iteration: " << k;
		cerr << "(rho: " << rho << " remExpense: " << remainingExpensive << ")\t";
		cerr << "Objectives: " << r2 << ' ' << s2 << ' ' << objective << endl;
#endif
		
		k++;

	} while (r2 > epspri and s2 > epsdual and k < maxK and (!expensive or remainingExpensive > 0));

#ifdef DEBUG_CONV
	cerr << FNAME << " Converged after " << k << " iterations ";
	if (k >= maxK) cerr << "by MAXIMUM ITERATIONS";
	else if (expensive && remainingExpensive <= 0) cerr << "by MAXIMUM EXPENSIVE ITERATIONS";
	else if (r2 <= epspri) cerr << "by PRIMAL CRITERION";
	else if (s2 <= epsdual) cerr << "by DUAL CRITERION";
	cerr << endl;
#endif

	// Here we project Z again onto the simplex
	// Note that this might increase the rank :(
	//Z = X[2];
	//Z = Y;
	Z = Ybest;
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

// A (hopefully) accelerated version of ADMM for WA
// With truncated SVD
// With a budget of "expensive steps" (uses same ratio as for PNFA)
// And backtracking (at the end)
// Possible with over-relaxation
#define OVER_RELAXED
#define OR_ALPHA 1.5
#define FNAME "(accel_trunc_admm_wa)"
void accel_trunc_admm_wa(MatrixXd & H, MatrixXd & H_Sigma, MatrixXd & A_Sigma, double tau, int maxK, int maxRank) {
	
	// For output of debugging
	cerr.precision(8);

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
	cerr << " maxRank = " << maxRank << endl;
#endif
	
	// Number of terms in the ADMM objective split
	const unsigned int NT = 2;

	// Placeholders for primal and dual variables in ADMM
	MatrixXd X = MatrixXd::Zero(rA, cA);
	MatrixXd Y = MatrixXd::Zero(rA, cA);
	MatrixXd Yprev = MatrixXd::Zero(rA, cA);
	MatrixXd U = MatrixXd::Zero(rA, cA);

	// Parameters controlling convergence
	double rho = 1.0;
	double r2 = 0.0;
	double s2 = 0.0;
	double epspri = 0.0;
	double epsdual = 0.0;

	// Parameters controlling SVD truncation
	int svdTrunc = (int) ceil(1.5 * (double) maxRank);
	//int svdExpensive = 10*svdTrunc;
	int svdExpensive = RANK_EXPENSE_FACTOR*svdTrunc;
	int mraca = min(rA,cA);
	int remainingExpensive = (int) ceil(EXPENSIVE_RATIO * (double) maxK);
	bool expensive = false;

#ifdef DEBUG_INPUT
	cerr << FNAME << " SVD SETUP --- svdTrunc: " << svdTrunc;
	cerr << " svdExpensive: " << svdExpensive;
	cerr << " (full? " << (svdExpensive >= mraca) << ")";
	cerr << " expensiveIterations: " << remainingExpensive << endl;
#endif

	// Parameters controlling "backtracking"
	double objective = H_Sigma.squaredNorm();
	double objectiveBest = objective;
	MatrixXd Ybest = MatrixXd::Zero(rA, cA);

	// Precomputations for proximal operator of
	// loss function ||X*H - H_Sigma||_F
	MatrixXd Aux1 = 2*H_Sigma*H.transpose();
	MatrixXd Aux2 = 2*H*H.transpose();
	MatrixXd Aux3 = MatrixXd::Identity(pH,pH);

	int k = 0;
	do {
		// Save variables from previous iteration
		Yprev = Y;

		// Update primal variables
		// Proximal of loss
		X = (Aux1 + rho*(Y - U))*((Aux2 + rho*Aux3).inverse());
		// Proximal of trace norm
#ifdef OVER_RELAXED
		Y = OR_ALPHA*X + (1.0-OR_ALPHA)*Yprev + U;
#else
		Y = X + U;
#endif
		// Check whether we want expensive projection of not
		if (expensive) {
			expensive = false;
			remainingExpensive--;
			if (svdExpensive < mraca) {
				fast_nn_projection(Y, tau/rho, svdExpensive);
			} else {
				nn_projection(Y, tau/rho);
			}
		} else {
			fast_nn_projection(Y, tau/rho, svdTrunc);
		}

		// Update dual variables
#ifdef OVER_RELAXED
		U = U + (OR_ALPHA*X + (1.0-OR_ALPHA)*Yprev - Y);
#else
		U = U + (X - Y);
#endif

		// If objective is not improving (surely because of approximated SVD)
		// Then schedule an expensive iteration
		objective = (Y*H-H_Sigma).squaredNorm();
		if (objective > objectiveBest) {
			// Backtracking??
			//objective = objectiveBest;
			//Y = Ybest;
			//for (int i = 0; i < NT; i++) U[i] = Ubest[i];
			expensive = true;
		} else {
			// Save current point for backtracking
			objectiveBest = objective;
			Ybest = Y;
			//for (int i = 0; i < NT; i++) Ubest[i] = U[i];
		}

		// Compute residuals (to measure convergence)
		r2 = (X-Y).norm();
		s2 = rho*(Y-Yprev).norm();
		
		// Update rho
		if (r2 > 10*s2) rho = 2*rho;
		else if (s2 > 10*r2) rho = rho/2;

		// Update stopping criterion
		epspri = sqrt(rA*cA)*EPSABS + EPSREL*max(Y.norm(),X.norm());
		epsdual = sqrt(rA*cA)*EPSABS + EPSREL*rho*U.norm();

#ifdef DEBUG_IT
		cerr << FNAME << " Iteration: " << k;
		cerr << "(rho: " << rho << " remExpense: " << remainingExpensive << ")\t";
		cerr << "Objectives: " << r2 << ' ' << s2 << ' ' << objective << endl;
#endif
		
		k++;

	} while (r2 > epspri and s2 > epsdual and k < maxK and (!expensive or remainingExpensive > 0));

#ifdef DEBUG_CONV
	cerr << FNAME << " Converged after " << k << " iterations ";
	if (k >= maxK) cerr << "by MAXIMUM ITERATIONS";
	else if (expensive && remainingExpensive <= 0) cerr << "by MAXIMUM EXPENSIVE ITERATIONS";
	else if (r2 <= epspri) cerr << "by PRIMAL CRITERION";
	else if (s2 <= epsdual) cerr << "by DUAL CRITERION";
	cerr << endl;
#endif

	// Assign to output
	A_Sigma = Ybest;
}

// A (hopefully) accelerated version of ADMM for WA
// With a budget of "expensive steps" (uses same ratio as for PNFA)
// And backtracking (at the end)
// Possible with over-relaxation
#define OVER_RELAXED
#define OR_ALPHA 1.5
#define FNAME "(accel_admm_wa)"
void accel_admm_wa(MatrixXd & H, MatrixXd & H_Sigma, MatrixXd & A_Sigma, double tau, int maxK, bool Ainitial) {
	
	// For output of debugging
	cerr.precision(8);

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
	if (Ainitial) cerr << "[initial] ";
	cerr << "tau = " << tau;
	cerr << " maxK = " << maxK;
	cerr << endl;
#endif
	
	// Number of terms in the ADMM objective split
	const unsigned int NT = 2;

	// Placeholders for primal and dual variables in ADMM
	MatrixXd X = MatrixXd::Zero(rA, cA);
	MatrixXd Y = MatrixXd::Zero(rA, cA);
	MatrixXd Yprev = MatrixXd::Zero(rA, cA);
	MatrixXd U = MatrixXd::Zero(rA, cA);

	if (Ainitial) {
		Y = A_Sigma;
	}

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
		// Save variables from previous iteration
		Yprev = Y;

		// Update primal variables
		// Proximal of loss
		X = (Aux1 + rho*(Y - U))*((Aux2 + rho*Aux3).inverse());
		// Proximal of trace norm
#ifdef OVER_RELAXED
		Y = OR_ALPHA*X + (1.0-OR_ALPHA)*Yprev + U;
#else
		Y = X + U;
#endif
		nn_projection(Y, tau/rho);

		// Update dual variables
#ifdef OVER_RELAXED
		U = U + (OR_ALPHA*X + (1.0-OR_ALPHA)*Yprev - Y);
#else
		U = U + (X - Y);
#endif

		// Compute residuals (to measure convergence)
		r2 = (X-Y).norm();
		s2 = rho*(Y-Yprev).norm();
		
		// Update rho
		if (r2 > 10*s2) rho = 2*rho;
		else if (s2 > 10*r2) rho = rho/2;

		// Update stopping criterion
		epspri = sqrt(rA*cA)*EPSABS + EPSREL*max(Y.norm(),X.norm());
		epsdual = sqrt(rA*cA)*EPSABS + EPSREL*rho*U.norm();

#ifdef DEBUG_IT
		cerr << FNAME << " Iteration: " << k;
		cerr << "(rho: " << rho << ")\t";
		cerr << "Objectives: " << r2 << ' ' << s2 << ' ' << endl;
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
	A_Sigma = Y;
}
