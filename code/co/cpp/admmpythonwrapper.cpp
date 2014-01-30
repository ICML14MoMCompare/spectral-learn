#include "include/admmpythonwrapper.h"


void admm_pnfa_learn(double* H, double* H_Sigma, double* ainf, double* A_Sigma, int basisSize, int numSymbols, double tau, int maxK, int maxRank)
{
	Map<MatrixXd> H_2(H, basisSize, basisSize);
	Map<MatrixXd> H_Sigma_2(H_Sigma, basisSize*numSymbols,basisSize);
	Map<VectorXd> ainf_2(ainf,basisSize);
	Map<MatrixXd> A_Sigma_2(A_Sigma,basisSize*numSymbols,basisSize);

	MatrixXd H_Mat(basisSize,basisSize);
	H_Mat = H_2;
	MatrixXd H_Sigma_Mat(basisSize*numSymbols,basisSize);
	H_Sigma_Mat = H_Sigma_2;
	MatrixXd A_Sigma_Mat(basisSize*numSymbols,basisSize);
	A_Sigma_Mat = A_Sigma_2;
	VectorXd ainf_Vec(basisSize);
    ainf_Vec = ainf_2;

   

	accel_trunc_admm(H_Mat, H_Sigma_Mat, ainf_Vec, A_Sigma_Mat, tau, maxK, maxRank);
	accel_admm(H_Mat,H_Sigma_Mat,ainf_Vec,A_Sigma_Mat,tau,4,true);

	A_Sigma_2 = A_Sigma_Mat;
	ainf_2 = ainf_Vec;
}

void admm_wfa_learn(double* H, double* H_Sigma, double* A_Sigma, int basisSize, int numSymbols, double tau, int maxK, int maxRank)
{
	Map<MatrixXd> H_2(H, basisSize, basisSize);
	Map<MatrixXd> H_Sigma_2(H_Sigma, basisSize*numSymbols,basisSize);
	Map<MatrixXd> A_Sigma_2(A_Sigma,basisSize*numSymbols,basisSize);

	MatrixXd H_Mat(basisSize,basisSize);
	H_Mat = H_2;
	MatrixXd H_Sigma_Mat(basisSize*numSymbols,basisSize);
	H_Sigma_Mat = H_Sigma_2;
	MatrixXd A_Sigma_Mat(basisSize*numSymbols,basisSize);
	A_Sigma_Mat = A_Sigma_2;

	accel_trunc_admm_wa(H_Mat, H_Sigma_Mat, A_Sigma_Mat, tau, maxK, maxRank);
	accel_admm_wa(H_Mat,H_Sigma_Mat,A_Sigma_Mat,tau,4,true);

	A_Sigma_2 = A_Sigma_Mat;
}