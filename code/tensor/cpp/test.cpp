#include "include/test.h"

void test_decomposition() {
	cout << "======= TEST_DECOMPOSITION =======" << endl;
	// Dimension
	int d = 100;
	// Number of mixed components
	int m = 5;
	// Number of evecs/evals
	int n = 10;
	// Obtain a random symmetric tensor (as a mixture of several "rank 1" symmetric tensors)
	Lvector v(d);
	Ltensor t(d,d,d);
        Index<'i'> i;
        Index<'j'> j;
        Index<'k'> k;
	for (int it = 0; it < m; it++) {
		random_spherical_vector(v,d);
		t(i,j,k) += v(i)*v(j)*v(k);
	}
	// Obtain the decomposition
	Ltensor R;
	Lmatrix evec;
	Lvector eval;
	double res;
	res = tensor_decomp(t,evec,eval,n,R);
	cout << "Residual: " << res << endl;
}

void test_symmetries() {
	cout << "======= TEST_SYMMETRIES =======" << endl;
	// Dimension
	int d = 50;
	// Number of mixed components
	int c = 5;
	// Number of evecs/evals
	int n = 10;
	// Obtain a random symmetric tensor (as a mixture of several "rank 1" symmetric tensors)
	Lvector v(d);
	Lmatrix m(d,d);
	Ltensor t(d,d,d);
        Index<'i'> i;
        Index<'j'> j;
        Index<'k'> k;
	for (int it = 0; it < c; it++) {
		random_spherical_vector(v,d);
		t(i,j,k) += v(i)*v(j)*v(k);
		m(i,j) += v(i)*v(j);
	}
	// Check whether m and t pass the symmetry test
	bool check;
	MatrixXd em(d,d);
	lmatrix2eigen(m,em);
	check = is_matrix_symmetric(em,true);
	cout << "Matrix symmetric? " << check << endl;
	MatrixXd sm;
	sm = symmetrize_matrix(em);
	check = is_matrix_symmetric(sm,true);
	cout << "Matrix symmetric? " << check << endl;
	check = is_tensor_symmetric(t,true);
	cout << "Tensor symmetric? " << check << endl;
	Ltensor st;
	st = symmetrize_tensor(t);
	check = is_tensor_symmetric(st,true);
	cout << "Tensor symmetric? " << check << endl;
}

void test_eigendecomposition() {
	cout << "======= TEST_EIGENDECOMPOSITION =======" << endl;
	// Dimension
	int d = 400;
	// Number of mixed components
	int c = 10;
	// Obtain a random symmetric tensor (as a mixture of several "rank 1" symmetric tensors)
	Lvector v(d);
	Lmatrix m(d,d);
        Index<'i'> i;
        Index<'j'> j;
	for (int it = 0; it < c; it++) {
		random_spherical_vector(v,d);
		m(i,j) += v(i)*v(j);
	}
	// Convert to Eigen and amplify size of matrix
	MatrixXd em(d,d);
	lmatrix2eigen(m,em);
	cout << "Matrix norm (before amplification): " << em.norm() << endl;
	em *= d;
	cout << "Matrix norm (after amplification): " << em.norm() << endl;
	// Check whether m and t pass the symmetry test
	bool check;
	check = is_matrix_symmetric(em,true);
	cout << "Matrix symmetric? " << check << endl;
	// Run eigendecomposition for symmetric matrices
	SelfAdjointEigenSolver<MatrixXd> saes(em);
	cout << "Eigenvalues: " << saes.eigenvalues().transpose().reverse().head(20) << " ..." << endl;
	// Check if eigenvectors are orthonormalized
	double defect = (MatrixXd::Identity(d,d) - saes.eigenvectors()*saes.eigenvectors().transpose()).norm();
	cout << "||I - V*V^T|| = " << defect << endl;
	// Compute low rank reconstruction
	MatrixXd Vhat(d,c);
	Vhat = saes.eigenvectors().rightCols(c);
	MatrixXd Dhat(c,c);
	Dhat = saes.eigenvalues().tail(c).asDiagonal();
	MatrixXd Mhat(d,d);
	Mhat = Vhat*Dhat*Vhat.transpose();
	cout << "||M - Mhat|| = " << (em-Mhat).norm() << endl;
	// See how far from complete whitening we get
	MatrixXd Wm(d,c);
	Wm = Vhat*Dhat.inverse().cwiseSqrt();
	cout << "||I - W^T*M*W|| = " << (MatrixXd::Identity(c,c) - Wm.transpose()*em*Wm).norm() << endl;
	// Compute the pseudo-inverse of the transpose of Wm
	MatrixXd ptWm(d,c);
	ptWm = Vhat*Dhat.cwiseSqrt();
	cout << "||I - W^T*ptW|| = " << (MatrixXd::Identity(c,c) - Wm.transpose()*ptWm).norm() << endl;
}

void test_learnOtilde() {
	cout << "======= TEST_LEARNOTILDE =======" << endl;
	// Dimension
	int d = 10;
	// Number of mixed components
	int c = 3;
	// Build symmetric matrices and tensors from same vectors
	Lvector v(d);
	Lmatrix m(d,d);
	Ltensor t(d,d,d);
        Index<'i'> i;
        Index<'j'> j;
        Index<'k'> k;
	for (int it = 0; it < c; it++) {
		random_spherical_vector(v,d);
		m(i,j) += v(i)*v(j);
		t(i,j,k) += v(i)*v(j)*v(k);
	}
	// Convert to Eigen
	MatrixXd em(d,d);
	lmatrix2eigen(m,em);
	// Call the learner
	MatrixXd O;
	VectorXd gamma;
	learn_Otilde(em,t,c,O,gamma);
	// Some info and sanity checks
	cout << "Eigenvalues: " << gamma.transpose() << endl;
	// See if we recover the matrix (using the tensor's eigenvectors!)
	MatrixXd M(d,d);
	M = MatrixXd::Zero(d,d);
	for (int it = 0; it < c; it++) {
		M += (O.col(it)*O.col(it).transpose())/(gamma(it)*gamma(it));
	}
	cout << "||M - em|| = " << (em-M).norm() << endl;
}

void test_data_symmetry() {
	cout << "======= TEST_DATA_SYMMETRY =======" << endl;
	// Define model data
	int nsymbols = 4;
	int nstates = 3;
	// Generate initial probabilities
	VectorXd a0(nstates);
	a0 = VectorXd::Random(nstates).cwiseAbs();
	a0 = a0 / a0.sum();
	cout << "a0 : " << a0.transpose() << endl;
	// Generate final probabilities
	VectorXd ainf(nstates);
	ainf = VectorXd::Random(nstates).cwiseAbs();
	ainf = 0.2*ainf / ainf.sum();
	cout << "ainf : " << ainf.transpose() << endl;
	// Generate operators for PNFA
	vector<MatrixXd> pnfa = vector<MatrixXd>();
	MatrixXd A;
	A = MatrixXd::Random(nsymbols*nstates,nstates).cwiseAbs();
	pnfa_projection(A,ainf);
	for (int it = 0; it < nsymbols; it++) {
		MatrixXd B;
		B = A.block(nstates*it,0,nstates,nstates);
		pnfa.push_back(B);
	}
	// Generator operators for HMM
	vector<MatrixXd> hmm = vector<MatrixXd>();
	MatrixXd T, O;
	T = MatrixXd::Random(nstates,nstates).cwiseAbs();
	pnfa_projection(T,VectorXd::Zero(nstates));
	O = MatrixXd::Random(nstates,nsymbols).cwiseAbs();
	for (int it = 0; it < nstates; it++) {
		VectorXd o = O.row(it).transpose();
		simplex_projection(o,1.0 - ainf(it));
		O.row(it) = o.transpose();
	}
	for (int it = 0; it < nsymbols; it++) {
		hmm.push_back(O.col(it).asDiagonal()*T);
	}
	// Check everything is okay
	MatrixXd Mp, Mh;
	Mp = MatrixXd::Zero(nstates,nstates);
	Mh = MatrixXd::Zero(nstates,nstates);
	for (int it = 0; it < nsymbols; it++) {
		Mp += pnfa[it];
		Mh += hmm[it];
	}
	cout << "Sum PNFA operators: ";
	cout << (ainf + Mp*VectorXd::Ones(nstates)).transpose();
	cout << endl;
	cout << "Sum HMM operators: ";
	cout << (ainf + Mh*VectorXd::Ones(nstates)).transpose();
	cout << endl;
	// Apply to some strings
	vector<int> s;
	for (int it = 0; it < 3; it++) s.push_back(it);
	cout << "PNFA : " << evalwa(a0,ainf,pnfa,s) << endl;
	cout << "HMM : " << evalwa(a0,ainf,hmm,s) << endl;
	// Build Hankel matrices and tensors
	MatrixXd pHpa, pHas, hHpa, hHas;
	pHpa = MatrixXd::Zero(nsymbols,nsymbols);
	pHas = MatrixXd::Zero(nsymbols,nsymbols);
	hHpa = MatrixXd::Zero(nsymbols,nsymbols);
	hHas = MatrixXd::Zero(nsymbols,nsymbols);
	// First the ones were integration wrt to operators is easy
	VectorXd pap, pas, hap, has;
	pap = (a0.transpose()*Mp).transpose();
	hap = (a0.transpose()*Mh).transpose();
	pas = Mp*ainf;
	has = Mh*ainf;
	s = vector<int>(2,0);
	for (int it = 0; it < nsymbols; it++) {
		s[0] = it;
		for (int jt = 0; jt < nsymbols; jt++) {
			s[1] = jt;
			pHpa(it,jt) = evalwa(a0,pas,pnfa,s);
			hHpa(it,jt) = evalwa(a0,has,hmm,s);
			pHas(it,jt) = evalwa(pap,ainf,pnfa,s);
			hHas(it,jt) = evalwa(hap,ainf,hmm,s);
		}
	}
	// Now the tensor and the middle integration
	MatrixXd pHps, hHps;
	pHps = MatrixXd::Zero(nsymbols,nsymbols);
	hHps = MatrixXd::Zero(nsymbols,nsymbols);
	Ltensor pHpas(nsymbols,nsymbols,nsymbols), hHpas(nsymbols,nsymbols,nsymbols);
	s = vector<int>(3,0);
	for (int it = 0; it < nsymbols; it++) {
		s[0] = it;
		for (int jt = 0; jt < nsymbols; jt++) {
			s[1] = jt;
			for (int kt = 0; kt < nsymbols; kt++) {
				s[2] = kt;
				double p = evalwa(a0,ainf,pnfa,s);
				double h = evalwa(a0,ainf,hmm,s);
				pHpas(it,jt,kt) = p;
				hHpas(it,jt,kt) = h;
				pHps(it,kt) += p;
				hHps(it,kt) += h;
			}
		}
	}
	// Check sums (must all be equal)
	cout << "Sums PNFA ||";
	cout << " Hpa: " << pHpa.sum();
	cout << " Hps: " << pHps.sum();
	cout << " Has: " << pHas.sum();
	cout << " Hpas: " << sum_ltensor(pHpas);
	cout << endl;
	cout << "Sums HMM ||";
	cout << " Hpa: " << hHpa.sum();
	cout << " Hps: " << hHps.sum();
	cout << " Has: " << hHas.sum();
	cout << " Hpas: " << sum_ltensor(hHpas);
	cout << endl;
	// Do the (alleged) symmetrization of matrix
	MatrixXd pN, hN;
	computeNfromLowRankInverse(pHps,nstates,pN);
	computeNfromLowRankInverse(hHps,nstates,hN);
	MatrixXd pX, hX;
	pX = pHas*pN*pHpa;
	hX = hHas*hN*hHpa;
	// Do the (alleged) symmetrization of tensor
	Ltensor pY(nsymbols,nsymbols,nsymbols), hY(nsymbols,nsymbols,nsymbols);
	Lmatrix pSp(nsymbols,nsymbols), pSs(nsymbols,nsymbols), hSp(nsymbols,nsymbols), hSs(nsymbols,nsymbols);
	eigen2lmatrix((pHas*pN).transpose(),pSp);
	eigen2lmatrix((hHas*hN).transpose(),hSp);
	eigen2lmatrix(pN*pHpa,pSs);
	eigen2lmatrix(hN*hHpa,hSs);
        Index<'i'> i;
        Index<'j'> j;
        Index<'k'> k;
	Index<'a'> a;
	Index<'b'> b;
	pY(i,j,k) = pHpas(a,j,b)*pSp(a,i)*pSs(b,k);
	hY(i,j,k) = hHpas(a,j,b)*hSp(a,i)*hSs(b,k);
	// Check the symmetries
	cout << "(PNFA) matrix X symmetric? ";
	is_matrix_symmetric(pX);
	cout << "(PNFA) tensor Y symmetric? ";
	is_tensor_symmetric(pY);
	cout << "(HMM) matrix X symmetric? ";
	is_matrix_symmetric(hX);
	cout << "(HMM) tensor Y symmetric? ";
	is_tensor_symmetric(hY);
	// Call the learner to see what happens
	MatrixXd pO, hO;
	VectorXd pgamma, hgamma;
	double err;
	err = learn_Otilde(pX,pY,nstates,pO,pgamma);
	cout << "(PNFA) With # of states: " << nstates << endl;
	cout << "(PNFA) Eigenvalues: " << pgamma.transpose() << endl;
	cout << "(PNFA) Residual from tensor decomposition: " << err << endl;
	// This breaks because though pX is not symmetic, it has rank = nstates
	//err = learn_Otilde(pX,pY,nsymbols,pO,pgamma);
	//cout << "(PNFA) With # of states: " << nsymbols << endl;
	//cout << "(PNFA) Eigenvalues: " << pgamma.transpose() << endl;
	//cout << "(PNFA) Residual from tensor decomposition: " << err << endl;
	err = learn_Otilde(hX,hY,nstates,hO,hgamma);
	cout << "(HMM) With # of states: " << nstates << endl;
	cout << "(HMM) Eigenvalues: " << hgamma.transpose() << endl;
	cout << "(HMM) Residual from tensor decomposition: " << err << endl;
}	


void run_tests() {
	test_decomposition();
	test_symmetries();
	test_eigendecomposition();
	test_learnOtilde();
	test_data_symmetry();
}

