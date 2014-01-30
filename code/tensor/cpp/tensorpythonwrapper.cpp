#include "include/tensorpythonwrapper.h"


void learn_tensor(double* Xdata, double* Ydata, int n, double * Odata, double* gammadata, int num_symbols, bool symmetrizeY)
{
	Map<MatrixXd> X_temp(Xdata, num_symbols, num_symbols);
	MatrixXd X(num_symbols,num_symbols);
	X = X_temp;
	Map<MatrixXd> O_temp(Odata, num_symbols, n);
	MatrixXd Otilde(num_symbols, n);
	Otilde = O_temp;
	Map<VectorXd> gamma_temp(gammadata, n);
	VectorXd gamma(n);
	gamma = gamma_temp;
	Ltensor Y(num_symbols, num_symbols, num_symbols);
	Y.fromCArray(Ydata, num_symbols, num_symbols, num_symbols);

	learn_Otilde(X, Y, n, Otilde, gamma, symmetrizeY);
	O_temp = Otilde;
	gamma_temp = gamma;
}

bool test_transfer_tensor(double* tensorData)
{
	cout << "======= TEST TENSOR TRANSFER =======" << endl;
	Ltensor testTensor(2,2,2);
	testTensor.fromCArray(tensorData, 2,2,2);

	int d1 = testTensor.get_dim1();
	int d2 = testTensor.get_dim2();
	int d3 = testTensor.get_dim3();
	if (d1 != d2 or d2 != d3 or d1 != 2) {
		cerr << "Wrong dimensions. Wanted 2x2x2 and got: " << d1 << "x" << d2 << "x" << d3 << endl;
		return false ;
	}

	int i = 0, j=0, k=0, counter=0;
	bool properTransfer = true;
	for(i=0;i<2;i++)
	{
		for(j=0;j<2;j++)
		{
			for(k=0;k<2;k++)
			{
				cout << "i: " << i << " j: " << j << " k: " << k << "  Value: " << testTensor(i,j,k) << endl;
				if(testTensor(i,j,k) != counter)
				{
						properTransfer = false;
				}
				counter++;
			}
		}
	}

	cout << "Proper transfer? " << properTransfer << endl;
	return properTransfer;
}

bool test_Y_symmetry(double *Ydata, int num_symbols)
{
	cout << "======= TEST Y Symmetry =======" << endl;
	Ltensor Y(num_symbols, num_symbols, num_symbols);
	Y.fromCArray(Ydata, num_symbols, num_symbols, num_symbols);	
	bool symcheck = is_tensor_symmetric(Y,true);
	return symcheck;
}

void do_simplex_projection(double* vectorData, double norm, int length)
{
	Map<VectorXd> vec_temp(vectorData, length);
	VectorXd vec(length);
	vec = vec_temp;

	simplex_projection(vec, norm);
	vec_temp = vec;
}

void do_probmat_projection(double* matrixData, double* normData, int m, int n)
{
	Map<MatrixXd> mat_temp(matrixData, m, n);
	MatrixXd mat(m,n);
	mat = mat_temp;
	Map<VectorXd> vec_temp(normData, m);
	VectorXd vec(m);
	vec = vec_temp;

	probmat_projection(mat, vec);
	mat_temp = mat;
}




void get_real_HMM_hankels(double* Hpsdata, double* Hpsigdata, double* Hsigsdata, double* tensorData, double *HpsnontensorData, double* HppsData)
{
	// Define model data
	int nsymbols = 4;
	int nstates = 3;

	Map<MatrixXd> H_ps_temp(Hpsdata, nsymbols, nsymbols);
	MatrixXd hHps(nsymbols,nsymbols);
	hHps = H_ps_temp;

	Map<MatrixXd> H_psig_temp(Hpsigdata, nsymbols, nsymbols);
	MatrixXd hHpa(nsymbols,nsymbols);
	hHpa = H_psig_temp;
	
	Map<MatrixXd> H_sigs_temp(Hsigsdata, nsymbols, nsymbols);
	MatrixXd hHas(nsymbols,nsymbols);
	hHas = H_sigs_temp;

	Map<MatrixXd> H_psnontense_temp(HpsnontensorData, nsymbols, nsymbols);
	MatrixXd Hps(nsymbols,nsymbols);
	Hps = H_psnontense_temp;

	Map<MatrixXd> H_pps_temp(HppsData, nsymbols, nsymbols);
	MatrixXd Hpps(nsymbols,nsymbols);
	Hpps = H_pps_temp;

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
	MatrixXd Mh;
	Mh = MatrixXd::Zero(nstates,nstates);
	for (int it = 0; it < nsymbols; it++) {
		Mh += hmm[it];
	}

	cout << "Sum HMM operators: ";
	cout << (ainf + Mh*VectorXd::Ones(nstates)).transpose();
	cout << endl;
	// Apply to some strings
	vector<int> s;
	for (int it = 0; it < 3; it++) s.push_back(it);
	cout << "HMM : " << evalwa(a0,ainf,hmm,s) << endl;
	
	// First the ones were integration wrt to operators is easy
	VectorXd hap, has;
	hap = (a0.transpose()*Mh).transpose();
	has = Mh*ainf;
	s = vector<int>(2,0);
	for (int it = 0; it < nsymbols; it++) {
		s[0] = it;
		for (int jt = 0; jt < nsymbols; jt++) {
			s[1] = jt;
			hHpa(it,jt) = evalwa(a0,has,hmm,s);
			hHas(it,jt) = evalwa(hap,ainf,hmm,s);
		}
	}

	VectorXd prefainf(nstates);
	prefainf = VectorXd::Ones(nstates);
	hap = (a0.transpose()*Mh).transpose();
	s = vector<int>(2,0);
	for (int it = 0; it < nsymbols; it++) {
		s[0] = it;
		for (int jt = 0; jt < nsymbols; jt++) {
			s[1] = jt;
			Hps(it,jt) = evalwa(a0,ainf,hmm,s);
			Hpps(it,jt) = evalwa(a0,prefainf,hmm,s);
		}
	}

	// Now the tensor and the middle integration
	Ltensor hHpas(nsymbols,nsymbols,nsymbols);
	s = vector<int>(3,0);
	for (int it = 0; it < nsymbols; it++) {
		s[0] = it;
		for (int jt = 0; jt < nsymbols; jt++) {
			s[1] = jt;
			for (int kt = 0; kt < nsymbols; kt++) {
				s[2] = kt;
				double h = evalwa(a0,ainf,hmm,s);
				hHpas(it,jt,kt) = h;
				hHps(it,kt) += h;
			}
		}
	}

	cout << "Sums HMM ||";
	cout << " Hpa: " << hHpa.sum();
	cout << " Hps: " << hHps.sum();
	cout << " Has: " << hHas.sum();
	cout << " Hpas: " << sum_ltensor(hHpas);
	cout << endl;
	cout << "Hpa in C++: " << hHpa << endl;
	cout << "Hps in C++: " << hHps << endl;
	cout << "Has in C++: " << hHas << endl;

	H_ps_temp = hHps;
	H_psig_temp = hHpa;
	H_sigs_temp = hHas;
	H_pps_temp = Hpps;
	H_psnontense_temp = Hps;
	hHpas.toCArray(tensorData);

		// Do the (alleged) symmetrization of matrix
	MatrixXd hN;
	computeNfromLowRankInverse(hHps,nstates,hN);
	MatrixXd hX;
	hX = hHas*hN*hHpa;
	cout << "X in CPP: " << hX << endl;
	// Do the (alleged) symmetrization of tensor
	Ltensor hY(nsymbols,nsymbols,nsymbols);
	Lmatrix  hSp(nsymbols,nsymbols), hSs(nsymbols,nsymbols);
	eigen2lmatrix((hHas*hN).transpose(),hSp);
	eigen2lmatrix(hN*hHpa,hSs);
        Index<'i'> i;
        Index<'j'> j;
        Index<'k'> k;
	Index<'a'> a;
	Index<'b'> b;
	hY(i,j,k) = hHpas(a,j,b)*hSp(a,i)*hSs(b,k);
	cout << "(HMM) matrix X symmetric? ";
	is_matrix_symmetric(hX);
	cout << "(HMM) tensor Y symmetric? ";
	is_tensor_symmetric(hY);
}