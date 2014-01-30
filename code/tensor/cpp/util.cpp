#include "include/util.h"

boost::random::mt19937 rng;

// Generate a random spherical vector
void random_spherical_vector(Lvector & v, int d) {
	boost::random::normal_distribution<> std_normal;
	v = Lvector(d);
	for (int i = 0; i < d; i++) v(i) = std_normal(rng);
	v *= (1.0 / v.getNormN(2));
}

/*
// Convert eigen -> lvector
void eigen2lvector(VectorXd & e, Lvector & l) {

}

*/

// Convert eigen -> lmatrix
#define FNAME "(eigen2lmatrix)"
void eigen2lmatrix(const MatrixXd & e, Lmatrix & l) {
	int rl = l.get_dim1();
	int cl = l.get_dim2();
	int re = e.rows();
	int ce = e.cols();
	if (rl != re or cl != ce) {
		cerr << FNAME << " Dimensions mismatch!" << endl;
		throw;
	}
	for (int i = 0; i < re; i++) {
		for (int j = 0; j < ce; j++) {
			l(i,j) = e(i,j);
		}
	}
}

// Convert lmatrix -> eigen
#define FNAME "(lmatrix2eigen)"
void lmatrix2eigen(Lmatrix & l, MatrixXd & e) {
	int rl = l.get_dim1();
	int cl = l.get_dim2();
	int re = e.rows();
	int ce = e.cols();
	if (rl != re or cl != ce) {
		cerr << FNAME << " Dimensions mismatch!\t";
		cerr << "l (" << rl << "*" << cl << ") ";      
		cerr << "e (" << re << "*" << ce << ")";      
		cerr << endl;
		throw;
	}
	for (int i = 0; i < rl; i++) {
		for (int j = 0; j < cl; j++) {
			e(i,j) = l(i,j);
		}
	}
}

// Convert lvector -> eigen
#define FNAME "(lvector2eigen)"
void lvector2eigen(Lvector & l, VectorXd & e) {
	int dl = l.get_dim1();
	int de = e.rows();
	// Check dimensions
	if (dl != de) {
		cerr << FNAME << " Dimensions mismatch!" << endl;
		throw;
	}
	for (int i = 0; i < dl; i++) e(i) = l(i);
}

// Flatten ltensor -> eigen matrix
void ltensor2eigenmatrix(Ltensor & t, MatrixXd & m) {
	int d1 = t.get_dim1();
	int d2 = t.get_dim2();
	int d3 = t.get_dim3();
	int r = d1*d2;
	int c = d3;
	m = MatrixXd(r,c);
	t.toCArray(m.data());
}

// Compute norm of ltensor
double frob_norm_ltensor(Ltensor & t) {
	/*
	MatrixXd m;
	ltensor2eigenmatrix(t, m);
	return m.norm();
	*/
        Index<'i'> i;
        Index<'j'> j;
        Index<'k'> k;
	return sqrt(t(i,j,k)*t(i,j,k));
}

// Compute norm of ltensor
double sum_ltensor(Ltensor & t) {
        Index<'i'> i;
        Index<'j'> j;
        Index<'k'> k;
	Ltensor ones(t.get_dim1(),t.get_dim2(),t.get_dim3(),1.0);
	return t(i,j,k)*ones(i,j,k);
}

#define FNAME "(symmetrize_matrix)"
MatrixXd symmetrize_matrix(MatrixXd & m) {
	if (m.cols() != m.rows()) {
		cerr << FNAME << " Dimensions mismatch" << endl;
		throw;
	}
	MatrixXd sm;
	sm = (m + m.transpose()) / 2.0;
	return sm;
}

// Check symmetry of a matrix
#define MAT_SYM_THR 1E-6
#define FNAME "(is_matrix_symmetric)"
bool is_matrix_symmetric(MatrixXd & m, bool debug) {
	double defect = (m - m.transpose()).norm();
	if (debug) cerr << FNAME << " Defect: " << defect << endl;
	return (defect < MAT_SYM_THR);
}

#define FNAME "(symmetrize_tensor)"
Ltensor symmetrize_tensor(Ltensor & t) {
	int d1 = t.get_dim1();
	int d2 = t.get_dim2();
	int d3 = t.get_dim3();
	if (d1 != d2 or d2 != d3) {
		cerr << FNAME << " Dimensions mismatch" << endl;
		throw;
	}
	Ltensor p(d1,d2,d3);
        Index<'i'> i;
        Index<'j'> j;
        Index<'k'> k;
	p = t;
	p(i,j,k) += t(j,i,k);
	p(i,j,k) += t(j,k,i);
	p(i,j,k) += t(k,j,i);
	p(i,j,k) += t(k,i,j);
	p(i,j,k) += t(i,k,j);
	p *= 1.0 / 6.0;
	return p;
}

// Check symmetry of a tensor
#define TEN_SYM_THR 1E-6
#define FNAME "(is_tensor_symmetric)"
bool is_tensor_symmetric(Ltensor & t, bool debug) {
	int d1 = t.get_dim1();
	int d2 = t.get_dim2();
	int d3 = t.get_dim3();
	if (d1 != d2 or d2 != d3) {
		if (debug) cerr << FNAME << " Dimensions mismatch" << endl;
		return false;
	}
	Ltensor p(d1,d2,d3);
        Index<'i'> i;
        Index<'j'> j;
        Index<'k'> k;
	double defect;
	/*
	defect = 0.0;
	p(i,j,k) = t(i,j,k) - t(j,i,k);
	defect += frob_norm_ltensor(p);
	p(i,j,k) = t(i,j,k) - t(j,k,i);
	defect += frob_norm_ltensor(p);
	p(i,j,k) = t(i,j,k) - t(k,j,i);
	defect += frob_norm_ltensor(p);
	p(i,j,k) = t(i,j,k) - t(k,i,j);
	defect += frob_norm_ltensor(p);
	p(i,j,k) = t(i,j,k) - t(i,k,j);
	defect += frob_norm_ltensor(p);
	defect = defect / 5;
	*/
	p = t;
	p(i,j,k) += t(j,i,k);
	p(i,j,k) += t(j,k,i);
	p(i,j,k) += t(k,j,i);
	p(i,j,k) += t(k,i,j);
	p(i,j,k) += t(i,k,j);
	p *= 1.0 / 6.0;
	p -= t;
	defect = frob_norm_ltensor(p);

	if (debug) cerr << FNAME << " Defect: " << defect << endl;
	return (defect < TEN_SYM_THR);
}

double evalwa(const VectorXd & a0, const VectorXd & ainf, const vector<MatrixXd> & wa, const vector<int> s) {
	RowVectorXd a = a0.transpose();
	for (vector<int>::const_iterator it = s.begin(); it != s.end(); it++) {
		a *= wa[*it];
	}
	return a*ainf;
}

