//==========================================================
// Author : Lixin He, Mohan Chen
// Update : Peize Lin
// Last Update : 2018-09-04
//==========================================================
#ifndef COMPLEXMATRIX_H
#define COMPLEXMATRIX_H

#include <complex>
using namespace std;

#include "matrix.h"
#include "utils_io.h"

//#ifdef _MCD_CHECK
//#include "src_parallel/mcd.h"
//#endif

class ComplexMatrix
{

public:

	int nr=0;
	int nc=0;
	int size=0;
	complex<double> *c=nullptr;

	ComplexMatrix(): nr(0), nc(0), size(0), c(nullptr){}
	ComplexMatrix(const int nrows,const int ncols,const bool flag_zero=true);		// Peize Lin add flag_zero 2019-05-13
	ComplexMatrix(const ComplexMatrix &m1);
	ComplexMatrix(ComplexMatrix && m1);						// Peize Lin add 2016-08-05
	explicit ComplexMatrix(const matrix &m);							// Peize Lin add 2017-03-29
	~ComplexMatrix();
	
	void create(const int nrow,const int ncol,const bool flag_zero=true);		// Peize Lin add flag_zero 2019-05-13
	ComplexMatrix& operator=(const ComplexMatrix &m);
	ComplexMatrix& operator=(ComplexMatrix && m);			// Peize Lin add 2016-08-05

	//============
	// Operators
	//============
	complex<double> &operator()(const int ir,const int ic)
	{
        if (ir >= nr) LIBRPA::utils::lib_printf("ir %d nr %d\n", ir, nr);
        if (ic >= nc) LIBRPA::utils::lib_printf("ic %d nc %d\n", ic, nc);
		assert(ir>=0);	assert(ir<nr);	assert(ic>=0);	assert(ic<nc);
		return c[ir*nc+ic];//mohan modify in-line 2007-10-1
	}
	const complex<double> &operator()(const int ir,const int ic)const
	{
        if (ir >= nr) LIBRPA::utils::lib_printf("ir %d nr %d\n", ir, nr);
        if (ic >= nc) LIBRPA::utils::lib_printf("ic %d nc %d\n", ic, nc);
		assert(ir>=0);	assert(ir<nr);	assert(ic>=0);	assert(ic<nc);
		return c[ir*nc+ic];//mohan modify in-line 2007-10-13
	}
	
	ComplexMatrix& operator*=(const complex<double> &s);
	ComplexMatrix& operator+=(const ComplexMatrix &m);
	ComplexMatrix& operator+=(const complex<double> &s); // minez add 2022-06-06
	ComplexMatrix& operator-=(const ComplexMatrix &m);
	matrix real() const;						// Peize Lin add 2017-03-29
	matrix imag() const;
	
	//==================
	// member function:
	//==================
	void zero_out(void);
	void set_as_identity_matrix(void);
    //! Get the max real value of matrix element
    double get_max_real() const; // minyez add 2022-05-05
    double get_max_real(int &ir, int &ic) const; // minyez add 2022-05-05
    //! Get the max imaginary value of matrix element
    double get_max_imag() const; // minyez add 2022-10-18
    double get_max_abs() const; // minyez add 2024-03-13
    double get_max_abs_imag() const; // minyez add 2022-10-26
    double get_max_abs_offdiag() const;
    bool is_diagonal(const double &thres = 1e-14) const;
};

ComplexMatrix operator+(const ComplexMatrix &m1,  const ComplexMatrix &m2);
ComplexMatrix operator-(const ComplexMatrix &m1,  const ComplexMatrix &m2);
ComplexMatrix operator-(const complex<double> &s, const ComplexMatrix &m); // minyez add 2022-06-06
ComplexMatrix operator-(const ComplexMatrix &m,   const complex<double> &s); // minyez add 2022-06-06
ComplexMatrix operator*(const ComplexMatrix &m1,  const ComplexMatrix &m2);
ComplexMatrix operator*(const complex<double> &s, const ComplexMatrix &m);
ComplexMatrix operator*(const ComplexMatrix &m,   const complex<double> &s);
ComplexMatrix operator*(const double &s,          const ComplexMatrix &m);
ComplexMatrix operator*(const ComplexMatrix &m,   const double &s);
	
complex<double> trace(const ComplexMatrix &m);

double abs2_row(const ComplexMatrix &m,const int ir);		// mohan add 2008-7-1
double abs2_column(const ComplexMatrix &m,const int ic);	// mohan add 2008-7-1
double abs2(const ComplexMatrix &m);
double abs2(const int nmat, ComplexMatrix **m);

ComplexMatrix transpose(const ComplexMatrix &m, const bool &conjugate);
ComplexMatrix conj(const ComplexMatrix &m);						// Peize Lin add 2019-05-13

void scale_accumulate(
		const complex<double> &s, 
		const ComplexMatrix &min, 
		ComplexMatrix &mout);

void scale_accumulate(
		const int &nmat, 
		const complex<double> &s, 
		ComplexMatrix **min, 
		ComplexMatrix **mout);

void scaled_sum(
		const complex<double> &s1, 
		const ComplexMatrix &m1, 
		const complex<double> &s2, 
		const ComplexMatrix &m2, 
		ComplexMatrix &mout);

void scaled_sum(
		const int &nmat,
		const complex<double> &s1, 
		ComplexMatrix **m1, 
		const complex<double> &s2,
		ComplexMatrix **m2, 
		ComplexMatrix **mout);

//! compute the power of Hermitian matrix. Note that the hermicity not checked itself
/*
 * @param cmat: the complex matrix to compute power
 * @param power: the power to compute
 * @param keep_ev: whether to keep the eigenvectors in cmat
 * @param filter_original: whether to filter the small values when recovering the original matrix
 * @param threshold: eigenvalue threshold to filter out in computing the power matrix
 */
ComplexMatrix power_hemat(ComplexMatrix &cmat, double power, bool keep_ev = false, bool filter_original = false,
                          double threshold = -1e16); // Minye Zhang add 2022-06-04

//! Does the same as power_hemat, but save the powered matrix in the original matrix
/*
 * @param cmat: the complex matrix to compute power
 * @param power: the power to compute
 * @param threshold: eigenvalue threshold to filter out in computing the power matrix
 */
void power_hemat_onsite(ComplexMatrix &cmat, double power, double threshold = -1e16); // Minye Zhang add 2022-07-07

void print_complex_matrix(const char *desc, const ComplexMatrix &mat);

void print_complex_matrix_file(const char *desc, const ComplexMatrix &mat, ofstream &fs, bool use_scientific);
void print_complex_matrix_mm(const ComplexMatrix &mat, ofstream &fs, double threshold = 1e-15, bool row_first = true);
void print_complex_matrix_file(const char *desc, const ComplexMatrix &mat, const string &fn, bool use_scientific);
void print_complex_matrix_mm(const ComplexMatrix &mat, const string &fn, double threshold = 1e-15, bool row_first = true);
void print_complex_real_matrix(const char* desc, const ComplexMatrix &mat );
#endif
