/*!
 @file matrix.h
 @brief utilies to handle square matrix and related operations
 @note mostly from Peize Lin's code, with some modification by minyez.
 @date 2018-07-29 update by Peize Lin
 @date 2022-05-06 update by minyez
 */
#ifndef MATRIX_H
#define MATRIX_H

// Peize Lin update 2018-07-29
/*
#ifdef _MCD_CHECK
#include "mcd.h"
#endif
*/
#include<ostream>
#include<cassert>

#include<fstream>// test

class matrix
{
	/* data */
public:
    //! Threshold to treat the matrix as equal, added by minyez 2022-05-06
    static const double EQUAL_THRES;
	int nr=0;
	int nc=0;   /* Number of rows and columns */
    int size=0;
	double *c=nullptr;    /* Holds the data */

	/* Constructors and destructor */
	matrix(): nr(0), nc(0), size(0), c(nullptr){}
	matrix( const int nrows, const int ncols, const bool flag_zero=true );		// Peize Lin add flag_zero 2018-07-02
	matrix( const matrix &m1 ); /* copy constructor */
	matrix( matrix && m1 );			// Peize Lin add 2016-08-05
	~matrix();

	void create( const int nrow, const int ncol, const bool flag_zero=true );			// Peize Lin add flag_zero 2018-07-02
	matrix& operator=(const matrix &m1); // Peize Lin change 2018-03-12
	matrix& operator=( matrix && m1 );	// Peize Lin add 2016-08-05

	double &operator()(const int ir,const int ic)
	{
//		assert(ir>=0);	assert(ir<nr);	assert(ic>=0);	assert(ic<nc);
		return c[ir*nc+ic];
	}

	const double &operator()(const int ir,const int ic) const
	{
//		assert(ir>=0);	assert(ir<nr);	assert(ic>=0);	assert(ic<nc);
		return c[ir*nc+ic];
	}

//	inline double &operator()(int i,int j) const
//	    { return c[nc*i+j]; }

	void operator*=(const double &s);
	void operator+=(const double &s);
	void operator+=(const matrix &m);
	void operator-=(const double &s);
	void operator-=(const matrix &m);
	bool operator==(const matrix &m);

    matrix operator+(double s) const;
    matrix operator-(double s) const;

	/* member function: */
	matrix Inverse(void);

	double det(void);

	// mohan add 2011-01-13
	double trace_on(void) const;

	/* zero out all the entries */
	void zero_out(void);

	void get_extreme_eigen_values(double &ev_lower, double &ev_upper) const;	// mohan add 2011-01-13

	void reshape( const double nr_new, const double nc_new );		// Peize Lin add 2017-05-27

	double max() const;					// Peize Lin add 2016-09-08
	double max(int &ir, int &ic) const; // minyez add 2022-05-05
	double min() const;					// Peize Lin add 2016-09-08
	double absmax() const;				// Peize Lin add 2018-07-02
    //! count the matrix element whose abs value is less equal than threshold
    unsigned count_absle(double thres) const; // minyez add 2022-05-05
    //! count the matrix element whose abs value is greater equal than threshold
    unsigned count_absge(double thres) const; // minyez add 2022-05-05
	//double norm() const;				// Peize Lin add 2018-08-12
};


matrix operator+(const matrix &m1, const matrix &m2);
matrix operator-(const matrix &m1, const matrix &m2);
matrix operator*(const matrix &m1, const matrix &m2);
matrix operator*(const double &s, const matrix &m);
matrix operator*(const matrix &m, const double &s);

matrix transpose(const matrix &m);
double trace_on(const matrix &A, const matrix &B);		// mohan add 2011-01-13
double mdot(const matrix &A, const matrix &B);			// mohan add 2011-01-13

std::ostream & operator<<( std::ostream & os, const matrix & m );		// Peize Lin add 2016-09-08

//! compute the power of a symmetric matrix. Symmetry not checked itself
matrix power_symat(matrix &mat, double power,
                   double threshold = -1e16); // Minye Zhang add 2022-06-04

// minyez copied from Shi Rong's cal_periodic_chi0.cpp, for test use
void print_matrix(const char *desc, const matrix &mat);

void print_matrix_mm(const matrix &mat, std::ostream &os, double threshold = 1e-15, bool row_first = true);
void print_matrix_mm(const matrix &mat, const std::string &fn, double threshold = 1e-15, bool row_first = true);
#endif // MATRIX_H
