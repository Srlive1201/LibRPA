/*
 * @file matrix_m.h
 * @brief implementing a matrix template class with support of column-major or row-major style
 */
#pragma once
#include <cassert>
#include <vector>
#include <cmath>
#include <ctime>
#include <random>
#include <valarray>
#include <memory>
#include <functional>
#include "base_utility.h"
#include "scalapack_connector.h"
#include "lapack_connector.h"
#include "vec.h"

// major dependent flatted index of 2D matrix
// row-major
inline int flatid_rm(const int &nr, const int &ir, const int &nc, const int &ic) { return ir*nc+ic; }
// column-major
inline int flatid_cm(const int &nr, const int &ir, const int &nc, const int &ic) { return ic*nr+ir; }

enum MAJOR { ROW, COL };

template <typename T>
void get_nr_nc_from_nested_vector(const std::vector<std::vector<T>> &nested_vec, int &nr, int &nc)
{
    nr = nested_vec.size();
    int nc_ = 0;
    for (const auto& v: nested_vec)
        nc_ = v.size() > nc_? v.size() : nc_;
    nc = nc_;
}

template <typename T>
void expand_nested_vector_to_pointer(const std::vector<std::vector<T>> &nested_vector, const int &nr, const int &nc, T* c, bool row_major = true)
{
    int nc_, nr_;
    get_nr_nc_from_nested_vector(nested_vector, nr_, nc_);
    if (nr < nr_) nr_ = nr;
    if (nc < nc_) nc_ = nc;

    if (nr&&nc)
    {
        if (row_major)
        {
            for (int ir = 0; ir < nr_; ir++)
                for (int ic = 0; ic < std::min(size_t(nc_), nested_vector[ir].size()); ic++)
                    c[ir*nc+ic] = nested_vector[ir][ic];
        }
        else
        {
            for (int ir = 0; ir < nr_; ir++)
                for (int ic = 0; ic < std::min(size_t(nc_), nested_vector[ir].size()); ic++)
                    c[ic*nr+ir] = nested_vector[ir][ic];
        }
    }
}

template <typename T>
class matrix_m
{
private:
    int mrank_;
    int size_;
    int ld_;
    int nr_;
    int nc_;
    MAJOR major_;
    std::function<int(const int& nr, const int& ir, const int& nc, const int& ic)> indx_picker_2d_;

public:
    std::shared_ptr<std::valarray<T>> data; // not used yet
    T* c;

private:
    void assign_value(const T &v)
    {
        for (int i = 0; i < size_; i++)
            c[i] = v;
    }
    void assign_value(const T* const pv)
    {
        // do not check out of bound
        // data are copied from pv as it is, despite the major
        for (int i = 0; i < size_; i++)
            c[i] = pv[i];
    }

    void allocate_storage()
    {
        if (nr_ && nc_)
        {
            c = new T[nr_ * nc_];
        }
    }
    void deallocate_storage()
    {
        if (c)
        {
            delete [] c;
            c = nullptr;
        }
    }
    void reallocate_storage(int nrows_new, int ncols_new)
    {
        const int size_new = nrows_new * ncols_new;
        if (size_new)
        {
            if (c)
            {
                if ( size_new != size() )
                {
                    delete [] c;
                    c = new T[size_new];
                }
            }
            else
                c = new T[size_new];
        }
        else
        {
            if(c) delete [] c;
            c = nullptr;
        }
    }

    void set_rank_size()
    {
        mrank_ = std::min(nr_, nc_);
        size_ = nr_ * nc_;
    }

    void set_indx_picker()
    {
        if (major_ == MAJOR::ROW)
        {
            ld_ = nc_;
            indx_picker_2d_ = flatid_rm;
        }
        else
        {
            ld_ = nr_;
            indx_picker_2d_ = flatid_cm;
        }
    }

public:
    using type = T;
    using real_t = typename to_real<T>::type;
    using cplx_t = typename to_cplx<T>::type;
    const bool is_complex = is_complex_t<T>::value;

    // constructors
    matrix_m() : nr_(0), nc_(0), major_(MAJOR::ROW), c(nullptr)
    {
        set_rank_size();
        set_indx_picker();
    }
    matrix_m(const int &nrows, const int &ncols, MAJOR major = MAJOR::ROW): nr_(nrows), nc_(ncols), major_(major), c(nullptr)
    {
        set_rank_size();
        set_indx_picker();
        allocate_storage();
        zero_out();
    }
    matrix_m(const int &nrows, const int &ncols, const T * const valarr, MAJOR major = MAJOR::ROW): nr_(nrows), nc_(ncols), major_(major), c(nullptr)
    {
        set_rank_size();
        set_indx_picker();
        allocate_storage();
        zero_out();
        assign_value(valarr);
    }
    // constructor from a nested vector
    matrix_m(const std::vector<std::vector<T>> &nested_vector, MAJOR major = MAJOR::ROW): major_(major), c(nullptr)
    {
        get_nr_nc_from_nested_vector(nested_vector, nr_, nc_);
        set_rank_size();
        set_indx_picker();
        allocate_storage();
        expand_nested_vector_to_pointer(nested_vector, nr_, nc_, c, is_row_major());
    }

    // copy constructor
    matrix_m(const matrix_m<T> &m): nr_(m.nr_), nc_(m.nc_), major_(m.major_), c(nullptr)
    {
        set_rank_size();
        set_indx_picker();
        allocate_storage();
        memcpy(c, m.c, nr_*nc_*sizeof(T));
    }

    matrix_m(matrix_m &&m) : nr_(m.nr_), nc_(m.nc_), major_(m.major_), c(nullptr)
    {
        c = m.c;
        mrank_ = m.mrank_;
        size_ = m.size_;
        indx_picker_2d_ = m.indx_picker_2d_;
        m.nr_ = m.nc_ = 0;
        m.c = nullptr;
    }
    // destructor
    ~matrix_m()
    {
        nc_ = nr_ = mrank_ = size_ = 0;
        deallocate_storage();
    }

    void zero_out()
    {
        assign_value(T(0));
    }

    void randomize(const T &lb = 0, const T &ub = 1)
    {
        std::default_random_engine e(time(0));
        std::uniform_real_distribution<real_t> dr(::get_real(lb), ::get_real(ub)), di(::get_imag(lb), ::get_imag(ub));
        if (is_complex)
        {
            for (int i = 0; i != size_; i++)
                this->c[i] = std::complex<T>{dr(e), di(e)};
        }
        else
        {
            for (int i = 0; i != size_; i++)
                this->c[i] = dr(e);
        }
    }

    // access to private variables
    int nr() const { return nr_; }
    int nc() const { return nc_; }
    // leading dimension
    int ld() const { return ld_; }
    MAJOR major() const { return major_; }
    int size() const { return size_; }

    bool is_row_major() const { return major_ == MAJOR::ROW; }
    bool is_col_major() const { return major_ == MAJOR::COL; }

    // indexing
    T &operator()(const int ir, const int ic) { return c[indx_picker_2d_(nr_, ir, nc_, ic)]; }
    const T &operator()(const int ir, const int ic) const { return this->at(ir, ic); }
    const T &at(const int ir, const int ic) const { return c[indx_picker_2d_(nr_, ir, nc_, ic)]; }

    // swap major
    void swap_to_row_major()
    {
        if (is_col_major())
        {
            int nr = nr_, nc = nc_;
            this->transpose();
            reshape(nr, nc);
            major_ = MAJOR::ROW;
            set_indx_picker();
        }
    };
    void swap_to_col_major()
    {
        if (is_row_major())
        {
            int nr = nr_, nc = nc_;
            this->transpose();
            reshape(nr, nc);
            major_ = MAJOR::COL;
            set_indx_picker();
        }
    };

    // assignment copy
    matrix_m<T> & operator=(const matrix_m<T> &m)
    {
        if (this == &m) return *this;
        resize(m.nr_, m.nc_);
        major_ = m.major_;
        set_indx_picker();
        memcpy(c, m.c, nr_*nc_*sizeof(T));
        return *this;
    }

    matrix_m<T> & operator=(matrix_m<T> &&m)
    {
        if (this == &m) return *this;
        nr_ = m.nr_;
        nc_ = m.nc_;
        mrank_ = m.mrank_;
        size_ = m.size_;
        major_ = m.major_;
        set_indx_picker();
        if(c) delete [] c;
        c = m.c;
        m.nr_ = m.nc_ = 0;
        m.c = nullptr;
        return *this;
    }

    // operator overload
    matrix_m<T> & operator=(const T &cnum)
    {
        for (int i = 0; i < size(); i++)
            c[i] = cnum;
        return *this;
    }

    template <typename T1>
    void operator+=(const matrix_m<T1> &m)
    {
        assert(size() == m.size() && major() == m.major());
        for (int i = 0; i < size(); i++)
            c[i] += m.c[i];
    }
    template <typename T1>
    void operator+=(const std::vector<T1> &v)
    {
        assert(nc_ == v.size());
        for (int i = 0; i < nr_; i++)
            for (int ic = 0; ic < nc_; ic++)
                this->at(i, ic) += v[ic];
    }
    template <typename T1>
    void operator+=(const T1 &cnum)
    {
        for (int i = 0; i < size(); i++)
            c[i] += cnum;
    }

    template <typename T1>
    void operator-=(const matrix_m<T1> &m)
    {
        assert(size() == m.size() && major() == m.major());
        for (int i = 0; i < size(); i++)
            c[i] -= m.c[i];
    }
    template <typename T1>
    void operator-=(const std::vector<T1> &v)
    {
        assert(nc_ == v.size());
        for (int i = 0; i < nr_; i++)
            for (int ic = 1; ic < nc_; ic++)
                this->at(i, ic) += v[ic];
    }
    template <typename T1>
    void operator-=(const T1 &cnum)
    {
        for (int i = 0; i < size(); i++)
            c[i] -= cnum;
    }
    template <typename T1>
    void operator*=(const T1 &cnum)
    {
        for (int i = 0; i < size(); i++)
            c[i] *= cnum;
    }
    template <typename T1>
    void operator/=(const T1 &cnum)
    {
        for (int i = 0; i < size(); i++)
            c[i] /= cnum;
    }

    void reshape(int nrows_new, int ncols_new)
    {
        assert ( size() == nrows_new * ncols_new);
        nr_ = nrows_new;
        nc_ = ncols_new;
        mrank_ = std::min(nr_, nc_);
    }
    void resize(int nrows_new, int ncols_new)
    {
        reallocate_storage(nrows_new, ncols_new);
        nr_ = nrows_new;
        nc_ = ncols_new;
        set_rank_size();
        zero_out();
    }
    void resize(int nrows_new, int ncols_new, MAJOR major)
    {
        this->resize(nrows_new, ncols_new);
        major_ = major;
        set_indx_picker();
    }

    void conj()
    {
        if (!is_complex) return;
        for (int i = 0; i < this->size(); i++)
            this->c[i] = get_conj(this->c[i]);
    };

    matrix_m<T> get_transpose(bool conjugate = false)
    {
        matrix_m<T> m_trans(nc_, nr_, major_);
        for (int i = 0; i < nr_; i++)
            for (int j = 0; j < nc_; j++)
                m_trans(j, i) = this->at(i, j);
        return m_trans;
    }

    void transpose(bool conjugate = false, bool onsite_when_square_mat = true)
    {
        if (nr_ == nc_ && onsite_when_square_mat)
        {
            // handling within the memory of c pointer for square matrix
            int n = nr_;
            for (int i = 0; i != n; i++)
                for (int j = i+1; j < n; j++)
                {
                    T temp = c[i*n+j];
                    c[i*n+j] = c[j*n+i];
                    c[j*n+i] = temp;
                }
        }
        else
        {
            // NOTE: may have memory issue for large matrix as extra memory of c is required
            T* c_new = nullptr;
            if (c)
            {
                c_new = new T [nr_*nc_];
                for (int i = 0; i < nr_; i++)
                    for (int j = 0; j < nc_; j++)
                        c_new[indx_picker_2d_(nc_, j, nr_, i)] = this->at(i, j);
                delete [] c;
            }
            c = c_new;
        }
        int temp = nc_;
        nr_ = temp;
        nc_ = nr_;
        if (conjugate) conj();
    }

    matrix_m<cplx_t> to_complex()
    {
        matrix_m<cplx_t> m(nr_, nc_, major_);
        for (int i = 0; i < size_; i++)
            m.c[i] = c[i];
        return m;
    }

    matrix_m<real_t> get_real()
    {
        matrix_m<real_t> m(nr_, nc_, major_);
        for (int i = 0; i < size_; i++)
            m.c[i] = ::get_real(c[i]);
        return m;
    }

    matrix_m<real_t> get_imag()
    {
        matrix_m<real_t> m(nr_, nc_, major_);
        for (int i = 0; i < size_; i++)
            m.c[i] = ::get_imag(c[i]);
        return m;
    }
};

template <typename T1, typename T2>
inline bool same_major(const matrix_m<T1> &m1, const matrix_m<T2> &m2)
{
    return m1.major() == m2.major();
}

template <typename T1, typename T2>
inline bool same_type(const matrix_m<T1> &m1, matrix_m<T2> &m2)
{
    return std::is_same<T1, T2>::value;
}

// copy matrix elements of type T1 src to matrix of a different type
template <typename T1, typename T2>
void copy(const matrix_m<T1> &src, matrix_m<T2> &dest)
{
    assert(src.size() == dest.size() && same_major(src, dest));
    for (int i = 0; i < src.size(); i++)
        dest.c[i] = src.c[i];
}

// copy matrix elements of src to matrix of the same type
template <typename T>
void copy(const matrix_m<T> &src, matrix_m<T> &dest)
{
    assert(src.size() == dest.size() && same_major(src, dest));
    memcpy(dest.c, src.c, src.size() * sizeof(T));
}

template <typename T1, typename T2>
matrix_m<T1> operator+(const matrix_m<T1> &m1, const matrix_m<T2> &m2)
{
    assert(m1.nc() == m2.nc() && m1.nr() == m2.nr() && same_major(m1, m2));
    matrix_m<T1> sum = m1;
    sum += m2;
    return sum;
}

template <typename T1, typename T2>
inline matrix_m<T1> operator+(const matrix_m<T1> &m, const std::vector<T2> &v)
{
    assert(m.nc() == v.size());
    matrix_m<T1> mnew = m;
    mnew += v;
    return mnew;
}

template <typename T, typename T2>
inline matrix_m<T> operator+(const matrix_m<T> &m, const T2 &cnum)
{
    matrix_m<T> sum = m;
    sum += cnum;
    return sum;
}

template <typename T1, typename T2>
inline matrix_m<T2> operator+(const T1 &cnum, const matrix_m<T2> &m)
{
    return m + cnum;
}

template <typename T1, typename T2>
inline matrix_m<T1> operator-(const matrix_m<T1> &m1, const matrix_m<T2> &m2)
{
    assert(m1.nc() == m2.nc() && m1.nr() == m2.nr());
    matrix_m<T1> mnew = m1;
    mnew -= m2;
    return mnew;
}

template <typename T1, typename T2>
inline matrix_m<T1> operator-(const matrix_m<T1> &m, const std::vector<T2> &v)
{
    assert(m.nc() == v.size());
    matrix_m<T1> mnew = m;
    mnew -= v;
    return mnew;
}

template <typename T>
inline matrix_m<T> operator*(const matrix_m<T> &m1, const matrix_m<T> &m2)
{
    assert(m1.nc() == m2.nr());
    assert(m1.major() == m2.major());

    auto major = m1.major();

    matrix_m<T> prod(m1.nr(), m2.nc(), major);
    if (m1.is_row_major())
    {
        LapackConnector::gemm('N', 'N', m1.nr(), m2.nc(), m1.nc(),
                              1.0, m1.c, m1.nc(), m2.c, m2.nc(), 0.0, prod.c, prod.nc());
    }
    else
    {
        LapackConnector::gemm_f('N', 'N', m1.nr(), m2.nc(), m1.nc(),
                                1.0, m1.c, m1.nr(), m2.c, m2.nr(), 0.0, prod.c, prod.nr());
    }

    return prod;
}

// a naive implementation for integer matrix, e.g. Spglib rotation matrices
template <>
inline matrix_m<int> operator*(const matrix_m<int> &m1, const matrix_m<int> &m2)
{
    assert(m1.nc() == m2.nr());
    matrix_m<int> prod(m1.nr(), m2.nc(), m1.major());
    for (int ir = 0; ir < m1.nr(); ir++)
        for (int ik = 0; ik < m1.nc(); ik++)
            for (int ic = 0; ic < m2.nc(); ic++)
                prod(ir, ic) += m1(ir, ik) * m2(ik, ic);
    return prod;
}

template <typename T1, typename T2>
inline matrix_m<T1> operator*(const matrix_m<T1> &m, const T2 &cnum)
{
    matrix_m<T1> sum = m;
    sum *= cnum;
    return sum;
}

template <typename T1, typename T2>
inline matrix_m<T2> operator*(const T1 &cnum, const matrix_m<T2> &m)
{
    return m * cnum;
}

//! invert matrix m and store in m_inv
template <typename T>
void inverse(const matrix_m<T> &m, matrix_m<T> &m_inv)
{
    assert(m.major() == m_inv.major());
    int lwork = m.nr();
    int info = 0;
    T work[lwork];
    int ipiv[std::min(m.nr(), m.nc())];
    m_inv.resize(m.nr(), m.nc());
    copy(m, m_inv);
    // debug
    // std::cout << m_inv.size() << " " << m_inv(0, 0) << " " << m_inv(m.nr-1, m.nc-1) << std::endl;
    if (m.is_row_major())
    {
        LapackConnector::getrf(m_inv.nr(), m_inv.nc(), m_inv.c, m_inv.nr(), ipiv, info);
        LapackConnector::getri(m_inv.nr(), m_inv.c, m_inv.nr(), ipiv, work, lwork, info);
    }
    else
    {
        LapackConnector::getrf_f(m_inv.nr(), m_inv.nc(), m_inv.c, m_inv.nr(), ipiv, info);
        LapackConnector::getri_f(m_inv.nr(), m_inv.c, m_inv.nr(), ipiv, work, lwork, info);
    }
    m_inv.reshape(m_inv.nc(), m_inv.nr());
}

//! get the transpose of matrix
template <typename T>
inline matrix_m<T> transpose(const matrix_m<T> &m, bool conjugate = false)
{
    return m.get_transpose(conjugate);
}

//! compute the determinant of a matrix
template <typename T>
T get_determinant(const matrix_m<T> &m)
{
    matrix_m<T> m_copy(m);
    int lwork = m_copy.nr();
    int info = 0;
    int mrank = std::min(m_copy.nr(), m_copy.nc());
    T work[lwork];
    int ipiv[mrank];
    // debug
    if (m_copy.is_row_major())
        LapackConnector::getrf(m_copy.nr(), m_copy.nc(), m_copy.c, m_copy.nr(), ipiv, info);
    else
        LapackConnector::getrf_f(m_copy.nr(), m_copy.nc(), m_copy.c, m_copy.nr(), ipiv, info);
    T det = 1;
    for (int i = 0; i < mrank; i++)
    {
        /* std::cout << i << " " << ipiv[i] << " " << m.c[i*m.nc+i] << " "; */
        det *= (2*int(ipiv[i] == (i+1))-1) * m_copy(i, i);
    }
    return det;
}

//! convert matrix into a string
template <typename T>
std::string str(const matrix_m<T> &m)
{
    std::string s = "";
    for (int i = 0; i < m.nr(); i++)
    {
        if (m.nc() > 0) s = s + std::to_string(m(i, 0));
        for (int j = 1; j < m.nc(); j++)
            s = s + " " + std::to_string(m(i, j));
        s = s + "\n";
    }
    return s;
}

//! convert matrix into a string
template <typename T>
std::string str(const matrix_m<std::complex<T>> &m)
{
    std::string s = "";
    for (int i = 0; i != m.nr(); i++)
    {
        if (m.nc() > 0)
            s = s + "(" + std::to_string(m(i, 0).real()) + "," + std::to_string(m(i, 0).imag()) + ")";
        for (int j = 1; j != m.nc(); j++)
            s = s + " (" + std::to_string(m(i, j).real()) + "," + std::to_string(m(i, j).imag()) + ")";
        s = s + "\n";
    }
    return s;
}

//! generate a random matrix
template <typename T>
inline matrix_m<T> random(int nr, int nc, const T &lb, const T &ub, MAJOR major = MAJOR::ROW)
{
    matrix_m<T> m(nr, nc, major);
    m.randomize(lb, ub);
    return m;
}

//! generate a random symmetric matrix
template <typename T>
inline matrix_m<T> random_sy(int n, const T &lb, const T &ub)
{
    matrix_m<T> m(n, n);
    std::default_random_engine e(time(0));
    std::uniform_real_distribution<T> d(lb, ub);
    for (int i = 0; i != n; i++)
    {
        for (int j = i; j != n; j++)
        {
            m(i, j) = m(j, i) = d(e);
        }
    }
    return m;
}

//! generate a random Hermitian matrix
template <typename T>
inline matrix_m<std::complex<T>> random_he(int n, const std::complex<T> &lb, const std::complex<T> &ub)
{
    matrix_m<std::complex<T>> m(n, n);
    std::default_random_engine e(time(0));
    std::uniform_real_distribution<T> dr(lb.real(), ub.real()), di(lb.imag(), lb.imag());
    for (int i = 0; i != n; i++)
    {
        m(i, i) = dr(e);
        for (int j = i + 1; j != n; j++)
        {
            m(i, j) = std::conj(m(j, i) = std::complex<T>{dr(e), di(e)});
        }
    }
    return m;
}
