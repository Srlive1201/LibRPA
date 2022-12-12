/*
 * @file matrix_m.h
 * @brief implementing a matrix template class with support of column-major or row-major style
 */
#pragma once
#include <cassert>
#include <iomanip>
#include <vector>
#include <cmath>
#include <ctime>
#include <random>
#include <valarray>
#include <memory>
#include <functional>
#include <ostream>
#include "base_utility.h"
#include "lapack_connector.h"
#include "vec.h"

//! type alias for functors to map 2D array index to flatten array index
typedef std::function<int(const int& nr, const int& ir, const int& nc, const int& ic)> Indx_picker_2d;
// major dependent flatted index of 2D matrix
//! row-major 2D index picker
inline int flatid_rm(const int &nr, const int &ir, const int &nc, const int &ic) { return ir*nc+ic; }
//! column-major 2D index picker
inline int flatid_cm(const int &nr, const int &ir, const int &nc, const int &ic) { return ic*nr+ir; }

enum MAJOR
{
    ROW = 0,
    COL = 1
};
const Indx_picker_2d Indx_pickers_2d[]
{
    flatid_rm,
    flatid_cm
};

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

    const Indx_picker_2d picker = Indx_pickers_2d[int(!row_major)];

    if (nr&&nc)
    {
        for (int ir = 0; ir < nr_; ir++)
            for (int ic = 0; ic < std::min(size_t(nc_), nested_vector[ir].size()); ic++)
                c[picker(nr, ir, nc, ic)] = nested_vector[ir][ic];
    }
}

template <typename T>
class matrix_m
{
public:
    //! Flag of whether the instantialized matrix class is complex-typed
    static const bool is_complex = is_complex_t<T>::value;
    static const bool is_double = std::is_same<T, double>::value;
private:
    //! The number of rows
    int nr_;
    //! The number of columns
    int nc_;
    //! The leading dimension of the matrix
    int ld_;
    //! Major of matrix data storage, either row- or column-major
    MAJOR major_;
    //! The maximal rank of the matrix, i.e. min(nr, nc)
    int mrank_;
    //! The number of matrix elements, i.e. nr * nc
    int size_;
    //! Function to extract the data in the memory from 2D indexing
    Indx_picker_2d indx_picker_2d_;

public:
    // std::shared_ptr<std::valarray<T>> data; // to conform with LibRI tensor
    T* c;

private:
    void assign_value(const T &v)
    {
        for (int i = 0; i < size_; i++)
            c[i] = v;
    }
    void assign_value(const T* const pv, MAJOR major_pv)
    {
        if (major_ == major_pv)
        {
            for (int i = 0; i < size_; i++) c[i] = pv[i];
        }
        else
        {
            const Indx_picker_2d pv_picker = Indx_pickers_2d[major_pv];
            for (int ir = 0; ir < nr_; ir++)
                for (int ic = 0; ic < nc_; ic++)
                    this->at(ir, ic) = pv[pv_picker(nr_, ir, nc_, ic)];
        }
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
        indx_picker_2d_ = Indx_pickers_2d[major_];
        ld_ = (major_ == MAJOR::ROW ? nc_ : nr_);
    }

public:
    using type = T;
    using real_t = typename to_real<T>::type;
    using cplx_t = typename to_cplx<T>::type;

    // constructors
    matrix_m() : nr_(0), nc_(0), ld_(0), major_(MAJOR::ROW), c(nullptr)
    {
        set_rank_size();
        set_indx_picker();
    }
    matrix_m(const int &nrows, const int &ncols, MAJOR major = MAJOR::ROW): nr_(nrows), nc_(ncols), ld_(nrows), major_(major), c(nullptr)
    {
        set_rank_size();
        set_indx_picker();
        allocate_storage();
        zero_out();
    }
    matrix_m(const int &nrows, const int &ncols, const T * const valarr, MAJOR major = MAJOR::ROW, MAJOR major_valarr = MAJOR::ROW): nr_(nrows), nc_(ncols), ld_(nrows), major_(major), c(nullptr)
    {
        set_rank_size();
        set_indx_picker();
        allocate_storage();
        zero_out();
        assign_value(valarr, major_valarr);
    }
    // constructor from a nested vector
    matrix_m(const std::vector<std::vector<T>> &nested_vector, MAJOR major = MAJOR::ROW): nr_(0), nc_(0), ld_(0), major_(major), mrank_(0), size_(0), c(nullptr)
    {
        get_nr_nc_from_nested_vector(nested_vector, nr_, nc_);
        set_rank_size();
        set_indx_picker();
        allocate_storage();
        zero_out();
        expand_nested_vector_to_pointer(nested_vector, nr_, nc_, c, is_row_major());
    }

    // copy constructor
    matrix_m(const matrix_m<T> &m): nr_(m.nr_), nc_(m.nc_), ld_(m.ld_), major_(m.major_), mrank_(0), size_(0), c(nullptr)
    {
        set_rank_size();
        set_indx_picker();
        allocate_storage();
        memcpy(c, m.c, nr_*nc_*sizeof(T));
    }

    matrix_m(matrix_m &&m) : nr_(m.nr_), nc_(m.nc_), ld_(m.ld_), major_(m.major_), mrank_(m.mrank_), size_(m.size_), indx_picker_2d_(m.indx_picker_2d_), c(m.c)
    {
        m.nr_ = m.nc_ = m.ld_ = m.mrank_ = m.size_ = 0;
        m.c = nullptr;
    }
    // destructor
    ~matrix_m()
    {
        nc_ = nr_ = ld_ = mrank_ = size_ = 0;
        deallocate_storage();
    }

    void zero_out()
    {
        assign_value(T(0));
    }

    void clear()
    {
        this->resize(0, 0);
    }

    void set_diag(const T& v)
    {
        for (int i = 0; i < mrank_; i++)
            this->at(i, i) = v;
    }

    void scale_row(int irow, const T &scale)
    {
        for (int ic = 0; ic < nc_; ic++)
            this->at(irow, ic) *= scale;
    }

    void scale_col(int icol, const T &scale)
    {
        for (int ir = 0; ir < nr_; ir++)
            this->at(ir, icol) *= scale;
    }

    //! Randomize the matrix elements with lower and upper bound and symmetry constraint
    void randomize(const T &lb = 0, const T &ub = 1, bool symmetrize = false, bool hermitian = false)
    {
        static std::default_random_engine e(time(0));
        std::uniform_real_distribution<real_t> dr(::get_real(lb), ::get_real(ub));
        if (matrix_m<T>::is_complex)
        {
            std::uniform_real_distribution<real_t> di(::get_imag(lb), ::get_imag(ub));
            if (symmetrize)
            {
                for (int i = 0; i != mrank_; i++)
                {
                    join_re_im((*this)(i, i), dr(e), di(e));
                    for (int j = i + 1; j < mrank_; j++)
                    {
                        join_re_im((*this)(i, j), dr(e), di(e));
                        (*this)(j, i) = (*this)(i, j);
                    }
                }
            }
            else if (hermitian)
            {
                for (int i = 0; i != mrank_; i++)
                {
                    join_re_im((*this)(i, i), dr(e), real_t(0.0));
                    for (int j = i + 1; j < mrank_; j++)
                    {
                        join_re_im((*this)(i, j), dr(e), di(e));
                        (*this)(j, i) = ::get_conj((*this)(i, j));
                    }
                }
            }
            else
                for (int i = 0; i != size_; i++)
                    join_re_im(this->c[i], dr(e), di(e));
        }
        else
        {
            if (symmetrize)
                for (int i = 0; i != mrank_; i++)
                {
                    this->at(i, i) = dr(e);
                    for (int j = i + 1; j < mrank_; j++)
                    {
                        this->at(i, j) = dr(e);
                        this->at(j, i) = this->at(i, j);
                    }
                }
            else
            {
                for (int i = 0; i != size_; i++)
                    this->c[i] = dr(e);
            }
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
    T &operator()(const int ir, const int ic) { return this->at(ir, ic); }
    const T &operator()(const int ir, const int ic) const { return this->at(ir, ic); }
    T &at(const int ir, const int ic) { return c[indx_picker_2d_(nr_, ir, nc_, ic)]; }
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
    matrix_m<T>& resize(int nrows_new, int ncols_new)
    {
        reallocate_storage(nrows_new, ncols_new);
        nr_ = nrows_new;
        nc_ = ncols_new;
        set_rank_size();
        zero_out();
        return *this;
    }

    matrix_m<T>& resize(int nrows_new, int ncols_new, MAJOR major)
    {
        major_ = major;
        set_indx_picker();
        return this->resize(nrows_new, ncols_new);
    }

    void conj()
    {
        if (!is_complex) return;
        for (int i = 0; i < this->size(); i++)
            this->c[i] = get_conj(this->c[i]);
    };

    matrix_m<T> get_transpose(bool conjugate = false) const
    {
        matrix_m<T> m_trans(nc_, nr_, major_);
        for (int i = 0; i < nr_; i++)
            for (int j = 0; j < nc_; j++)
                m_trans(j, i) = this->at(i, j);
        if (conjugate) m_trans.conj();
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
                              T(1.0), m1.c, m1.nc(), m2.c, m2.nc(), T(0.0), prod.c, prod.nc());
    }
    else
    {
        LapackConnector::gemm_f('N', 'N', m1.nr(), m2.nc(), m1.nc(),
                                T(1.0), m1.c, m1.nr(), m2.c, m2.nr(), T(0.0), prod.c, prod.nr());
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

template <typename T>
std::ostream& operator<<(std::ostream& os, const matrix_m<T> &m)
{
    os << str(m);
    return os;
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
inline matrix_m<T> random_sy(int n, const T &lb, const T &ub, MAJOR major = MAJOR::ROW)
{
    matrix_m<T> m(n, n, major);
    m.randomize(lb, ub, true, false);
    return m;
}

//! generate a random Hermitian matrix
template <typename T>
inline matrix_m<std::complex<T>> random_he(int n, const std::complex<T> &lb, const std::complex<T> &ub, MAJOR major = MAJOR::ROW)
{
    matrix_m<std::complex<T>> m(n, n, major);
    m.randomize(lb, ub, false, true);
    return m;
}

//! generate a random unitary matrix, implemented as using ?heev upon hermitian matrix
template <typename T>
inline matrix_m<std::complex<T>> random_unitary(int n, MAJOR major = MAJOR::ROW)
{
    int info;
    // always use column major here to ensure that the eigvec construction is correct
    auto mat = random_he<T>(n, {-1.0, -1.0}, {1.0, 1.0}, MAJOR::COL);
    const char *name = mat.is_double? "ZHETRD" : "CHETRD";
    // printf("%s\n", name);
    int nb = LapackConnector::ilaenv(1, name, "VU", n, -1, -1, -1);
    if (nb <= 1) nb = std::max(1, n);
    const int lwork = n * (nb + 1);

    T *w = new T[n];
    std::complex<T> *work = new std::complex<T>[lwork];
    T *rwork = new T[std::max(1, 3 * n - 2)];

    LapackConnector::heev_f('V', 'U', n, mat.c, mat.nr(), w, work, lwork, rwork, info);
    delete [] rwork, w, work;
    return matrix_m<std::complex<T>>(n, n, mat.c, major, MAJOR::COL);
}

//! generate a random Hermitian matrix with selected eigenvalues
template <typename T>
inline matrix_m<std::complex<T>> random_he_selected_ev(int n, const vector<T> &evs, MAJOR major = MAJOR::ROW)
{
    const int nvec = evs.size();
    if (nvec == 0) return matrix_m<std::complex<T>>{0, 0, major};

    auto eigvec = random_unitary<T>(n, major);
    auto eigvec_H = eigvec.get_transpose(true);
    for (int ie = 0; ie != n; ie++)
        eigvec_H.scale_row(ie, ie < nvec? evs[ie]: +0.0);
    return eigvec * eigvec_H;
}

//! print the matrix in matrix-market style
template <typename T, typename Treal = typename to_real<T>::type>
void print_matrix_mm(const matrix_m<T> &mat, ostream &os, Treal threshold = 1e-15, bool row_first = true)
{
    const int nr = mat.nr(), nc = mat.nc();
    size_t nnz = 0;
    for (int i = 0; i != mat.size(); i++)
    {
        if (fabs(mat.c[i]) > threshold)
            nnz++;
    }
    os << "%%MatrixMarket matrix coordinate "
       << (is_complex<T>()? "complex" : "real") << " general" << endl
       << "%" << endl;
    os << nr << " " << nc << " " << nnz << endl;
        if (row_first)
        {
            for (int i = 0; i != nr; i++)
                for (int j = 0; j != nc; j++)
                {
                    auto &v = mat(i, j);
                    if (fabs(v) > threshold)
                    {
                        os << i + 1 << " " << j + 1 << " " << showpoint << scientific << setprecision(15) << get_real(v);
                        if (is_complex<T>())
                            os << " " << showpoint << scientific << setprecision(15) << get_imag(v);
                        os << endl;
                    }
                }
        }
        else
        {
            for (int j = 0; j != nc; j++)
                for (int i = 0; i != nr; i++)
                {
                    auto &v = mat(i, j);
                    if (fabs(v) > threshold)
                    {
                        os << i + 1 << " " << j + 1 << " " << showpoint << scientific << setprecision(15) << get_real(v);
                        if (is_complex<T>())
                            os << " " << showpoint << scientific << setprecision(15) << get_imag(v);
                        os << endl;
                    }
                }
        }
}
