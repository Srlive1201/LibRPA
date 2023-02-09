/*
 * @file matrix_m.h
 * @brief implementing a matrix template class with support of column-major or row-major style
 */
#pragma once
#include <array>
#include <cassert>
#include <iomanip>
#include <vector>
#include <cmath>
#include <ctime>
#include <random>
#include <valarray>
#include <memory>
#include <functional>
#include <utility>
#include <ostream>
#include "base_utility.h"
#include "lapack_connector.h"
#include "vec.h"

//! type alias for functors to map 2D array index to flatten array index
typedef std::function<int(const int& nr, const int& ir, const int& nc, const int& ic)> Indx_picker_2d;
typedef std::function<void(const int& indx, const int& nr, int& ir, const int& nc, int& ic)> Indx_folder_2d;
// major dependent flatted index of 2D matrix
//! row-major 2D index picker
inline int flatid_rm(const int &nr, const int &ir, const int &nc, const int &ic) { return ir*nc+ic; }
//! column-major 2D index picker
inline int flatid_cm(const int &nr, const int &ir, const int &nc, const int &ic) { return ic*nr+ir; }

// major dependent folding of 1D index to 2D index of matrix
//! row-major 1D index folder to 2D
inline void foldid_rm(const int& indx, const int &nr, int &ir, const int &nc, int &ic) { ir = indx / nc; ic = indx % nc; }
//! column-major 1D index folder to 2D
inline void foldid_cm(const int& indx, const int &nr, int &ir, const int &nc, int &ic) { ir = indx % nr; ic = indx / nr; }


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
const Indx_folder_2d Indx_folders_2d[]
{
    foldid_rm,
    foldid_cm
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

//! struct to store the shared data and dimension
template <typename T>
struct matrix_data
{
public:
    std::shared_ptr<std::array<int, 2>> dim = nullptr;
    std::shared_ptr<std::valarray<T>> data = nullptr;

private:
    int size_;
    int mrank_;
    void set_size_mrank()
    {
        mrank_ = std::min((*dim)[0], (*dim)[1]);
        size_ = (*dim)[0] * (*dim)[1];
    }

public:
    matrix_data()
    {
        dim = std::make_shared<std::array<int, 2>>();
        mrank_ = size_ = 0;
    }
    matrix_data(const std::array<int, 2> &dim_in)
    {
        dim = std::make_shared<std::array<int, 2>>(dim_in);
        set_size_mrank();
        if (size())
            data = std::make_shared<std::valarray<T>>(size());
    }
    matrix_data(const std::array<int, 2> &dim_in, const std::valarray<T> &data_in)
    {
        dim = std::make_shared<std::array<int, 2>>(dim_in);
        data = std::make_shared<std::valarray<T>>(data_in);
        set_size_mrank();
    }
    matrix_data(const std::array<int, 2> &dim_in, const T* const pdata)
    {
        dim = std::make_shared<std::array<int, 2>>(dim_in);
        set_size_mrank();
        if (size() > 0)
            data = std::make_shared<std::valarray<T>>(pdata, size());
    }
    matrix_data(const std::array<int, 2> &dim_in,
                const std::shared_ptr<std::valarray<T>> &data_in)
        : data(data_in)
    {
        dim = std::make_shared<std::array<int, 2>>(dim_in);
        set_size_mrank();
    }
    matrix_data(const std::shared_ptr<std::array<int, 2>> &dim_in,
                const std::shared_ptr<std::valarray<T>> &data_in)
        : dim(dim_in), data(data_in)
    {
        set_size_mrank();
    }

    inline int nr() const { return (*dim)[0]; }
    inline int nc() const { return (*dim)[1]; }
    inline int size() const { return size_; }
    inline int mrank() const { return mrank_; }

    // data access
    inline T &operator[](const int i) { return (*data)[i]; }
    inline const T &operator[](const int i) const { return (*data)[i]; }

    inline T* ptr() { return &(*this->data)[0]; }
    inline const T* ptr() const { return &(*this->data)[0]; }
    inline std::shared_ptr<std::valarray<T>> sptr() { return data; }
    inline const std::shared_ptr<std::valarray<T>> sptr() const { return data; }

    void resize(const std::array<int, 2>& dim_new)
    {
        *dim = dim_new;
        set_size_mrank();
        if (data == nullptr)
        {
            if (size() > 0)
                data = std::make_shared<std::valarray<T>>(0, size());
        }
        else
        {
            if (size() > 0)
                data->resize(size());
            else
                data->resize(0);
        }
    }

    void reshape(const std::array<int, 2>& dim_new)
    {
        assert(dim_new[0] * dim_new[1] == size());
        *dim = dim_new;
        set_size_mrank();
    }

    void operator+=(const matrix_data<T> &mdata)
    {
        if (data != nullptr && mdata.data != nullptr)
            *data += *(mdata.data);
    }

    void operator-=(const matrix_data<T> &mdata)
    {
        if (data != nullptr && mdata.data != nullptr)
            *data -= *(mdata.data);
    }
};

template <typename T>
bool operator==(const matrix_data<T> &mdata1, const matrix_data<T> &mdata2)
{
    if (mdata1.nr() != mdata2.nr() || mdata1.nc() != mdata2.nc())
        return false;
    if (mdata1.sptr() == nullptr || mdata2.sptr() == nullptr)
        return false;
    using real_t = typename to_real<T>::type;
    const bool is_double = std::is_same<T, double>::value;
    const real_t thres = is_double ? 1e-14 : 1e-6;
    auto diff = std::abs(*(mdata1.data) - *(mdata2.data));
    for (int i = 0; i < mdata1.size(); i++)
        if (::get_real(diff[i]) > thres)
            return false;
    return true;
}

template <typename T>
class matrix_m
{
public:
    //! Flag of whether the instantialized matrix class is complex-typed
    static const bool is_complex = is_complex_t<T>::value;
    static const bool is_double = std::is_same<T, double>::value;
private:
    //! Major of matrix data storage, either row- or column-major
    MAJOR major_;
    //! The leading dimension of the matrix
    int ld_;
    //! Function to extract the data in the memory from 2D indexing
    Indx_picker_2d indx_picker_2d_;
    //! Function to convert 1D index to 2D index
    Indx_folder_2d indx_folder_2d_;

public:
    matrix_data<T> dataobj;

private:
    void assign_value(const T &v)
    {
        if (size() > 0)
            *(dataobj.data) = v;
    }
    void assign_value(const T* const pv, MAJOR major_pv)
    {
        if (major_ == major_pv)
        {
            for (int i = 0; i < size(); i++) dataobj[i] = pv[i];
        }
        else
        {
            const Indx_picker_2d pv_picker = Indx_pickers_2d[major_pv];
            for (int ir = 0; ir < nr(); ir++)
                for (int ic = 0; ic < nc(); ic++)
                    this->at(ir, ic) = pv[pv_picker(nr(), ir, nc(), ic)];
        }
    }

    void set_ld()
    {
        ld_ = (major_ == MAJOR::ROW ? nc() : nr());
    }

    void set_indx_picker_folder()
    {
        indx_picker_2d_ = Indx_pickers_2d[major_];
        indx_folder_2d_ = Indx_folders_2d[major_];
    }

public:
    using type = T;
    using real_t = typename to_real<T>::type;
    using cplx_t = typename to_cplx<T>::type;

    // constructors
    matrix_m() : major_(MAJOR::ROW), ld_(0), dataobj()
    {
        set_ld();
        set_indx_picker_folder();
    }
    matrix_m(const int &nrows, const int &ncols, MAJOR major = MAJOR::ROW): major_(major), ld_(0), dataobj({nrows, ncols})
    {
        set_ld();
        set_indx_picker_folder();
        zero_out();
    }
    matrix_m(const int &nrows, const int &ncols, const T * const parr, MAJOR major = MAJOR::ROW): major_(major), ld_(0), dataobj({nrows, ncols}, parr)
    {
        set_ld();
        set_indx_picker_folder();
    }
    matrix_m(const int &nrows, const int &ncols, const T * const parr, MAJOR major_arr, MAJOR major): major_(major), ld_(0), dataobj({nrows, ncols})
    {
        set_ld();
        set_indx_picker_folder();
        zero_out();
        assign_value(parr, major_arr);
    }
    matrix_m(const int &nrows, const int &ncols, const std::shared_ptr<std::valarray<T>> &valarr, MAJOR major_valarr): major_(major_valarr), ld_(0), dataobj({nrows, ncols}, valarr)
    {
        set_ld();
        set_indx_picker_folder();
    }
    // constructor from a nested vector
    matrix_m(const std::vector<std::vector<T>> &nested_vector, MAJOR major = MAJOR::ROW): major_(major), ld_(0), dataobj()
    {
        int nr_, nc_;
        get_nr_nc_from_nested_vector(nested_vector, nr_, nc_);
        dataobj.resize({nr_, nc_});
        set_ld();
        set_indx_picker_folder();
        zero_out();
        expand_nested_vector_to_pointer(nested_vector, nr_, nc_, dataobj.ptr(), is_row_major());
    }

    // copy constructor
    matrix_m(const matrix_m<T> &m) = default;

    matrix_m(matrix_m &&m) = default;

    // destructor
    ~matrix_m() {}

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
        for (int i = 0; i < mrank(); i++)
            this->at(i, i) = v;
    }

    void scale_row(int irow, const T &scale)
    {
        for (int ic = 0; ic < nc(); ic++)
            this->at(irow, ic) *= scale;
    }

    void scale_col(int icol, const T &scale)
    {
        for (int ir = 0; ir < nr(); ir++)
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
                for (int i = 0; i != mrank(); i++)
                {
                    join_re_im((*this)(i, i), dr(e), di(e));
                    for (int j = i + 1; j < mrank(); j++)
                    {
                        join_re_im((*this)(i, j), dr(e), di(e));
                        (*this)(j, i) = (*this)(i, j);
                    }
                }
            }
            else if (hermitian)
            {
                for (int i = 0; i != mrank(); i++)
                {
                    join_re_im((*this)(i, i), dr(e), real_t(0.0));
                    for (int j = i + 1; j < mrank(); j++)
                    {
                        join_re_im((*this)(i, j), dr(e), di(e));
                        (*this)(j, i) = ::get_conj((*this)(i, j));
                    }
                }
            }
            else
                for (int i = 0; i != size(); i++)
                    join_re_im(this->dataobj[i], dr(e), di(e));
        }
        else
        {
            if (symmetrize)
                for (int i = 0; i != mrank(); i++)
                {
                    this->at(i, i) = dr(e);
                    for (int j = i + 1; j < mrank(); j++)
                    {
                        this->at(i, j) = dr(e);
                        this->at(j, i) = this->at(i, j);
                    }
                }
            else
            {
                for (int i = 0; i != size(); i++)
                    this->dataobj[i] = dr(e);
            }
        }
    }

    // access to private variables
    int nr() const { return dataobj.nr(); }
    int nc() const { return dataobj.nc(); }
    // leading dimension
    int ld() const { return ld_; }
    MAJOR major() const { return major_; }
    int size() const { return dataobj.size(); }
    int mrank() const { return dataobj.mrank(); }

    bool is_row_major() const { return major_ == MAJOR::ROW; }
    bool is_col_major() const { return major_ == MAJOR::COL; }

    // indexing
    T &operator()(const int ir, const int ic) { return this->at(ir, ic); }
    const T &operator()(const int ir, const int ic) const { return this->at(ir, ic); }
    T &at(const int ir, const int ic) { return dataobj[indx_picker_2d_(nr(), ir, nc(), ic)]; }
    const T &at(const int ir, const int ic) const { return dataobj[indx_picker_2d_(nr(), ir, nc(), ic)]; }

    T* ptr() { return dataobj.ptr(); }
    const T* ptr() const { return dataobj.ptr(); }
    std::shared_ptr<std::valarray<T>> sptr() { return dataobj.sptr(); }
    const std::shared_ptr<std::valarray<T>> sptr() const { return dataobj.sptr(); }

    matrix_m<T> copy() const
    {
        matrix_m<T> m_copy(nr(), nc(), ptr(), major_);
        return m_copy;
    }

    // swap major
    void swap_to_row_major()
    {
        if (is_col_major())
        {
            int nr_old = nr(), nc_old = nc();
            this->transpose();
            reshape(nr_old, nc_old);
            major_ = MAJOR::ROW;
            set_ld();
            set_indx_picker_folder();
        }
    }
    void swap_to_col_major()
    {
        if (is_row_major())
        {
            int nr_old = nr(), nc_old = nc();
            this->transpose();
            reshape(nr_old, nc_old);
            major_ = MAJOR::COL;
            set_ld();
            set_indx_picker_folder();
        }
    }

    // assignment copy
    matrix_m<T> & operator=(const matrix_m<T> &m) = default;
    matrix_m<T> & operator=(matrix_m<T> &&m) = default;

    // operator overload
    matrix_m<T> & operator=(const T &cnum)
    {
        (dataobj.data)->fill(cnum);
        return *this;
    }

    template <typename T1>
    void operator+=(const matrix_m<T1> &m)
    {
        assert(size() == m.size() && major() == m.major());
        for (int i = 0; i < size(); i++)
            dataobj[i] += m.dataobj[i];
    }
    template <typename T1>
    void operator+=(const std::vector<T1> &v)
    {
        assert(nc() == v.size());
        for (int i = 0; i < nr(); i++)
            for (int ic = 0; ic < nc(); ic++)
                this->at(i, ic) += v[ic];
    }
    template <typename T1>
    void operator+=(const T1 &cnum)
    {
        for (int i = 0; i < size(); i++)
            dataobj[i] += cnum;
    }

    template <typename T1>
    void operator-=(const matrix_m<T1> &m)
    {
        assert(size() == m.size() && major() == m.major());
        for (int i = 0; i < size(); i++)
            dataobj[i] -= m.dataobj[i];
    }
    template <typename T1>
    void operator-=(const std::vector<T1> &v)
    {
        assert(nc() == v.size());
        for (int i = 0; i < nr(); i++)
            for (int ic = 1; ic < nc(); ic++)
                this->at(i, ic) += v[ic];
    }
    template <typename T1>
    void operator-=(const T1 &cnum)
    {
        for (int i = 0; i < size(); i++)
            dataobj[i] -= cnum;
    }
    template <typename T1>
    void operator*=(const T1 &cnum)
    {
        for (int i = 0; i < size(); i++)
            dataobj[i] *= cnum;
    }
    template <typename T1>
    void operator/=(const T1 &cnum)
    {
        for (int i = 0; i < size(); i++)
            dataobj[i] /= cnum;
    }

    void reshape(int nrows_new, int ncols_new)
    {
        assert(size() == nrows_new * ncols_new);
        dataobj.reshape({nrows_new, ncols_new});
    }
    matrix_m<T>& resize(int nrows_new, int ncols_new)
    {
        dataobj.resize({nrows_new, ncols_new});
        set_ld();
        zero_out();
        return *this;
    }

    matrix_m<T>& resize(int nrows_new, int ncols_new, MAJOR major)
    {
        major_ = major;
        set_indx_picker_folder();
        return this->resize(nrows_new, ncols_new);
    }

    matrix_m<T>& conj()
    {
        if (!matrix_m<T>::is_complex) return *this;
        for (int i = 0; i < this->size(); i++)
            dataobj[i] = get_conj(dataobj[i]);
        return *this;
    };

    matrix_m<T> get_transpose(bool conjugate = false) const
    {
        matrix_m<T> m_trans(nc(), nr(), major_);
        for (int i = 0; i < nr(); i++)
            for (int j = 0; j < nc(); j++)
                m_trans(j, i) = this->at(i, j);
        if (conjugate) m_trans.conj();
        return m_trans;
    }

    matrix_m<T>& transpose(bool conjugate = false)
    {
        if (nr() == nc())
        {
            // handling within the memory of c pointer for square matrix
            int n = nr();
            for (int i = 0; i != n; i++)
                for (int j = i + 1; j < n; j++)
                {
                    T temp = dataobj[i * n + j];
                    dataobj[i * n + j] = dataobj[j * n + i];
                    dataobj[j * n + i] = temp;
                }
        }
        else
        {
            // NOTE: may have memory issue for large matrix as extra memory is requested
            if (size() > 0)
            {
                std::valarray<T> a_new(size());
                for (int i = 0; i < nr(); i++)
                    for (int j = 0; j < nc(); j++)
                        a_new[indx_picker_2d_(nc(), j, nr(), i)] = this->at(i, j);
                *(dataobj.data) = a_new;
            }
        }
        if (conjugate) conj();
        return *this;
    }

    matrix_m<cplx_t> to_complex() const
    {
        matrix_m<cplx_t> m(nr(), nc(), major_);
        for (int i = 0; i < size(); i++)
            m.dataobj[i] = dataobj[i];
        return m;
    }

    matrix_m<real_t> get_real() const
    {
        matrix_m<real_t> m(nr(), nc(), major_);
        for (int i = 0; i < size(); i++)
            m.dataobj[i] = ::get_real(dataobj[i]);
        return m;
    }

    matrix_m<real_t> get_imag() const
    {
        matrix_m<real_t> m(nr(), nc(), major_);
        for (int i = 0; i < size(); i++)
            m.dataobj[i] = ::get_imag(dataobj[i]);
        return m;
    }

    int count_absle(real_t thres) const
    {
        int count = 0;
        for (int i = 0; i < size(); i++)
            if (std::abs(dataobj[i]) <= thres) count++;
        return count;
    }

    int count_absge(real_t thres) const
    {
        int count = 0;
        for (int i = 0; i < size(); i++)
            if (std::abs(dataobj[i]) >= thres) count++;
        return count;
    }

    void argmin(int &ir, int &ic) const
    {
        // for complex type, compare the absolute value
        const std::function<real_t(T)> filter =
            matrix_m<T>::is_complex ? [](T a) { return std::abs(a); } : [](T a) { return ::get_real(a); };
        real_t value = std::numeric_limits<real_t>::max();
        int indx = 0;
        for (int i = 0; i < size(); i++)
            if (filter(dataobj[i]) < value)
            {
                indx = i;
                value = filter(dataobj[i]);
            }
        indx_folder_2d_(indx, nr(), ir, nc(), ic);
    }

    void argmax(int &ir, int &ic) const
    {
        // for complex type, compare the absolute value
        const std::function<real_t(T)> filter =
            matrix_m<T>::is_complex ? [](T a) { return std::abs(a); } : [](T a) { return ::get_real(a); };
        real_t value = std::numeric_limits<real_t>::min();
        int indx = 0;
        for (int i = 0; i < size(); i++)
            if (filter(dataobj[i]) > value)
            {
                indx = i;
                value = filter(dataobj[i]);
            }
        indx_folder_2d_(indx, nr(), ir, nc(), ic);
    }

    real_t min() const
    {
        // for complex type, compare the absolute value
        const std::function<real_t(T)> filter =
            matrix_m<T>::is_complex ? [](T a) { return std::abs(a); } : [](T a) { return ::get_real(a); };
        real_t value = std::numeric_limits<real_t>::max();
        for (int i = 0; i < size(); i++)
            if (filter(dataobj[i]) < value)
            {
                value = filter(dataobj[i]);
            }
        return value;
    }

    real_t max() const
    {
        // for complex type, compare the absolute value
        const std::function<real_t(T)> filter =
            matrix_m<T>::is_complex ? [](T a) { return std::abs(a); } : [](T a) { return ::get_real(a); };
        real_t value = std::numeric_limits<real_t>::min();
        for (int i = 0; i < size(); i++)
            if (filter(dataobj[i]) > value)
            {
                value = filter(dataobj[i]);
            }
        return value;
    }

    real_t absmin() const
    {
        real_t value = 0;
        for (int i = 0; i < size(); i++)
            value = std::min(value, std::abs(dataobj[i]));
        return value;
    }

    real_t absmax() const
    {
        real_t value = 0;
        for (int i = 0; i < size(); i++)
            value = std::max(value, std::abs(dataobj[i]));
        return value;
    }
};

//! typedef of double matrix
typedef matrix_m<double> Matd;
//! typedef of complex double matrix
typedef matrix_m<std::complex<double>> Matz;

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
        dest.dataobj[i] = src.dataobj[i];
}

// copy matrix elements of src to matrix of the same type
template <typename T>
void copy(const matrix_m<T> &src, matrix_m<T> &dest)
{
    assert(src.size() == dest.size() && same_major(src, dest));
    memcpy(dest.ptr(), src.ptr(), src.size() * sizeof(T));
}

template <typename T1, typename T2>
matrix_m<T1> operator+(const matrix_m<T1> &m1, const matrix_m<T2> &m2)
{
    assert(m1.nc() == m2.nc() && m1.nr() == m2.nr() && same_major(m1, m2));
    auto sum = m1.copy();
    sum += m2;
    return sum;
}

template <typename T1, typename T2>
inline matrix_m<T1> operator+(const matrix_m<T1> &m, const std::vector<T2> &v)
{
    assert(m.nc() == v.size());
    auto mnew = m.copy();
    mnew += v;
    return mnew;
}

template <typename T, typename T2>
inline matrix_m<T> operator+(const matrix_m<T> &m, const T2 &cnum)
{
    auto sum = m.copy();
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
    auto mnew = m1.copy();
    mnew -= m2;
    return mnew;
}

template <typename T1, typename T2>
inline matrix_m<T1> operator-(const matrix_m<T1> &m, const std::vector<T2> &v)
{
    assert(m.nc() == v.size());
    auto mnew = m.copy();
    mnew -= v;
    return mnew;
}

//! memory-retaining matrix multiply
template <typename T>
void matmul(const matrix_m<T> &m1, const matrix_m<T> &m2, matrix_m<T> &prod)
{
    assert(m1.nc() == m2.nr());
    assert(m1.nr() == prod.nr());
    assert(m2.nc() == prod.nc());
    assert(m1.major() == m2.major() && m1.major() == prod.major());

    if (m1.is_row_major())
    {
        LapackConnector::gemm('N', 'N', m1.nr(), m2.nc(), m1.nc(),
                              T(1.0), m1.ptr(), m1.nc(), m2.ptr(), m2.nc(), T(0.0), prod.ptr(), prod.nc());
    }
    else
    {
        LapackConnector::gemm_f('N', 'N', m1.nr(), m2.nc(), m1.nc(),
                                T(1.0), m1.ptr(), m1.nr(), m2.ptr(), m2.nr(), T(0.0), prod.ptr(), prod.nr());
    }
}

template <typename T>
inline matrix_m<T> operator*(const matrix_m<T> &m1, const matrix_m<T> &m2)
{
    assert(m1.nc() == m2.nr());
    assert(m1.major() == m2.major());

    auto major = m1.major();

    matrix_m<T> prod(m1.nr(), m2.nc(), major);
    matmul<T>(m1, m2, prod);
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
    auto sum = m.copy();
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
        LapackConnector::getrf(m_inv.nr(), m_inv.nc(), m_inv.ptr(), m_inv.nr(), ipiv, info);
        LapackConnector::getri(m_inv.nr(), m_inv.ptr(), m_inv.nr(), ipiv, work, lwork, info);
    }
    else
    {
        LapackConnector::getrf_f(m_inv.nr(), m_inv.nc(), m_inv.ptr(), m_inv.nr(), ipiv, info);
        LapackConnector::getri_f(m_inv.nr(), m_inv.ptr(), m_inv.nr(), ipiv, work, lwork, info);
    }
    m_inv.reshape(m_inv.nc(), m_inv.nr());
}

//! get the transpose of matrix
template <typename T>
inline matrix_m<T> transpose(const matrix_m<T> &m, bool conjugate = false)
{
    return m.get_transpose(conjugate);
}

//! get the conjugate of matrix
template <typename T>
inline matrix_m<T> conj(const matrix_m<T> &m)
{
    return m.copy().conj();
}

//! compute the determinant of a matrix
template <typename T>
T get_determinant(const matrix_m<T> &m)
{
    matrix_m<T> m_copy = m.copy();
    int lwork = m_copy.nr();
    int info = 0;
    int mrank = std::min(m_copy.nr(), m_copy.nc());
    T work[lwork];
    int ipiv[mrank];
    // debug
    if (m_copy.is_row_major())
        LapackConnector::getrf(m_copy.nr(), m_copy.nc(), m_copy.ptr(), m_copy.nr(), ipiv, info);
    else
        LapackConnector::getrf_f(m_copy.nr(), m_copy.nc(), m_copy.ptr(), m_copy.nr(), ipiv, info);
    T det = 1;
    for (int i = 0; i < mrank; i++)
    {
        /* std::cout << i << " " << ipiv[i] << " " << m.c[i*m.nc+i] << " "; */
        det *= (2*int(ipiv[i] == (i+1))-1) * m_copy(i, i);
    }
    return det;
}

template <typename T>
matrix_m<std::complex<T>> power_hemat(matrix_m<std::complex<T>> &mat,
                                      T power, bool filter_original,
                                      const T &threshold = -1.e5)
{
    assert (mat.nr() == mat.nc());
    const char jobz = 'V';
    const char uplo = 'U';
    const int n = mat.nr();
    int nb;
    if (matrix_m<T>::is_double)
        nb = LapackConnector::ilaenv(1, "zheev", "VU", n, -1, -1, -1);
    else
        nb = LapackConnector::ilaenv(1, "cheev", "VU", n, -1, -1, -1);

    int lwork = mat.nc() * (nb+1);
    int info = 0;
    T w[n], wpow[n];
    T rwork[3*n-2];
    std::complex<T> work[lwork];
    if (mat.is_row_major())
        LapackConnector::heev(jobz, uplo, n, mat.ptr(), n, w, work, lwork, rwork, info);
    else
        LapackConnector::heev_f(jobz, uplo, n, mat.ptr(), n, w, work, lwork, rwork, info);
    bool is_int_power = fabs(power - int(power)) < 1e-8;
    // debug
    // for ( int i = 0; i != n; i++ )
    // {
    //     cout << w[i] << " ";
    // }
    // cout << endl;
    // end debug
    for ( int i = 0; i != n; i++ )
    {
        if (w[i] < 0 && w[i] > threshold && !is_int_power)
            printf("Warning! kept negative eigenvalue with non-integer power: # %d ev = %f , pow = %f\n", i, w[i], power);
        if (fabs(w[i]) < 1e-10 && power < 0)
            printf("Warning! nearly-zero eigenvalue with negative power: # %d ev = %f , pow = %f\n", i, w[i], power);
        if (w[i] < threshold)
        {
            wpow[i] = 0;
            if (filter_original)
                w[i] = 0;
        }
        else
            wpow[i] = w[i];
        wpow[i] = pow(wpow[i], power);
    }
    auto evconj = transpose(mat, true);
    auto temp = evconj.copy();
    for (int i = 0; i != n; i++)
        temp.scale_row(i, wpow[i]);
    auto pmat = mat * temp;
    // recover the original matrix here
    temp = mat.copy();
    for (int i = 0; i != n; i++)
        evconj.scale_row(i, w[i]);
    matmul(temp, evconj, mat);
    return pmat;
}

template <typename T>
bool operator==(const matrix_m<T> &m1, const matrix_m<T> &m2)
{
    return m1.dataobj == m2.dataobj;
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

    LapackConnector::heev_f('V', 'U', n, mat.ptr(), mat.nr(), w, work, lwork, rwork, info);
    delete [] rwork, w, work;
    return matrix_m<std::complex<T>>(n, n, mat.ptr(), MAJOR::COL, major);
}

//! generate a random Hermitian matrix with selected eigenvalues.
//! when the number of provided eigenvalues is smaller than the dimension,
//! padding_zero will be used as the missing eigenvalues
template <typename T>
inline matrix_m<std::complex<T>> random_he_selected_ev(int n, const vector<T> &evs, MAJOR major = MAJOR::ROW, const T& padding_zero = 1e-14)
{
    const int nvec = evs.size();
    if (nvec == 0) return matrix_m<std::complex<T>>{0, 0, major};

    auto eigvec = random_unitary<T>(n, major);
    auto eigvec_H = eigvec.get_transpose(true);
    for (int ie = 0; ie != n; ie++)
        eigvec_H.scale_row(ie, ie < nvec? evs[ie]: padding_zero);
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
        if (fabs(mat.dataobj[i]) > threshold)
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

template <typename T, typename Treal = typename to_real<T>::type>
void print_whole_matrix(const char *desc, const matrix_m<T> &mat)
{
    int nr = mat.nr();
    int nc = mat.nc();
    printf("\n %s\n", desc);
    printf("nr = %d, nc = %d\n", nr, nc);
    if(is_complex<T>())
    {
        for (int i = 0; i < nr; i++)
        {
            for (int j = 0; j < nc; j++)
                printf("%10.6f,%9.6f ", mat.c[i * nc + j].real(), mat.c[i * nc + j].imag());
            printf("\n");
        }
    }
    else
    {
        for (int i = 0; i < nr; i++)
        {
            for (int j = 0; j < nc; j++)
                printf("%10.6f", mat.c[i * nc + j].real());
            printf("\n");
        }
    }
}

template <typename T, typename Treal = typename to_real<T>::type>
void print_matrix_mm_file(const matrix_m<T> &mat, const char *fn, Treal threshold = 1e-15, bool row_first = true)
{
    ofstream fs;
    fs.open(fn);
    print_matrix_mm(mat, fs, threshold, row_first);
}
