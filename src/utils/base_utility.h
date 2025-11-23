#pragma once
#include <algorithm>
#include <vector>
#include <complex>

constexpr static const double DOUBLE_EQUAL_THRES = 1e-10;

typedef std::complex<double> cplxdb;

// check if passed template type argument is a complex value
// from https://stackoverflow.com/questions/41438493/how-to-identifying-whether-a-template-argument-is-stdcomplex
template <typename T>
struct is_complex_t: public std::false_type {};

template <typename T>
struct is_complex_t<std::complex<T>>: public std::true_type {};

template <typename T>
bool is_complex() { return is_complex_t<T>::value; }

template <typename T>
struct to_real { using type = T; };
template <typename T>
struct to_real<std::complex<T>> { using type = T; };

template <typename T>
struct to_cplx { using type = std::complex<T>; };
template <typename T>
struct to_cplx<std::complex<T>> { using type = T; };

template <typename T>
inline void join_re_im(T &pt, const T &re, const T &im) { pt = re; }
template <typename T>
inline void join_re_im(std::complex<T> &pt, const T &re, const T &im) { pt = std::complex<T>{re, im}; }

template <typename T>
inline T get_real(const T& v) { return v; }
template <typename T>
inline T get_real(const std::complex<T>& v) { return v.real(); }
template <typename T>
inline T get_imag(const T& v) { return T(0); }
template <typename T>
inline T get_imag(const std::complex<T>& v) { return v.imag(); }
template <typename T>
inline T get_conj(const T& v) { return v; }
template <typename T>
inline std::complex<T> get_conj(const std::complex<T>& v) { return std::conj(v); }

template <typename T, typename T1, typename T2, typename TR = typename to_real<T>::type>
TR norm(const T v[], const T1 &n, const T2 &power)
{
    T norm = 0;
    for (int i = 0; i < n; i++)
        norm += std::pow(std::fabs(v[i]), power);
    return std::pow(norm, 1./power);
}

template <typename TI>
std::vector<TI> flatten_2d_indices(const std::vector<std::pair<TI, TI>> &indices_2d,
                                   const TI &rows, const TI &cols,
                                   bool row_major)
{
    std::vector<TI> indices_1d;
    indices_1d.reserve(indices_2d.size());

    if (row_major)
    {
        for (const auto &index_2d: indices_2d)
        {
            indices_1d.push_back(index_2d.first * cols + index_2d.second);
        }
    }
    else // column major
    {
        for (const auto &index_2d: indices_2d)
        {
            indices_1d.push_back(index_2d.first + index_2d.second * rows);
        }
    }
    return indices_1d;
}

template <typename T>
constexpr T ceil_div(T num, T den)
{
    return (num + den - 1) / den;  // assumes num >= 0, den > 0
}

inline std::size_t as_size(int x) noexcept
{
    return static_cast<std::size_t>(x);
}

template <typename T>
inline int as_int(T x) noexcept
{
    return static_cast<int>(x);
}

template <class Pair>
struct FastLess
{
    bool row_fast;  // true: (second, first); false: (first, second)
    bool operator()(const Pair &a, const Pair &b) const noexcept
    {
        if (row_fast)
        {
            return std::tie(a.second, a.first) < std::tie(b.second, b.first);
        }
        else
        {
            return std::tie(a.first, a.second) < std::tie(b.first, b.second);
        }
    }
};
