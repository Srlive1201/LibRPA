#pragma once
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

