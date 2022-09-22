#include <complex>
template <typename T>
inline bool fequal(const T &a, const T &b, double thres = 1e-14)
{
    return std::abs(a - b) < thres;
}

template <typename T>
bool fequal_array(int n, const T *a, const T *b, bool print = false, double thres = 1e-14)
{
    int ndiff = 0;
    for ( int i = 0; i != n; i++ )
    {
        if (print)
            printf("%d : a %f , b %f\n", i, *(a+i), *(b+i));
        if ( !fequal(*(a+i), *(b+i), thres) )
        {
            ndiff++;
            if (!print)
                return false;
        }
    }
    return ndiff == 0;
}

template <>
bool fequal_array<std::complex<double>>(int n, const std::complex<double> *a,
                                        const std::complex<double> *b, bool print, double thres)
{
    int ndiff = 0;
    for ( int i = 0; i != n; i++ )
    {
        if (print)
            printf("%d : a (%f,%f), b (%f,%f)\n", i, (*(a+i)).real(),
                                                     (*(a+i)).imag(),
                                                     (*(b+i)).real(),
                                                     (*(b+i)).imag());
        if ( !fequal(*(a+i), *(b+i), thres) )
        {
            ndiff++;
            if (!print)
                return false;
        }
    }
    return ndiff == 0;
}

inline std::complex<double> cmplx(double re, double im)
{
    return std::complex<double>(re, im);
}

template<typename T>
bool is_mat_A_equal_B(const int m, const int n, const T *A, const T* B, bool transpose, bool print,
                      double thres = 1e-14, const char stra[]  = "A", const char strb[]  = "B")
{
    int ndiff = 0;
    for ( int im = 0; im != m; im++)
        for ( int in = 0; in != n; in ++)
        {
            T a = A[im*n+in];
            T b;
            if (transpose)
                b = B[in*m+im];
            else
                b = B[im*n+in];
            if ( print )
            {
                printf("%d %d: %s = %f, %s = %f\n", im, in, stra, a, strb, b);
            }
            if ( !fequal(a, b, thres) )
            {
                ndiff++;
                if (!print) return false;
            }
        }
    return ndiff == 0;
}

template<>
bool is_mat_A_equal_B<std::complex<double>>(const int m, const int n, const std::complex<double> *A,
                                            const std::complex<double> *B, bool transpose, bool print,
                                            double thres, const char stra[], const char strb[])
{
    int ndiff = 0;
    for ( int im = 0; im != m; im++)
        for ( int in = 0; in != n; in ++)
        {
            std::complex<double> a = A[im*n+in];
            std::complex<double> b;
            if (transpose)
                b = B[in*m+im];
            else
                b = B[im*n+in];
            if ( print )
            {
                printf("%d %d: %s = (%f,%f), %s = (%f,%f)\n", im, in, stra, a.real(), a.imag(), strb, b.real(), b.imag());
            }
            if ( !fequal(a, b, thres) )
            {
                ndiff++;
                if (!print) return false;
            }
        }
    return ndiff == 0;
}

template <typename T>
bool is_matmul_AB_equal_C(const int m, const int n, const int k,
                          const T *A, const T *B, const T *C, bool transpose_C, bool print, double thres = 1e-14)
{
    T AB[m][n];
    for ( int im = 0; im != m; im++)
        for ( int in = 0; in != n; in++)
        {
            AB[im][in] = 0;
            for ( int ik = 0; ik != k; ik++)
                AB[im][in] += A[im*k+ik] * B[ik*n+in];
        }
    return is_mat_A_equal_B(m, n, *AB, C, transpose_C, print, thres, "AB", "C");
}

