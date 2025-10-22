#include <complex>
#include <iostream>
#include <vector>
#include <map>
#include <utility>


#ifdef LIBRPA_DEBUG
constexpr bool print_mat_equal = true;
#else
constexpr bool print_mat_equal = false;
#endif

template <typename T>
inline bool fequal(const T &a, const T &b, const T &thres = 1e-14)
{
    return std::abs(a - b) < std::abs(thres);
}

template <typename T>
bool fequal_array(int n, const std::complex<T> *a,
                  const std::complex<T> *b, bool print, const std::complex<T> thres = 1e-14)
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

template <typename T>
bool equal_array(int n, const T *a, const T *b, bool print = false)
{
    int ndiff = 0;
    for ( int i = 0; i != n; i++ )
    {
        if (print)
        {
            std::cout << i << " : a " << *(a+i) << " , " << *(b+i) << std::endl;
        }
        if ( *(a+i) != *(b+i) )
        {
            ndiff++;
            if (!print)
                return false;
        }
    }
    return ndiff == 0;
}

template <typename T>
bool equal_vector(const std::vector<T> &v1, const std::vector<T> &v2, bool print = false)
{
    if (v1.size() != v2.size()) return false;
    return equal_array(v1.size(), v1.data(), v2.data(), print);
}

template <typename T1, typename T2>
bool fequal_pair(const std::pair<T1, T2> &p1, const std::pair<T1, T2> &p2,
                 const T1 &thres1 = 1e-14, const T1 &thres2 = 1e-14)
{
    return fequal(p1.first, p2.first, thres1) && fequal(p1.second, p2.second, thres2);
}

template <typename T1, typename T2>
bool equal_pair(const std::pair<T1, T2> &p1, const std::pair<T1, T2> &p2)
{
    return p1.first == p2.first && p1.second == p2.second;
}

template <typename Tk, typename Tv1, typename Tv2>
bool equal_map_vector_pv(const std::map<Tk, std::vector<std::pair<Tv1, std::vector<Tv2>>>> &mvpv1,
                         const std::map<Tk, std::vector<std::pair<Tv1, std::vector<Tv2>>>> &mvpv2,
                         bool print = false)
{
    if (mvpv1.size() != mvpv2.size())
    {
        if (print) std::cout << "Map size differ: "
                             << mvpv1.size() << " != " << mvpv2.size() << std::endl;
        return false;
    }

    bool flag = true;
    for (const auto &kv1: mvpv1)
    {
        const auto &key = kv1.first;
        const auto &v1 = kv1.second;
        if (mvpv2.count(kv1.first) == 0)
        {
            if (print) std::cout << "Missing key in map2: " << key << std::endl;
            return false;
        }
        const auto &v2 = mvpv2.at(key);
        if (v1.size() != v2.size()) return false;
        if (print) std::cout << "Key " << key << std::endl;
        for (size_t i = 0; i < v1.size(); i++)
        {
            const auto &v1p = v1[i].first;
            const auto &v2p = v2[i].first;
            const auto &v1v = v1[i].second;
            const auto &v2v = v2[i].second;
            if (!equal_pair(v1p, v2p))
            {
                if (print) std::cout << "Pair identifier differ at vector index " << i << " : " << v1p << " != " << v2p << std::endl;
                return false;
            }
            if (v1v.size() != v2v.size())
            {
                if (print) std::cout << "Vector for the pair differ in size: " << v1v.size() << " != " << v2v.size() << std::endl;
                return false;
            }
            if (print) std::cout << "Vector for the pair: " << v1v.size() << " != " << v2v.size() << std::endl;
            flag = flag && equal_array(v1v.size(), v1v.data(), v2v.data(), print);
        }
    }
    return flag;
}

template <typename Tk, typename Tv>
bool equal_map_vector(const std::map<Tk, std::vector<Tv>> &mv1,
                      const std::map<Tk, std::vector<Tv>> &mv2,
                      bool print = false)
{
    if (mv1.size() != mv2.size())
    {
        if (print) std::cout << "Map size differ: "
                             << mv1.size() << " != " << mv2.size() << std::endl;
        return false;
    }
    bool flag = true;
    for (const auto &kv1: mv1)
    {
        const auto &key = kv1.first;
        const auto &v1 = kv1.second;
        if (mv2.count(kv1.first) == 0)
        {
            if (print) std::cout << "Missing key in map2: " << key << std::endl;
            return false;
        }
        const auto &v2 = mv2.at(key);
        if (v1.size() != v2.size()) return false;
        if (print) std::cout << "Key " << key << std::endl;
        flag = flag && equal_array(v1.size(), v1.data(), v2.data(), print);
    }
    return flag;
}

template <typename T>
bool fequal_array(int n, const T *a, const T *b, bool print = false, T thres = 1e-14)
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

template <typename T>
inline std::complex<T> cmplx(T re, T im)
{
    return std::complex<T>(re, im);
}

template<typename T>
bool is_mat_A_equal_B(const int m, const int n, const T *A, const T* B, bool transpose, bool print,
                      const T &thres = 1e-14, const char stra[] = "A", const char strb[] = "B")
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

template <typename T>
bool is_mat_A_equal_B(const int m, const int n, const std::complex<T> *A,
                      const std::complex<T> *B, bool transpose, bool print,
                      const std::complex<T> &thres, const char stra[] = "A", const char strb[] = "B")
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
                          const T *A, const T *B, const T *C, bool transpose_C = false, bool print = false, T thres = 1e-14)
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
