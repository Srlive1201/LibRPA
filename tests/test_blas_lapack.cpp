#include "blas_connector.h"
#include "lapack_connector.h"
#include <cassert>
#include <cstdlib>
#include <ctime>

template<typename T>
bool is_mat_A_equal_B(const int m, const int n, const T *A, const T* B, bool transpose, bool print,
                      const char stra[]  = "A", const char strb[]  = "B")
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
            if ( fabs(a-b) > 1e-15 )
            {
                ndiff++;
                if (!print) return false;
            }
        }
    return ndiff == 0;
}

template <typename T>
bool is_matmul_AB_equal_C(const int m, const int n, const int k,
                          const T *A, const T *B, const T *C, bool transpose_C, bool print)
{
    T AB[m][n];
    for ( int im = 0; im != m; im++)
        for ( int in = 0; in != n; in++)
        {
            AB[im][in] = 0;
            for ( int ik = 0; ik != k; ik++)
                AB[im][in] += A[im*k+ik] * B[ik*n+in];
        }
    return is_mat_A_equal_B(m, n, *AB, C, transpose_C, print, "AB", "C");
}

void test_blas_float() {}

void test_blas_double()
{
    const int RDIM = 6;
    const int IDIM = 2;
    const int JDIM = 3;
    assert( IDIM * JDIM == RDIM );
    const int CDIM = 2;
    double v1[RDIM] { -1, 2, 3, 1, 4, 6};
    double v2[RDIM] { 4, 0, 1, -1, 2, 3};
    double mat[RDIM][CDIM] {
        {1, -1}, {-1, 2}, {3, 0}, {0, -1}, {2, 3}, {0, 1}
    };
    int incr = 1;
    // vector dot
    assert(ddot_(&RDIM, v1, &incr, v2, &incr) == 24);
    assert(LapackConnector::dot(RDIM, v1, incr, v2, incr) == 24);

    // matrix-vector
    double y[CDIM] {};
    const char transn = 'N';
    const double alpha = 1.0, beta = 0.0;
    dgemv_(&transn, &CDIM, &RDIM, &alpha, *mat, &CDIM, v1, &incr, &beta, y, &incr);
    assert( y[0] == 14.0 && y[1] == 22.0 );
    LapackConnector::gemv('T', RDIM, CDIM, 1.0, *mat, CDIM, v1, 1, 0.5, y, 1);
    assert( y[0] == 21.0 && y[1] == 33.0 );

    srand(time(0));

    // matrix-matrix
    double mat2[CDIM][RDIM] = {};
    for ( int ic = 0; ic != CDIM; ic++ )
        for ( int ir = 0; ir != RDIM; ir++ )
            mat2[ic][ir] = double(rand()) / RAND_MAX;

    double mat_or[RDIM][RDIM] = {};
    double mat_oc[CDIM][CDIM] = {};
    const char transt = 'T';
    // mat[RDIM][CDIM] * mat2[CDIM][RIM] = m[RDIM][RDIM]
    dgemm_(&transn, &transn, &RDIM, &RDIM, &CDIM, &alpha, *mat2, &RDIM, *mat, &CDIM, &beta, *mat_or, &RDIM);
    assert( is_matmul_AB_equal_C(RDIM, RDIM, CDIM, *mat, *mat2, *mat_or, false, false));
    LapackConnector::gemm(transn, transn, RDIM, RDIM, CDIM, 1.0, *mat, CDIM, *mat2, RDIM, 0.0, *mat_or, RDIM);
    assert( is_matmul_AB_equal_C(RDIM, RDIM, CDIM, *mat, *mat2, *mat_or, false, false));

    // Check Cs like per-mu gemm
    const size_t n_mu = 3, n_i = 4, n_k = 12, n_j = 4;
    double Cs[n_mu][n_i*n_k];
    double G[n_k][n_j];
    double res[n_mu][n_i*n_j] = {};
    double ref[n_mu][n_i*n_j] = {};
    // initialize
    for ( int k = 0; k != n_k; k++ )
    {
        for ( int mu = 0; mu != n_mu; mu++ )
            for ( int i = 0; i != n_i; i++ )
            {
                Cs[mu][i*n_k+k] = double(rand()) / RAND_MAX;
            }
        for ( int j = 0; j != n_j; j++ )
            G[k][j] = double(rand()) / RAND_MAX;
    }
    // calculate reference in direct way
    for ( int mu = 0; mu != n_mu; mu++ )
        for ( int i = 0; i != n_i; i++ )
            for ( int k = 0; k != n_k; k++ )
                for ( int j = 0; j != n_j; j++ )
                {
                    ref[mu][i*n_j+j] += Cs[mu][i*n_k+k] * G[k][j];
                }
    for ( int mu = 0; mu != n_mu; mu++ )
        LapackConnector::gemm('N', 'N', n_i, n_j, n_k, 1.0,
                              Cs[mu], n_k,
                              *G, n_j,
                              1.0, res[mu], n_j);
    assert(is_mat_A_equal_B(n_mu, n_i*n_j, *res, *ref, false, false));
}

void test_blas_complex() {}

int main (int argc, char *argv[])
{
    test_blas_float();
    test_blas_double();
    test_blas_complex();
    return 0;
}
