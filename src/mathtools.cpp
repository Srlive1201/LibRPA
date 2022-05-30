/*!
 * @file mathtools.cpp
 * @note Many functions are directly copied from Gauss_Quadrature.cpp
 */
#include "mathtools.h"
#include "constants.h"
#include <cmath>
#include <stdexcept>

void GaussChebyshevI_unit(size_t N, double *nodes, double *weights)
{
    for ( int i = 0; i != (N+1)/2; i++)
    {
        nodes[i] = cos(PI*(i+0.5)/N);
        weights[i] = sqrt(1.0 - nodes[i] * nodes[i] ) * PI / N;
        nodes[N-1-i] = -nodes[i];
        weights[N-1-i] = -weights[i];
    }
}

void GaussChebyshevII_unit(size_t N, double *nodes, double *weights)
{
    size_t Np1 = N + 1;
    for ( int i = 0; i != Np1/2; i++)
    {
        nodes[i] = cos(PI*(i+1)/Np1);
        weights[i] = 1.0 / sqrt(1.0 - nodes[i] * nodes[i] ) * PI / Np1
            * pow(sin(PI*(i+1)/Np1), 2);
        nodes[N-1-i] = -nodes[i];
        weights[N-1-i] = -weights[i];
    }
}

void GaussLegendre_unit(size_t N, double *nodes, double *weights)
{
    const double thres = 1e-16;
    const size_t maxiter = 1000;
    double x, y, dy, ratio;
    size_t halfNp1;
    halfNp1 = (N+1)/2;
    for ( int i = 0; i != halfNp1; i++ )
    {
        int ite = 0;
        x = cos( PI *( i+1 -0.25 ) / ( N +0.5 ) );
        while ( ++ite != maxiter )
        {
            LegendreP_order_N(N, x, y, dy);
            ratio = y / dy;
            x -= ratio;
            if ( fabs(ratio) < thres )
                break;
        }
        if ( ite == maxiter) throw std::logic_error("Legendre Newton-Raphson not converged");
        // apply symmetry
        nodes[N-1-i] = -nodes[i];
        weights[N-1-i] = -weights[i];
    }
}

void transform_GaussQuad_unit2x0inf(double x0, size_t N, double *nodes, double *weights)
{
    for ( int i = 0; i != N; i++ )
    {
        weights[i] = weights[i] /( ( 1.0- nodes[i] )*( 1.0- nodes[i] ) );
        nodes[i] = 0.5* ( nodes[i] +1.0 ) / ( 1.0- nodes[i] ) + x0;
    }
}

void transform_GaussQuad_unit2x0x1(double x0, double x1, size_t N, double *nodes, double *weights)
{
    double midp = (x0+x1)/2;
    double hwid = (x1-x0)/2;

    for( int i = 0; i != N; i++ )
    {
        nodes[i] = hwid * nodes[i] +midp;
        weights[i] = hwid * weights[i];
    }
}

void LegendreP_order_N(int N, double x, double &y, double &dy)
{
    double p0, p, dp, rtmp, xtmp;
    int n;

    p0 = 1.0;    p = x;
    xtmp = 1.0/ ( x*x -1.0 );
    /* Recurrence relation:  for n = 1, 2, ..., N,
     *    nL_n(x)  = (2n-1)xL_{n-1}(x) - (n-1)L_{n-2}(x)
     *    (x^2-1)L'_n(x) = nxL_n(x) - nL_{n-1}(x)
     */
    for( n = 1; n != N; ++n )
    {
        rtmp = p0;    p0 = p;
        p = ( ( 2*n +1 ) *x *p0 - n *rtmp ) /(n+1);
        dp = (n+1) * xtmp * ( x*p - p0 );
    }
    y = p; dy = dp;
}

