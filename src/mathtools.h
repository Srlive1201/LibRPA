#pragma once
#include <cstddef>

//! Get the Gauss-Chebyshev quadrature based on Chebyshev points of the first type (kind)
void GaussChebyshevI_unit(size_t N, double *nodes, double *weights);
//! Get the Gauss-Chebyshev quadrature based on Chebyshev points of the second type (kind)
void GaussChebyshevII_unit(size_t N, double *nodes, double *weights);

//! Compute N-point Gauss-Legendre grids in [-1, 1];
void GaussLegendre_unit(size_t N, double *nodes, double *weights);

//! Compute Legendre polynomial of order N at x in [-1, 1];
void LegendreP_order_N(int N, double x, double &y, double &dy);

//! Transform Gauss grids from [-1, 1] to [x0, infty]
void transform_GaussQuad_unit2x0inf(double x0, size_t N, double *nodes, double *weights);

//! Transform Gauss grids from [-1, 1] to [x0, x1]
void transform_GaussQuad_unit2x0x1(double x0, double x1, size_t N, double *nodes, double *weights);

//! Transform Gauss grids from [-1, 1] to [-infty, infty]
void transform_GaussQuad_unit2inf(size_t N, double *nodes, double *weights);

//! Transform Gauss grids from [-1, 1] to [-infty, x0]
void transform_GaussQuad_unit2minfx0(double x0, size_t N, double *nodes, double *weights);
