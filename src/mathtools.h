/*!
 * @file mathtools.h
 * @brief Basic mathematic tools for common tasks
 * @note Most functions are adapted from Gauss_Quadrature.cpp
 */
#pragma once
#include <cstddef>

//! Get the Gauss-Chebyshev quadrature based on Chebyshev points of the first type (kind)
/*!
 * @param N Number of grid points
 * @param nodes Locations of grid points
 * @param weights Weights of grid points
 */
void GaussChebyshevI_unit(size_t N, double *nodes, double *weights);

//! Get the Gauss-Chebyshev quadrature based on Chebyshev points of the second type (kind)
/*!
 * @param N Number of grid points
 * @param nodes Locations of grid points
 * @param weights Weights of grid points
 */
void GaussChebyshevII_unit(size_t N, double *nodes, double *weights);

//! Compute N-point Gauss-Legendre grids in [-1, 1];
/*!
 * @param N Number of grid points
 * @param nodes Locations of grid points
 * @param weights Weights of grid points
 */
void GaussLegendre_unit(size_t N, double *nodes, double *weights);

//! Compute Legendre polynomial of order N at x in [-1, 1];
/*!
 * @param N  Order of polynomial
 * @param x  Grid point to evaluate the polynomial
 * @param y  Value of the polynomial at x
 * @param dy value of the derivative of polynomial at x
 */
void LegendreP_order_N(int N, double x, double &y, double &dy);

//! Transform Gauss grids from [-1, 1] to [x0, infty]
/*!
 * @param x0 Left end of the targeted region
 * @param N Number of grid points
 * @param nodes Locations of grid points
 * @param weights Weights of grid points
 */
void transform_GaussQuad_unit2x0inf(double x0, size_t N, double *nodes, double *weights);

//! Transform Gauss grids from [-1, 1] to [x0, x1]
/*!
 * @param x0 Left end of the targeted region
 * @param x1 Right end of the targeted region
 * @param N Number of grid points
 * @param nodes Locations of grid points
 * @param weights Weights of grid points
 */
void transform_GaussQuad_unit2x0x1(double x0, double x1, size_t N, double *nodes, double *weights);

//! Transform Gauss grids from [-1, 1] to [-infty, infty]
/*!
 * @param N Number of grid points
 * @param nodes Locations of grid points
 * @param weights Weights of grid points
 */
void transform_GaussQuad_unit2inf(size_t N, double *nodes, double *weights);

//! Transform Gauss grids from [-1, 1] to [-infty, x0]
/*!
 * @param x0 right end of the targeted region
 * @param N Number of grid points
 * @param nodes Locations of grid points
 * @param weights Weights of grid points
 */
void transform_GaussQuad_unit2minfx0(double x0, size_t N, double *nodes, double *weights);
