#pragma once
#include <cstddef>

void GaussChebyshevI_unit(size_t N, double *nodes, double *weights);
void GaussChebyshevII_unit(size_t N, double *nodes, double *weights);
void GaussLegendre_unit(size_t N, double *nodes, double *weights);

void transform_GaussQuad_unit2x0inf(double x0, size_t N, double *nodes, double *weights);
void transform_GaussQuad_unit2x0x1(double x0, double x1, size_t N, double *nodes, double *weights);

//! Compute Legendre polynomial of order N at x in [-1, 1];
void LegendreP_order_N(int N, double x, double &y, double &dy);
