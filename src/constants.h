#include <complex>
#ifndef CONSTANTS_H
#define CONSTANTS_H

extern const double PINF;
extern const double NINF;
extern const double PI;
extern const double TWO_PI;
extern const double HA2EV;
extern const double RY2EV;
extern const double ANG2BOHR;
extern const double BOHR2ANG;

constexpr std::complex<double> C_ONE{1.0, 0.0};
constexpr std::complex<double> C_ZERO{0.0, 0.0};
#endif
