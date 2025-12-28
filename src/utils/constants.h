#include <complex>
#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace librpa_int {

// Scientific constants
extern const double PINF;
extern const double NINF;
extern const double PI;
extern const double TWO_PI;
// Non-const to behave as customizable constant
extern double HA2EV;
extern double RY2EV;
extern double BOHR2ANG;
extern double ANG2BOHR;

// Complex numbers
constexpr std::complex<double> C_ONE{1.0, 0.0};
constexpr std::complex<double> C_ZERO{0.0, 0.0};

constexpr size_t kbytes = 1024;
constexpr size_t mbytes = kbytes * 1024;
constexpr size_t gbytes = mbytes * 1024;

//! Set constants to those in FHI-aims.
//! Internally everything is computed in atomic unit.
void set_aims_constants();

}
#endif
