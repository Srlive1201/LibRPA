#include "constants.h"
#include <cmath>

namespace librpa_int {

const double PINF = INFINITY;
const double NINF = - INFINITY;
const double PI = M_PI;
const double TWO_PI = 2.0 * PI;

// NIST CODATA 2022
// https://physics.nist.gov/cgi-bin/cuu/Value?eqhrev
double HA2EV = 27.211386245988;
double RY2EV = HA2EV * 0.5;

// NIST CODATA 2022
// https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0
double BOHR2ANG = 0.529177210544;
double ANG2BOHR = 1.0 / BOHR2ANG;

void set_aims_constants()
{
    // FHI-aims uses the following value as default, e.g. "FHIaims_pre2024" as of 2025.
    // To compare the driver results in eV with FHI-aims, the following should be used.
    HA2EV = 27.2113845e0;
    RY2EV = HA2EV * 0.5;
    BOHR2ANG = 0.5291772108;
    ANG2BOHR = 1.0 / BOHR2ANG;
}

}
