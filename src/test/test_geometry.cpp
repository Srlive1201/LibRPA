#include <cassert>
#include <cmath>

#include "../core/geometry.h"
#include "../utils/constants.h"
#include "testutils.h"

static void check_cubic()
{
    using namespace librpa_int;

    Atoms atoms;

    atoms.set({}, {});
    atoms.set({0, 1}, {});
    assert(atoms.types[0] == 0);
    assert(atoms.types[1] == 1);
    atoms.set({},
              {{0.3, 3.4, 5.1}, {0.0, 0.3, -0.2}});
    assert(fequal(atoms.coords[0][0], 0.3));
    assert(fequal(atoms.coords[0][1], 3.4));
    assert(fequal(atoms.coords[0][2], 5.1));
    assert(fequal(atoms.coords[1][0], 0.0));
    assert(fequal(atoms.coords[1][1], 0.3));
    assert(fequal(atoms.coords[1][2], -0.2));

    Matrix3 latt{0.0, 0.5, 0.5, 0.5, 0.0, 0.5, 0.5, 0.5, 0.0};
    // Reset
    atoms.set({0, 1, 0},
              {{0.00, 0.00, 0.00},
               {0.25, 0.25, 0.25},
               {0.50, 0.50, 0.00}}, latt);
    assert(atoms.types[0] == 0);
    assert(atoms.types[1] == 1);
    assert(atoms.types[2] == 0);
    assert(fequal(atoms.coords_frac[0][0], 0.0));
    assert(fequal(atoms.coords_frac[0][1], 0.0));
    assert(fequal(atoms.coords_frac[0][2], 0.0));
    // std::cout << cubic.coords_frac[1][0] << " "
    //           << cubic.coords_frac[1][1] << " "
    //           << cubic.coords_frac[1][2] << std::endl;
    assert(fequal(atoms.coords_frac[1][0], 0.25));
    assert(fequal(atoms.coords_frac[1][1], 0.25));
    assert(fequal(atoms.coords_frac[1][2], 0.25));
    // std::cout << cubic.coords_frac[2][0] << " "
    //           << cubic.coords_frac[2][1] << " "
    //           << cubic.coords_frac[2][2] << std::endl;
    assert(fequal(atoms.coords_frac[2][0], 0.0));
    assert(fequal(atoms.coords_frac[2][1], 0.0));
    assert(fequal(atoms.coords_frac[2][2], 1.0));
}

static void check_hexagonal()
{
    using namespace librpa_int;

    Atoms atoms;
    Matrix3 latt{3.60, 0.0, 0.0, -1.80, 1.80 * std::sqrt(3.0), 0.0, 0.0, 0.0, 20.0};
    // Reset
    atoms.set({0, 1},
              {{0.00, 0.00, 10.00},
               {0.00, 3.6 / std::sqrt(3), 10.00}}, latt);
    assert(fequal(atoms.coords_frac[0][0], 0.0));
    assert(fequal(atoms.coords_frac[0][1], 0.0));
    assert(fequal(atoms.coords_frac[0][2], 0.5));
    // std::cout << cubic.coords_frac[1][0] << " "
    //           << cubic.coords_frac[1][1] << " "
    //           << cubic.coords_frac[1][2] << std::endl;
    assert(fequal(atoms.coords_frac[1][0], 1.0 / 3.0));
    assert(fequal(atoms.coords_frac[1][1], 2.0 / 3.0));
    assert(fequal(atoms.coords_frac[1][2], 0.5));
}

int main (int argc, char *argv[])
{
    check_cubic();
    check_hexagonal();
    return 0;
}
