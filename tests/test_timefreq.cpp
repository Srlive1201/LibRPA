#include "timefreq.h"
#include <stdexcept>

int main (int argc, char *argv[])
{
    TFGrids tfg(16);
    // Check minimax grids
    // emin and emax (in Hartree) from diamond PBE k2c
    double emin = 0.173388;
    double emax = 23.7738;
    tfg.generate_minimax(emin, emax);
    if (tfg.get_grid_type() != TFGrids::GRID_TYPES::Minimax)
        throw logic_error("internal type should be minimax grid");
    return 0;
}
