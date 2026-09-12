#include "../read_data.h"

#include "../../src/utils/constants.h"

#include <cassert>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <string>

int main()
{
    const auto test_dir =
        std::filesystem::temp_directory_path() / "librpa_test_velocity_reader";
    std::filesystem::remove_all(test_dir);
    std::filesystem::create_directories(test_dir);

    const auto velocity_file = test_dir / "velocity_matrix.txt";
    std::ofstream out(velocity_file);
    out << "1\n1\n1\n1\n";
    out << "1 1 1\n2.0 3.0\n";
    out << "2 1 1\n4.0 5.0\n";
    out << "3 1 1\n6.0 7.0\n";
    out.close();

    // A second prefix match must not shadow the canonical ABACUS text file.
    std::ofstream(test_dir / "velocity_matrix_h.txt") << "not a LibRPA velocity file\n";

    MeanField mf(1, 1, 1, 1);
    velocity_matrix_t velocity;
    read_velocity_abacus(mf, test_dir.string(), "velocity_matrix", velocity);

    const double scale = librpa_int::ANG2BOHR / librpa_int::HA2EV;
    const auto expected = scale * std::complex<double>(2.0, 3.0);
    assert(std::abs(velocity.at(0).at(0).at(0)(0, 0) - expected) < 1e-12);

    std::filesystem::remove_all(test_dir);
}
