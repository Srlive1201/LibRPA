#ifdef NDEBUG
#undef NDEBUG
#endif

#include "../qsgw/band_output.h"

#include "../utils/constants.h"

#include <cassert>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using librpa_int::HA2EV;
using librpa_int::MeanField;
using librpa_int::Matz;
using librpa_int::Vector3_Order;
using librpa_int::qsgw::SpinKMatrixMap;
using librpa_int::qsgw::write_qsgw_band_spin_tables;

namespace
{

template <typename Function>
void assert_throws(Function&& function)
{
    bool threw = false;
    try
    {
        function();
    }
    catch (const std::exception&)
    {
        threw = true;
    }
    assert(threw);
}

Matz diagonal(const double first, const double second)
{
    Matz result(2, 2);
    result(0, 0) = first;
    result(1, 1) = second;
    return result;
}

void test_legacy_band_columns_and_energy_components()
{
    MeanField reference(1, 2, 2, 2, 1);
    MeanField live = reference;
    for (int kpoint = 0; kpoint < 2; ++kpoint)
    {
        reference.get_weight()[0](kpoint, 0) = 1.0;
        reference.get_weight()[0](kpoint, 1) = 0.0;
        live.get_eigenvals()[0](kpoint, 0) = 0.1 + 0.1 * kpoint;
        live.get_eigenvals()[0](kpoint, 1) = 1.0 + 0.1 * kpoint;
    }

    SpinKMatrixMap h_ks;
    SpinKMatrixMap vxc;
    SpinKMatrixMap exx;
    h_ks[0][0] = diagonal(0.0, 0.8);
    h_ks[0][1] = diagonal(0.1, 0.9);
    vxc[0][0] = diagonal(0.2, 0.3);
    vxc[0][1] = diagonal(0.2, 0.3);
    exx[0][0] = diagonal(-0.1, -0.05);
    exx[0][1] = diagonal(-0.1, -0.05);
    exx[0][1](0, 0) += librpa_int::cplxdb(0.0, 0.25);
    const std::vector<Vector3_Order<double>> kpoints{
        {0.0, 0.0, 0.0}, {0.5, 0.0, 0.0}};

    std::ostringstream ks;
    std::ostringstream hf;
    std::ostringstream qsgw;
    write_qsgw_band_spin_tables(
        ks, hf, qsgw, live, reference, kpoints, h_ks, vxc, exx,
        0, 2.0);

    int index = 0;
    double kx = 0.0;
    double ky = 0.0;
    double kz = 0.0;
    double occ0 = 0.0;
    double energy0 = 0.0;
    double occ1 = 0.0;
    double energy1 = 0.0;
    std::istringstream ks_input(ks.str());
    ks_input >> index >> kx >> ky >> kz
             >> occ0 >> energy0 >> occ1 >> energy1;
    assert(index == 1);
    assert(std::abs(occ0 - 2.0) < 1.0e-12);
    assert(std::abs(occ1) < 1.0e-12);
    assert(std::abs(energy0) < 1.0e-12);
    assert(std::abs(energy1 - 0.8 * HA2EV) < 1.0e-4);

    std::istringstream hf_input(hf.str());
    hf_input >> index >> kx >> ky >> kz
             >> occ0 >> energy0 >> occ1 >> energy1;
    assert(std::abs(energy0 - (-0.3 * HA2EV)) < 1.0e-4);
    assert(std::abs(energy1 - (0.45 * HA2EV)) < 1.0e-4);

    std::istringstream qsgw_input(qsgw.str());
    qsgw_input >> index >> kx >> ky >> kz
               >> occ0 >> energy0 >> occ1 >> energy1;
    assert(std::abs(occ0 - 2.0) < 1.0e-12);
    assert(std::abs(occ1) < 1.0e-12);
    assert(std::abs(energy0 - 0.1 * HA2EV) < 1.0e-4);
    assert(std::abs(energy1 - 1.0 * HA2EV) < 1.0e-4);

    ks_input >> index >> kx >> ky >> kz
             >> occ0 >> energy0 >> occ1 >> energy1;
    assert(index == 2);
    assert(std::abs(kx - 0.5) < 1.0e-12);
    assert(std::abs(ky) < 1.0e-12);
    assert(std::abs(kz) < 1.0e-12);
    assert(std::abs(occ0 - 2.0) < 1.0e-12);
    assert(std::abs(occ1) < 1.0e-12);
    assert(std::abs(energy0 - 0.1 * HA2EV) < 1.0e-4);
    assert(std::abs(energy1 - 0.9 * HA2EV) < 1.0e-4);
    ks_input >> std::ws;
    assert(ks_input.peek() == std::char_traits<char>::eof());

    hf_input >> index >> kx >> ky >> kz
             >> occ0 >> energy0 >> occ1 >> energy1;
    assert(index == 2);
    assert(std::abs(energy0 - (-0.2 * HA2EV)) < 1.0e-4);
    assert(std::abs(energy1 - (0.55 * HA2EV)) < 1.0e-4);
    hf_input >> std::ws;
    assert(hf_input.peek() == std::char_traits<char>::eof());

    qsgw_input >> index >> kx >> ky >> kz
               >> occ0 >> energy0 >> occ1 >> energy1;
    assert(index == 2);
    assert(std::abs(occ0 - 2.0) < 1.0e-12);
    assert(std::abs(occ1) < 1.0e-12);
    assert(std::abs(energy0 - 0.2 * HA2EV) < 1.0e-4);
    assert(std::abs(energy1 - 1.1 * HA2EV) < 1.0e-4);
    qsgw_input >> std::ws;
    assert(qsgw_input.peek() == std::char_traits<char>::eof());
}

void test_invalid_band_output_contracts_are_rejected()
{
    MeanField reference(1, 1, 1, 1, 1);
    MeanField live = reference;
    reference.get_weight()[0](0, 0) = 2.0;
    live.get_eigenvals()[0](0, 0) = 0.0;
    SpinKMatrixMap matrices;
    matrices[0][0] = Matz(1, 1);
    const std::vector<Vector3_Order<double>> kpoints{{0.0, 0.0, 0.0}};

    assert_throws([&] {
        std::ostringstream a;
        std::ostringstream b;
        std::ostringstream c;
        write_qsgw_band_spin_tables(
            a, b, c, live, reference, {}, matrices, matrices,
            matrices, 0, 0.0);
    });
    assert_throws([&] {
        std::ostringstream a;
        std::ostringstream b;
        std::ostringstream c;
        write_qsgw_band_spin_tables(
            a, b, c, live, reference, kpoints, matrices, matrices,
            matrices, 1, 0.0);
    });
    assert_throws([&] {
        std::ostringstream a;
        std::ostringstream b;
        std::ostringstream c;
        write_qsgw_band_spin_tables(
            a, b, c, live, reference, kpoints, matrices, matrices,
            matrices, 0, std::numeric_limits<double>::quiet_NaN());
    });

    SpinKMatrixMap nonfinite = matrices;
    nonfinite.at(0).at(0)(0, 0) =
        std::numeric_limits<double>::infinity();
    assert_throws([&] {
        std::ostringstream a;
        std::ostringstream b;
        std::ostringstream c;
        write_qsgw_band_spin_tables(
            a, b, c, live, reference, kpoints, nonfinite, matrices,
            matrices, 0, 0.0);
    });

    nonfinite.at(0).at(0)(0, 0) = librpa_int::cplxdb(
        0.0, std::numeric_limits<double>::infinity());
    assert_throws([&] {
        std::ostringstream a;
        std::ostringstream b;
        std::ostringstream c;
        write_qsgw_band_spin_tables(
            a, b, c, live, reference, kpoints, nonfinite, matrices,
            matrices, 0, 0.0);
    });
}

} // namespace

int main()
{
    test_legacy_band_columns_and_energy_components();
    test_invalid_band_output_contracts_are_rejected();
    return 0;
}
