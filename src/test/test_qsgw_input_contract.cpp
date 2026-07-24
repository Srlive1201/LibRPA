#ifdef NDEBUG
#undef NDEBUG
#endif

#include "../qsgw/input_contract.h"

#include <cassert>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using librpa_int::qsgw::BandUpdateMode;
using librpa_int::qsgw::HeadwingGridMode;
using librpa_int::qsgw::HeadwingUpdateMode;
using librpa_int::qsgw::QsgwInputContract;
using librpa_int::qsgw::QsgwProducer;
using librpa_int::qsgw::resolve_band_reference_paths;
using librpa_int::qsgw::resolve_same_grid_velocity_paths;
using librpa_int::qsgw::validate_band_reference_binding;
using librpa_int::qsgw::validate_qsgw_execution_modes;
using librpa_int::qsgw::validate_scf_input_binding;

namespace
{

constexpr const char* abc_sha256 =
    "ba7816bf8f01cfea414140de5dae2223"
    "b00361a396177a9cb410ff61f20015ad";

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

std::string supported_contract(const bool head, const bool band)
{
    std::ostringstream text;
    text << "# librpa-qsgw-input-contract-v1\n"
         << "producer\tfhi-aims\n"
         << "internal_energy_units\thartree\n"
         << "mf0_basis\tstate_coefficients_in_nao\n"
         << "mf0_gauge\tproducer_state\n"
         << "n_spins\t1\n"
         << "n_bands\t4\n"
         << "n_aos\t4\n"
         << "n_scf_kpoints\t2\n"
         << "n_headwing_kpoints\t" << (head ? 2 : 0) << "\n"
         << "n_band_kpoints\t" << (band ? 2 : 0) << "\n"
         << "headwing_grid\t" << (head ? "scf" : "disabled") << "\n"
         << "headwing_update\t"
         << (head ? "fixed_reference" : "none") << "\n"
         << "hartree_update\toff\n"
         << "band_update\t"
         << (band ? "fixed_basis_rotation" : "off") << "\n"
         << "role\tsha256\tfile\n"
         << "mf0_eigenvalues\t" << abc_sha256
         << "\tmf0_eigenvalues.dat\n"
         << "mf0_wavefunctions\t" << abc_sha256
         << "\tmf0_wavefunctions.dat\n"
         << "scf_kpoints\t" << abc_sha256
         << "\tscf_kpoints.dat\n"
         << "vxc_scf_manifest\t" << abc_sha256
         << "\tvxc_scf.manifest\n"
         << "reader_static\t" << abc_sha256
         << "\treader_static.dat\n";
    if (head)
    {
        text << "velocity_mf0\t" << abc_sha256
             << "\tmommat_ks_kpt_000001.dat\n"
             << "velocity_mf0\t" << abc_sha256
             << "\tmommat_ks_kpt_000002.dat\n";
    }
    if (band)
    {
        text << "band_kpoints\t" << abc_sha256
             << "\tband_kpath_info\n"
             << "vxc_band_manifest\t" << abc_sha256
             << "\tvxc_band.manifest\n";
        for (int kpoint = 1; kpoint <= 2; ++kpoint)
        {
            std::ostringstream index;
            index << std::setw(5) << std::setfill('0') << kpoint;
            text << "band_mf0_eigenvalues\t" << abc_sha256
                 << "\tband_KS_eigenvalue_k_" << index.str() << ".txt\n"
                 << "band_mf0_wavefunctions\t" << abc_sha256
                 << "\tband_KS_eigenvector_k_" << index.str() << ".txt\n";
        }
    }
    return text.str();
}

QsgwInputContract parse(const std::string& text)
{
    std::istringstream input(text);
    return QsgwInputContract::parse(input, "test-contract");
}

void test_supported_head_only_band_contract()
{
    const QsgwInputContract contract = parse(supported_contract(true, true));
    assert(contract.producer() == QsgwProducer::FhiAims);
    assert(contract.headwing_grid() == HeadwingGridMode::ScfGrid);
    assert(contract.headwing_update() ==
           HeadwingUpdateMode::FixedReference);
    assert(contract.band_update() == BandUpdateMode::FixedBasisRotation);
    assert(contract.n_spins() == 1);
    assert(contract.n_bands() == 4);
    assert(contract.n_aos() == 4);
    assert(contract.n_scf_kpoints() == 2);
    assert(contract.n_headwing_kpoints() == 2);
    assert(contract.n_band_kpoints() == 2);
    assert(contract.files("velocity_mf0").size() == 2);

    validate_qsgw_execution_modes(
        contract, HeadwingGridMode::ScfGrid, true);
    assert_throws([&] {
        validate_qsgw_execution_modes(
            contract, HeadwingGridMode::Disabled, true);
    });
    assert_throws([&] {
        validate_qsgw_execution_modes(
            contract, HeadwingGridMode::ScfGrid, false);
    });
}

void test_disabled_head_and_band_contract()
{
    const QsgwInputContract contract = parse(supported_contract(false, false));
    assert(contract.headwing_grid() == HeadwingGridMode::Disabled);
    assert(contract.headwing_update() == HeadwingUpdateMode::None);
    assert(contract.band_update() == BandUpdateMode::Off);
    validate_qsgw_execution_modes(
        contract, HeadwingGridMode::Disabled, false);
}

void test_unsupported_modes_are_rejected()
{
    std::string independent = supported_contract(true, false);
    const auto grid = independent.find("headwing_grid\tscf");
    independent.replace(grid, std::string("headwing_grid\tscf").size(),
                        "headwing_grid\tindependent_full");
    assert_throws([&] { (void)parse(independent); });

    std::string live_fourier = supported_contract(true, false);
    const auto update =
        live_fourier.find("headwing_update\tfixed_reference");
    live_fourier.replace(
        update,
        std::string("headwing_update\tfixed_reference").size(),
        "headwing_update\tlive_ao_fourier");
    assert_throws([&] { (void)parse(live_fourier); });

    std::string hartree = supported_contract(true, false);
    const auto mode = hartree.find("hartree_update\toff");
    hartree.replace(mode, std::string("hartree_update\toff").size(),
                    "hartree_update\tdelta_density");
    assert_throws([&] { (void)parse(hartree); });
}

void test_reader_bindings_match_exact_files()
{
    const std::filesystem::path root =
        std::filesystem::absolute("test_qsgw_contract_binding.tmp")
            .lexically_normal();
    const QsgwInputContract contract = parse(supported_contract(true, true));

    validate_scf_input_binding(
        contract, root.string(), root / "mf0_eigenvalues.dat",
        {root / "mf0_wavefunctions.dat"}, root / "scf_kpoints.dat",
        {root / "reader_static.dat"});
    assert_throws([&] {
        validate_scf_input_binding(
            contract, root.string(), root / "mf0_eigenvalues.dat",
            {root / "wrong_wavefunctions.dat"}, root / "scf_kpoints.dat",
            {root / "reader_static.dat"});
    });

    validate_band_reference_binding(
        contract, root.string(), root.string(), "band_kpath_info");
    assert_throws([&] {
        validate_band_reference_binding(
            contract, root.string(), root.string(), "wrong_kpath_info");
    });
}

void test_reader_path_conventions()
{
    const std::filesystem::path root =
        "test_qsgw_reader_paths.tmp";
    const auto band = resolve_band_reference_paths(
        root.string(), "band_kpath_info", 2);
    const std::filesystem::path absolute_root =
        std::filesystem::absolute(root).lexically_normal();
    assert(band.kpoints == absolute_root / "band_kpath_info");
    assert(band.eigenvalues.front() ==
           absolute_root / "band_KS_eigenvalue_k_00001.txt");
    assert(band.wavefunctions.back() ==
           absolute_root / "band_KS_eigenvector_k_00002.txt");

    const auto aims = resolve_same_grid_velocity_paths(
        QsgwProducer::FhiAims, root.string(), 2);
    assert(aims.size() == 2);
    assert(aims.front() ==
           absolute_root / "mommat_ks_kpt_000001.dat");
    assert(aims.back() ==
           absolute_root / "mommat_ks_kpt_000002.dat");
}

void test_contract_hashes_are_checked()
{
    const std::filesystem::path root =
        "test_qsgw_contract_hashes.tmp";
    std::filesystem::remove_all(root);
    std::filesystem::create_directories(root);
    const QsgwInputContract contract = parse(supported_contract(true, true));
    for (const std::string& role : {
             "mf0_eigenvalues", "mf0_wavefunctions", "scf_kpoints",
             "vxc_scf_manifest", "reader_static", "velocity_mf0",
             "band_kpoints", "vxc_band_manifest",
             "band_mf0_eigenvalues", "band_mf0_wavefunctions"})
    {
        for (const auto& file : contract.files(role))
        {
            std::ofstream output(root / file.file, std::ios::binary);
            output << "abc";
        }
    }
    contract.validate_file_hashes(root.string());
    {
        std::ofstream output(root / "vxc_scf.manifest",
                             std::ios::binary | std::ios::trunc);
        output << "changed";
    }
    assert_throws([&] { contract.validate_file_hashes(root.string()); });
    std::filesystem::remove_all(root);
}

} // namespace

int main()
{
    test_supported_head_only_band_contract();
    test_disabled_head_and_band_contract();
    test_unsupported_modes_are_rejected();
    test_reader_bindings_match_exact_files();
    test_reader_path_conventions();
    test_contract_hashes_are_checked();
    std::cout << "test_qsgw_input_contract: all tests passed\n";
    return 0;
}
