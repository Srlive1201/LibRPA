#ifdef NDEBUG
#undef NDEBUG
#endif

#include "../qsgw/vxc_io.h"
#include "../qsgw/sha256.h"

#include <cassert>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <limits>
#include <cstdio>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using librpa_int::Vector3_Order;
using librpa_int::ComplexMatrix;
using librpa_int::MeanField;
using librpa_int::Matz;
using librpa_int::cplxdb;
using librpa_int::qsgw::VxcDatasetKind;
using librpa_int::qsgw::VxcBasis;
using librpa_int::qsgw::VxcManifest;
using librpa_int::qsgw::VxcGauge;
using librpa_int::qsgw::VxcUnits;
using librpa_int::qsgw::prepare_vxc_in_fixed_state_basis;
using librpa_int::qsgw::project_vxc_nao_to_fixed_basis;
using librpa_int::qsgw::read_abacus_vxc_ha;

namespace
{

void assert_close(const cplxdb actual, const cplxdb expected)
{
    assert(std::abs(actual - expected) < 1.0e-14);
}

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

void test_abacus_vxc_uses_fixed_rydberg_to_hartree_conversion()
{
    const std::string text =
        "# rows 2\n"
        "# columns 2\n"
        "Row 1\n"
        " (2.0,0.0) (1.0,2.0)\n"
        "Row 2\n"
        " (4.0,0.0)\n";

    for (const int nspins: {1, 2})
    {
        // nspins is intentionally not an input to the parser.  ABACUS Vxc is
        // in Ry for both scalar and spin-polarized calculations.
        (void)nspins;
        std::istringstream input(text);
        const auto matrix = read_abacus_vxc_ha(input, "synthetic-vxck");
        assert_close(matrix(0, 0), 1.0);
        assert_close(matrix(0, 1), cplxdb(0.5, 1.0));
        assert_close(matrix(1, 0), cplxdb(0.5, -1.0));
        assert_close(matrix(1, 1), 2.0);
    }
}

void test_abacus_vxc_accepts_legacy_dimension_first_triangular_layout()
{
    const std::string text =
        "3\n"
        " (2.0,0.0) (1.0,2.0) (3.0,-1.0)\n"
        " (4.0,0.0) (5.0,6.0)\n"
        " (8.0,0.0)\n";
    std::istringstream input(text);
    const auto matrix = read_abacus_vxc_ha(input, "legacy-vxcs");

    assert_close(matrix(0, 0), 1.0);
    assert_close(matrix(0, 1), cplxdb(0.5, 1.0));
    assert_close(matrix(1, 0), cplxdb(0.5, -1.0));
    assert_close(matrix(0, 2), cplxdb(1.5, -0.5));
    assert_close(matrix(2, 0), cplxdb(1.5, 0.5));
    assert_close(matrix(1, 1), 2.0);
    assert_close(matrix(1, 2), cplxdb(2.5, 3.0));
    assert_close(matrix(2, 1), cplxdb(2.5, -3.0));
    assert_close(matrix(2, 2), 4.0);

    std::istringstream incomplete(
        "2\n"
        " (2.0,0.0) (1.0,0.0)\n");
    assert_throws([&] {
        (void)read_abacus_vxc_ha(incomplete, "incomplete-legacy-vxcs");
    });
}

void test_abacus_nao_vxc_is_projected_to_the_fixed_state_basis()
{
    MeanField reference(1, 1, 2, 2, 1);
    const double inv_sqrt_two = 1.0 / std::sqrt(2.0);
    ComplexMatrix wfc(2, 2);
    wfc(0, 0) = inv_sqrt_two;
    wfc(0, 1) = cplxdb(0.0, inv_sqrt_two);
    wfc(1, 0) = cplxdb(0.0, inv_sqrt_two);
    wfc(1, 1) = inv_sqrt_two;
    reference.get_eigenvectors()[0][0][0] = wfc;

    Matz vxc_nao(2, 2);
    vxc_nao(0, 0) = 2.0;
    vxc_nao(0, 1) = cplxdb(1.0, 0.5);
    vxc_nao(1, 0) = cplxdb(1.0, -0.5);
    vxc_nao(1, 1) = 4.0;

    Matz coefficients(2, 2);
    for (int row = 0; row < 2; ++row)
    {
        for (int column = 0; column < 2; ++column)
        {
            coefficients(row, column) = wfc(row, column);
        }
    }
    const Matz expected =
        conj(coefficients) * vxc_nao * transpose(coefficients, false);
    const Matz actual =
        project_vxc_nao_to_fixed_basis(vxc_nao, reference, 0, 0);
    const Matz selected = prepare_vxc_in_fixed_state_basis(
        vxc_nao, VxcBasis::Nao, reference, 0, 0);

    for (int row = 0; row < 2; ++row)
    {
        for (int column = 0; column < 2; ++column)
        {
            assert_close(actual(row, column), expected(row, column));
            assert_close(selected(row, column), expected(row, column));
        }
    }
    assert(std::abs(actual(0, 1) - vxc_nao(0, 1)) > 1.0e-6);

    Matz vxc_state(2, 2);
    vxc_state(0, 0) = -1.0;
    vxc_state(0, 1) = cplxdb(0.2, 0.1);
    vxc_state(1, 0) = cplxdb(0.2, -0.1);
    vxc_state(1, 1) = 0.5;
    const Matz state_selected = prepare_vxc_in_fixed_state_basis(
        vxc_state, VxcBasis::State, reference, 0, 0);
    for (int row = 0; row < 2; ++row)
    {
        for (int column = 0; column < 2; ++column)
        {
            assert_close(state_selected(row, column), vxc_state(row, column));
        }
    }
}

void test_vxc_matrix_validation_rejects_invalid_numeric_data_and_layout()
{
    MeanField reference(1, 1, 2, 2, 1);
    ComplexMatrix identity(2, 2);
    identity(0, 0) = 1.0;
    identity(1, 1) = 1.0;
    reference.get_eigenvectors()[0][0][0] = identity;

    Matz wrong_shape(1, 1);
    wrong_shape(0, 0) = 1.0;
    assert_throws([&] {
        (void)prepare_vxc_in_fixed_state_basis(
            wrong_shape, VxcBasis::Nao, reference, 0, 0);
    });
    assert_throws([&] {
        (void)prepare_vxc_in_fixed_state_basis(
            wrong_shape, VxcBasis::State, reference, 0, 0);
    });

    Matz nonhermitian(2, 2);
    nonhermitian(0, 0) = 1.0;
    nonhermitian(0, 1) = cplxdb(0.0, 0.5);
    nonhermitian(1, 0) = cplxdb(0.0, 0.5);
    nonhermitian(1, 1) = 2.0;
    assert_throws([&] {
        (void)prepare_vxc_in_fixed_state_basis(
            nonhermitian, VxcBasis::State, reference, 0, 0);
    });

    Matz nonfinite(2, 2);
    nonfinite(0, 0) = std::numeric_limits<double>::quiet_NaN();
    nonfinite(1, 1) = 2.0;
    assert_throws([&] {
        (void)prepare_vxc_in_fixed_state_basis(
            nonfinite, VxcBasis::Nao, reference, 0, 0);
    });

    const std::string imaginary_diagonal =
        "# rows 2\n"
        "# columns 2\n"
        "Row 1\n"
        " (2.0,0.1) (1.0,0.0)\n"
        "Row 2\n"
        " (4.0,0.0)\n";
    std::istringstream triangular_input(imaginary_diagonal);
    assert_throws([&] {
        (void)read_abacus_vxc_ha(triangular_input, "imaginary-diagonal");
    });

    const std::string nonhermitian_dense =
        "# rows 2\n"
        "# columns 2\n"
        " (2.0,0.0) (1.0,1.0)\n"
        " (1.0,1.0) (4.0,0.0)\n";
    std::istringstream dense_input(nonhermitian_dense);
    assert_throws([&] {
        (void)read_abacus_vxc_ha(dense_input, "nonhermitian-dense");
    });

    const std::string nonfinite_text =
        "# rows 1\n"
        "# columns 1\n"
        "Row 1\n"
        " (nan,0.0)\n";
    std::istringstream nonfinite_input(nonfinite_text);
    assert_throws([&] {
        (void)read_abacus_vxc_ha(nonfinite_input, "nonfinite-vxck");
    });
}

std::string valid_grid_manifest()
{
    return
        "# librpa-qsgw-vxc-manifest-v2\n"
        "kind\tscf\n"
        "producer\tabacus\n"
        "units\tRy\n"
        "basis\tnao\n"
        "gauge\tao_bloch\n"
        "spin\tk_index\tkx\tky\tkz\trows\tcolumns\tsha256\tfile\n"
        "1\t1\t0.0\t0.0\t0.0\t2\t2\tba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad\tvxck1s1_nao.txt\n"
        "1\t2\t0.5\t0.5\t0.0\t2\t2\tba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad\tvxck2s1_nao.txt\n";
}

std::string valid_abacus_state_grid_manifest()
{
    return
        "# librpa-qsgw-vxc-manifest-v2\n"
        "kind\tscf\n"
        "producer\tabacus\n"
        "units\tRy\n"
        "basis\tstate\n"
        "gauge\tmf0_state\n"
        "spin\tk_index\tkx\tky\tkz\trows\tcolumns\tsha256\tfile\n"
        "1\t1\t0.0\t0.0\t0.0\t2\t2\tba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad\tvxck1s1_nao.txt\n";
}

std::string valid_band_manifest()
{
    return
        "# librpa-qsgw-vxc-manifest-v2\n"
        "kind\tband\n"
        "producer\tabacus\n"
        "units\tRy\n"
        "basis\tnao\n"
        "gauge\tao_bloch\n"
        "spin\tk_index\tkx\tky\tkz\trows\tcolumns\tsha256\tfile\n"
        "1\t1\t0.0\t0.0\t0.0\t2\t2\tba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad\tvxck1s1_nao.txt\n"
        "1\t2\t0.25\t0.25\t0.25\t2\t2\tba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad\tvxck2s1_nao.txt\n";
}

void test_manifest_binds_each_matrix_to_spin_index_and_k_coordinate()
{
    std::istringstream input(valid_grid_manifest());
    const VxcManifest manifest = VxcManifest::parse(input, "grid-manifest");
    const std::vector<Vector3_Order<double>> expected_kpoints{
        {0.0, 0.0, 0.0},
        {0.5, 0.5, 0.0},
    };
    manifest.validate(
        VxcDatasetKind::ScfGrid, 1, expected_kpoints, 2, 2, 1.0e-10);
    assert(manifest.producer() == "abacus");
    assert(manifest.units() == VxcUnits::Rydberg);
    assert(manifest.basis() == VxcBasis::Nao);
    assert(manifest.gauge() == VxcGauge::AoBloch);
    assert(manifest.at(0, 0).file == "vxck1s1_nao.txt");
    assert(manifest.at(0, 1).file == "vxck2s1_nao.txt");
}

void test_abacus_out_mat_xc_manifest_is_state_basis_without_projection()
{
    std::istringstream input(valid_abacus_state_grid_manifest());
    const VxcManifest manifest =
        VxcManifest::parse(input, "abacus-out-mat-xc-manifest");
    manifest.validate(
        VxcDatasetKind::ScfGrid, 1, {{0.0, 0.0, 0.0}}, 2, 2, 1.0e-10);
    assert(manifest.producer() == "abacus");
    assert(manifest.units() == VxcUnits::Rydberg);
    assert(manifest.basis() == VxcBasis::State);
    assert(manifest.gauge() == VxcGauge::Mf0State);

    MeanField reference(1, 1, 2, 2, 1);
    const double inv_sqrt_two = 1.0 / std::sqrt(2.0);
    ComplexMatrix wfc(2, 2);
    wfc(0, 0) = inv_sqrt_two;
    wfc(0, 1) = cplxdb(0.0, inv_sqrt_two);
    wfc(1, 0) = cplxdb(0.0, inv_sqrt_two);
    wfc(1, 1) = inv_sqrt_two;
    reference.get_eigenvectors()[0][0][0] = wfc;

    Matz vxc_state(2, 2);
    vxc_state(0, 0) = -1.0;
    vxc_state(0, 1) = cplxdb(0.2, 0.1);
    vxc_state(1, 0) = cplxdb(0.2, -0.1);
    vxc_state(1, 1) = 0.5;
    const Matz selected = prepare_vxc_in_fixed_state_basis(
        vxc_state, manifest.basis(), reference, 0, 0);
    for (int row = 0; row < 2; ++row)
    {
        for (int column = 0; column < 2; ++column)
        {
            assert_close(selected(row, column), vxc_state(row, column));
        }
    }
}

void test_fhi_aims_manifest_marks_hartree_state_basis_without_projection()
{
    const std::string text =
        "# librpa-qsgw-vxc-manifest-v2\n"
        "kind\tscf\n"
        "producer\tfhi-aims\n"
        "units\tHa\n"
        "basis\tstate\n"
        "gauge\tmf0_state\n"
        "spin\tk_index\tkx\tky\tkz\trows\tcolumns\tsha256\tfile\n"
        "1\t1\t0.0\t0.0\t0.0\t2\t2\tba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad\txc_matr_spin_1_kpt_000001.csc\n";
    std::istringstream input(text);
    const VxcManifest manifest = VxcManifest::parse(input, "aims-manifest");
    manifest.validate(
        VxcDatasetKind::ScfGrid, 1, {{0.0, 0.0, 0.0}}, 2, 2, 1.0e-10);
    assert(manifest.producer() == "fhi-aims");
    assert(manifest.units() == VxcUnits::Hartree);
    assert(manifest.basis() == VxcBasis::State);
    assert(manifest.gauge() == VxcGauge::Mf0State);
}

void test_manifest_rejects_nscf_files_presented_as_scf_grid()
{
    std::string text = valid_grid_manifest();
    const auto position = text.find("0.5\t0.5\t0.0");
    text.replace(position, std::string("0.5\t0.5\t0.0").size(),
                 "0.05\t0.0\t0.0");
    std::istringstream input(text);
    const VxcManifest manifest = VxcManifest::parse(input, "wrong-grid-manifest");
    const std::vector<Vector3_Order<double>> expected_kpoints{
        {0.0, 0.0, 0.0},
        {0.5, 0.5, 0.0},
    };
    assert_throws([&] {
        manifest.validate(
            VxcDatasetKind::ScfGrid, 1, expected_kpoints, 2, 2, 1.0e-10);
    });
}

void test_manifest_rejects_duplicate_spin_k_entries()
{
    std::string text = valid_grid_manifest();
    text += "1\t2\t0.5\t0.5\t0.0\t2\t2\tba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad\tduplicate.txt\n";
    std::istringstream input(text);
    assert_throws([&] { (void)VxcManifest::parse(input, "duplicate-manifest"); });
}

void test_band_manifest_is_distinct_from_scf_grid()
{
    std::istringstream input(valid_band_manifest());
    const VxcManifest manifest = VxcManifest::parse(input, "band-manifest");
    const std::vector<Vector3_Order<double>> band_kpoints{
        {0.0, 0.0, 0.0},
        {0.25, 0.25, 0.25},
    };
    manifest.validate(
        VxcDatasetKind::BandPath, 1, band_kpoints, 2, 2, 1.0e-10);
    assert_throws([&] {
        manifest.validate(
            VxcDatasetKind::ScfGrid, 1, band_kpoints, 2, 2, 1.0e-10);
    });
}

void test_manifest_rejects_missing_spin_or_kpoint_entries()
{
    std::istringstream input(valid_grid_manifest());
    const VxcManifest manifest = VxcManifest::parse(input, "incomplete-spin-manifest");
    const std::vector<Vector3_Order<double>> expected_kpoints{
        {0.0, 0.0, 0.0},
        {0.5, 0.5, 0.0},
    };
    assert_throws([&] {
        manifest.validate(
            VxcDatasetKind::ScfGrid, 2, expected_kpoints, 2, 2, 1.0e-10);
    });

    std::string missing_kpoint = valid_grid_manifest();
    const auto last_entry = missing_kpoint.rfind("1\t2\t");
    missing_kpoint.erase(last_entry);
    std::istringstream missing_input(missing_kpoint);
    const VxcManifest missing =
        VxcManifest::parse(missing_input, "incomplete-kpoint-manifest");
    assert_throws([&] {
        missing.validate(
            VxcDatasetKind::ScfGrid, 1, expected_kpoints, 2, 2, 1.0e-10);
    });
}

void test_manifest_rejects_unsupported_units()
{
    std::string text = valid_grid_manifest();
    const auto units = text.find("units\tRy");
    text.replace(units, std::string("units\tRy").size(), "units\teV");
    std::istringstream input(text);
    assert_throws([&] {
        (void)VxcManifest::parse(input, "wrong-units-manifest");
    });
}

void test_manifest_rejects_incompatible_producer_units_or_basis()
{
    const auto assert_bad_manifest = [](std::string text,
                                        const std::string& needle,
                                        const std::string& replacement) {
        const auto position = text.find(needle);
        assert(position != std::string::npos);
        text.replace(position, needle.size(), replacement);
        std::istringstream input(text);
        assert_throws([&] {
            (void)VxcManifest::parse(input, "incompatible-vxc-manifest");
        });
    };

    assert_bad_manifest(valid_grid_manifest(), "units\tRy", "units\tHa");
    assert_bad_manifest(valid_grid_manifest(), "basis\tnao", "basis\tstate");
    assert_bad_manifest(valid_grid_manifest(), "gauge\tao_bloch", "gauge\tmf0_state");

    const std::string aims_manifest =
        "# librpa-qsgw-vxc-manifest-v2\n"
        "kind\tscf\n"
        "producer\tfhi-aims\n"
        "units\tHa\n"
        "basis\tstate\n"
        "gauge\tmf0_state\n"
        "spin\tk_index\tkx\tky\tkz\trows\tcolumns\tsha256\tfile\n"
        "1\t1\t0.0\t0.0\t0.0\t2\t2\tba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad\txc_matr_spin_1_kpt_000001.csc\n";
    assert_bad_manifest(aims_manifest, "units\tHa", "units\tRy");
    assert_bad_manifest(aims_manifest, "basis\tstate", "basis\tnao");
}

void test_manifest_rejects_shape_hash_gauge_and_changed_files()
{
    const auto assert_bad = [](std::string text,
                               const std::string& needle,
                               const std::string& replacement) {
        const auto position = text.find(needle);
        assert(position != std::string::npos);
        text.replace(position, needle.size(), replacement);
        std::istringstream input(text);
        assert_throws([&] { (void)VxcManifest::parse(input, "bad-v2"); });
    };
    assert_bad(valid_grid_manifest(), "gauge\tao_bloch", "gauge\tmf0_state");
    assert_bad(valid_grid_manifest(),
               "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",
               "not-a-sha256");

    std::istringstream shape_input(valid_grid_manifest());
    const VxcManifest shape_manifest =
        VxcManifest::parse(shape_input, "shape-v2");
    assert_throws([&] {
        shape_manifest.validate(
            VxcDatasetKind::ScfGrid, 1,
            {{0.0, 0.0, 0.0}, {0.5, 0.5, 0.0}}, 3, 3, 1.0e-10);
    });

    const std::string path = "test_qsgw_vxc_hash_input.tmp";
    {
        std::ofstream output(path, std::ios::binary);
        output << "abc";
    }
    std::string one_file = valid_grid_manifest();
    one_file.erase(one_file.rfind("1\t2\t"));
    const auto file_position = one_file.find("vxck1s1_nao.txt");
    one_file.replace(file_position, std::string("vxck1s1_nao.txt").size(), path);
    std::istringstream file_input(one_file);
    const VxcManifest file_manifest =
        VxcManifest::parse(file_input, "hash-v2");
    file_manifest.validate_file_hashes(".");
    {
        std::ofstream output(path, std::ios::binary | std::ios::app);
        output << "changed";
    }
    assert_throws([&] { file_manifest.validate_file_hashes("."); });
    std::remove(path.c_str());
}

} // namespace

int main()
{
    test_abacus_vxc_uses_fixed_rydberg_to_hartree_conversion();
    test_abacus_vxc_accepts_legacy_dimension_first_triangular_layout();
    test_abacus_nao_vxc_is_projected_to_the_fixed_state_basis();
    test_vxc_matrix_validation_rejects_invalid_numeric_data_and_layout();
    test_manifest_binds_each_matrix_to_spin_index_and_k_coordinate();
    test_abacus_out_mat_xc_manifest_is_state_basis_without_projection();
    test_fhi_aims_manifest_marks_hartree_state_basis_without_projection();
    test_manifest_rejects_nscf_files_presented_as_scf_grid();
    test_manifest_rejects_duplicate_spin_k_entries();
    test_band_manifest_is_distinct_from_scf_grid();
    test_manifest_rejects_missing_spin_or_kpoint_entries();
    test_manifest_rejects_unsupported_units();
    test_manifest_rejects_incompatible_producer_units_or_basis();
    test_manifest_rejects_shape_hash_gauge_and_changed_files();
    std::cout << "test_qsgw_vxc_io: all tests passed\n";
    return 0;
}
