#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstddef>
#include <stdexcept>

#include "../core/pbc.h"
#include "../api/compute_helper.h"
#include "../utils/constants.h"

using namespace librpa_int;

static void test_is_gamma_point()
{
    assert(is_gamma_point(Vector3_Order<double>{0.0, 0.0, 0.0}));
    assert(!is_gamma_point(Vector3_Order<double>{0.3, 0.2, 0.1}));
    assert(is_gamma_point(Vector3_Order<int>{0, 0, 0}));
}

static void test_get_R_index()
{
    Vector3_Order<int> period{2, 2, 2};
    std::vector<Vector3_Order<int>> sc222 = construct_R_grid(period);
    for (size_t i = 0; i != sc222.size(); i++ )
        printf("i=%zu: x %d, y %d, z %d\n", i, sc222[i].x, sc222[i].y, sc222[i].z);
    Vector3_Order<int> R {-1, -1, -1};
    printf(" R: x %d, y %d, z %d\n", R.x, R.y, R.z);
    auto mR = -R;
    printf("-R: x %d, y %d, z %d\n", mR.x, mR.y, mR.z);
    assert(get_R_index(sc222, R) == 0);
    assert(get_R_index(sc222, Vector3_Order<int>{-1, 0, -1}) == 2);
    assert(get_R_index(sc222, Vector3_Order<int>{3, 3, -1}) < 0);
    assert(get_R_index(sc222, Vector3_Order<int>{3, 3, -1} % period) == 0);
}

static void test_periodic_boundary_data()
{
    PeriodicBoundaryData pbc;
    const Vector3_Order<int> gamma_period{1, 1, 1};
    const std::array<int, 3> gamma_period_array{1, 1, 1};
    assert(pbc.period == gamma_period);
    assert(pbc.period_array == gamma_period_array);
    assert(pbc.Rlist == construct_R_grid(pbc.period));

    pbc.set_period(2, 3, 1);
    const Vector3_Order<int> period{2, 3, 1};
    const std::array<int, 3> period_array{2, 3, 1};
    assert(pbc.period == period);
    assert(pbc.period_array == period_array);
    assert(pbc.Rlist == construct_R_grid(pbc.period));

    pbc.set_latvec_and_G({1, 2, 3, 4, 5, 6, 7, 8, 9},
                         {1, 2, 3, 4, 5, 6, 7, 8, 9});
}

static void test_kgrids_with_weighted_coulomb_mapping()
{
    PeriodicBoundaryData pbc;
    pbc.set_latvec({1, 0, 0, 0, 1, 0, 0, 0, 1});

    const std::vector<double> kvecs{
        0.0, 0.0, 0.0,
        librpa_int::TWO_PI / 3.0, 0.0, 0.0,
        librpa_int::TWO_PI * 2.0 / 3.0, 0.0, 0.0,
    };
    pbc.set_kgrids_kvec(3, 1, 1, kvecs);
    pbc.set_kq_mapping({0, 0, 2});

    assert(pbc.klist.size() == 3);
    assert(pbc.klist_full.size() == 3);
    assert(std::abs(pbc.weight_k[0] - 1.0 / 3.0) < 1e-12);
    assert(std::abs(pbc.weight_k[1] - 1.0 / 3.0) < 1e-12);
    assert(std::abs(pbc.weight_k[2] - 1.0 / 3.0) < 1e-12);
    assert(pbc.klist_coul.size() == 2);
    assert(std::abs(pbc.weight_q[0] - 2.0 / 3.0) < 1e-12);
    assert(std::abs(pbc.weight_q[1] - 1.0 / 3.0) < 1e-12);
    assert(pbc.map_irk_ks.at(pbc.klist_coul[0]).size() == 2);
    assert(pbc.map_irk_ks.at(pbc.klist_coul[1]).size() == 1);

    pbc.set_period(3, 1, 1);
    assert(pbc.klist.size() == 3);
    pbc.set_period(1, 1, 1);
    assert(pbc.klist.empty());
    assert(pbc.klist_full.empty());
    assert(pbc.klist_coul.empty());
}

static void test_full_scf_kgrids_keep_loaded_order()
{
    PeriodicBoundaryData pbc;
    pbc.set_latvec({1, 0, 0, 0, 1, 0, 0, 0, 1});

    const std::vector<double> kvecs{
        0.0, 0.0, 0.0,
        librpa_int::TWO_PI * 2.0 / 3.0, 0.0, 0.0,
        librpa_int::TWO_PI / 3.0, 0.0, 0.0,
    };
    pbc.set_kgrids_kvec(3, 1, 1, kvecs);

    assert(pbc.klist_full == pbc.klist);
    assert(pbc.kfrac_list_full == pbc.kfrac_list);
    assert(pbc.k_to_kfull == std::vector<int>({0, 1, 2}));
    assert(pbc.kfull_to_k == std::vector<int>({0, 1, 2}));
    assert(pbc.kfull_to_k_relation == std::vector<KFullToKRelation>({
        KFullToKRelation::DIRECT,
        KFullToKRelation::DIRECT,
        KFullToKRelation::DIRECT,
    }));
    assert(!pbc.kgrid_uses_time_reversal);
}

static void test_full_scf_kgrids_canonicalize_rounded_cartesian_coordinates()
{
    PeriodicBoundaryData pbc;
    pbc.set_latvec({1, 0, 0, 0, 1, 0, 0, 0, 1});

    std::vector<double> kvecs;
    for (int ik = 0; ik != 8; ++ik)
    {
        const double exact_kfrac = static_cast<double>(ik) / 8.0;
        const double rounded_kfrac = exact_kfrac - static_cast<double>(ik) * 7.0e-7;
        kvecs.insert(kvecs.end(), {TWO_PI * rounded_kfrac, 0.0, 0.0});
    }
    pbc.set_kgrids_kvec(8, 1, 1, kvecs);

    for (int ik = 0; ik != 8; ++ik)
    {
        const double exact_kfrac = static_cast<double>(ik) / 8.0;
        if (std::abs(pbc.kfrac_list[ik].x - exact_kfrac) >= 1e-15
            || std::abs(pbc.klist[ik].x - exact_kfrac) >= 1e-15)
        {
            throw std::runtime_error("BvK k-point was not canonicalized");
        }
    }
}

static void test_reduced_scf_kgrids()
{
    PeriodicBoundaryData pbc;
    pbc.set_latvec({1, 0, 0, 0, 1, 0, 0, 0, 1});

    const std::vector<double> kvecs{
        0.0, 0.0, 0.0,
        librpa_int::TWO_PI * (1.0 / 3.0 - 2.0e-6), 0.0, 0.0,
    };
    pbc.set_kgrids_kvec(3, 1, 1, kvecs, {1.0 / 3.0, 2.0 / 3.0});
    pbc.set_kq_mapping({0, 1});

    assert(pbc.klist.size() == 2);
    assert(pbc.klist_full.size() == 3);
    assert(pbc.k_to_kfull == std::vector<int>({0, 1}));
    assert(pbc.kfull_to_k == std::vector<int>({0, 1, 1}));
    assert(pbc.kfull_to_k_relation == std::vector<KFullToKRelation>({
        KFullToKRelation::DIRECT,
        KFullToKRelation::DIRECT,
        KFullToKRelation::TIME_REVERSAL,
    }));
    assert(pbc.kgrid_uses_time_reversal);
    assert(pbc.klist_coul.size() == 2);
    if (std::abs(pbc.kfrac_list[1].x - 1.0 / 3.0) >= 1e-12)
    {
        throw std::runtime_error("Reduced BvK k-point was not canonicalized");
    }
    assert(std::abs(pbc.weight_q[0] - 1.0 / 3.0) < 1e-12);
    assert(std::abs(pbc.weight_q[1] - 2.0 / 3.0) < 1e-12);
}

static void test_incomplete_time_reversal_reduced_scf_kgrids()
{
    PeriodicBoundaryData pbc;
    pbc.set_latvec({1, 0, 0, 0, 1, 0, 0, 0, 1});

    const std::vector<double> kvecs{
        0.0, 0.0, 0.0,
        librpa_int::TWO_PI * 0.25, 0.0, 0.0,
    };
    pbc.set_kgrids_kvec(4, 1, 1, kvecs, {0.25, 0.75});

    assert(pbc.kfull_to_k == std::vector<int>({0, 1, -1, 1}));
    assert(pbc.weight_k == std::vector<double>({0.25, 0.75}));
    assert(pbc.kfull_to_k_relation == std::vector<KFullToKRelation>({
        KFullToKRelation::DIRECT,
        KFullToKRelation::DIRECT,
        KFullToKRelation::NONE,
        KFullToKRelation::TIME_REVERSAL,
    }));
    assert(!pbc.kgrid_uses_time_reversal);

    pbc.set_kgrids_kvec(4, 1, 1, kvecs, {1.0, 3.0});
    assert(std::abs(pbc.weight_k[0] - 0.25) < 1e-12);
    assert(std::abs(pbc.weight_k[1] - 0.75) < 1e-12);

    PeriodicBoundaryData pbc_without_weights;
    pbc_without_weights.set_latvec({1, 0, 0, 0, 1, 0, 0, 0, 1});
    pbc_without_weights.set_kgrids_kvec(4, 1, 1, kvecs);
    assert(pbc_without_weights.klist_coul.size() == 2);
    assert(pbc_without_weights.weight_k.empty());
    assert(pbc_without_weights.weight_q.empty());
    assert(pbc_without_weights.map_q_weight.empty());
}

static void test_irreducible_kgrids_from_symmetry_stars()
{
    PeriodicBoundaryData pbc;
    pbc.set_latvec({2, 0, 0, 0, 3, 0, 0, 0, 4});

    const std::vector<double> kvecs_ibz{
        0.0, 0.0, 0.0,
        librpa_int::TWO_PI / 6.0, 0.0, 0.0,
    };
    const std::vector<std::vector<Vector3_Order<double>>> full_kstars{
        {{0.0, 0.0, 0.0}},
        {{1.0 / 6.0, 0.0, 0.0}, {-1.0 / 6.0, 0.0, 0.0}},
    };
    pbc.set_irreducible_kgrids_kvec(3, 1, 1, kvecs_ibz, full_kstars);

    assert(pbc.klist.size() == 2);
    assert(pbc.klist_full.size() == 3);
    assert(std::abs(pbc.kfrac_list[1].x - 1.0 / 3.0) < 1e-12);
    assert(std::abs(pbc.weight_q[0] - 1.0 / 3.0) < 1e-12);
    assert(std::abs(pbc.weight_q[1] - 2.0 / 3.0) < 1e-12);
    assert(pbc.map_irk_ks.at(pbc.klist_coul[1]).size() == 2);
    assert(std::abs(pbc.klist_full[1].x - 1.0 / 6.0) < 1e-12);
    assert(std::abs(pbc.klist_full[2].x + 1.0 / 6.0) < 1e-12);
}

static void test_atom_pair_bvk_remap()
{
    typedef std::size_t atom_t;
    typedef std::pair<atom_t, atom_t> atpair_t;

    const std::map<atom_t, Vector3<double>> coord_fracs{
        {0, {0.1, 0.0, 0.0}},
        {1, {0.9, 0.0, 0.0}},
    };
    const std::vector<Vector3_Order<int>> Rs{
        {0, 0, 0},
        {1, 0, 0},
    };
    const Vector3_Order<int> period{3, 3, 3};
    const Matrix3 latvec;

    const AtomPairBvKRemap<atom_t, atpair_t> remap(coord_fracs, Rs, period, latvec);
    const atpair_t pair00{0, 0};
    const atpair_t pair01{0, 1};
    const atpair_t pair10{1, 0};
    const Vector3_Order<int> R0{0, 0, 0};
    const Vector3_Order<int> R1{1, 0, 0};
    const Vector3_Order<int> Rm2{-2, 0, 0};

    const auto *R_bvk = remap.find_R_bvk(pair01, R1);
    assert(R_bvk != nullptr);
    assert(R_bvk->size() == 1);
    assert(R_bvk->front() == Rm2);
    assert(remap.find_R_bvk(pair00, R1) == nullptr);
    assert(remap.find_R_bvk(pair01, R0) == nullptr);
    assert(remap.find_R_bvk(pair10, R1) == nullptr);
    assert(remap.find_R_bvk({2, 0}, R1) == nullptr);

    const Vector3_Order<int> period_2{2, 2, 2};
    const Vector3_Order<int> Rm1{-1, 0, 0};
    const AtomPairBvKRemap<atom_t, atpair_t> remap_ws(coord_fracs, Rs, period_2, latvec, 1);
    const auto *R_bvks_ws = remap_ws.find_R_bvk(pair00, R1);
    assert(R_bvks_ws != nullptr);
    assert(R_bvks_ws->size() == 2);
    assert(std::find(R_bvks_ws->cbegin(), R_bvks_ws->cend(), R1) != R_bvks_ws->cend());
    assert(std::find(R_bvks_ws->cbegin(), R_bvks_ws->cend(), Rm1) != R_bvks_ws->cend());
    assert(remap_ws.find_R_bvk(pair00, R0) == nullptr);
}

using BandMatrix = std::array<std::complex<double>, 16>;

static BandMatrix hermitian_mesh_matrix(std::size_t ik, double shift)
{
    BandMatrix matrix{};
    for (std::size_t i = 0; i != 4; ++i)
    {
        for (std::size_t j = i; j != 4; ++j)
        {
            const auto real = std::sin((ik + 1.0 + shift) * (i + j + 1.0) * 0.37);
            const auto imag = i == j ? 0.0 : std::cos((ik + 0.5 + shift) * (j - i) * 0.61);
            matrix[4 * i + j] = {real, imag};
            matrix[4 * j + i] = std::conj(matrix[4 * i + j]);
        }
    }
    return matrix;
}

static BandMatrix interpolate_band_matrix(const std::vector<BandMatrix> &cells,
                                          const std::vector<Vector3_Order<int>> &Rs,
                                          const AtomPairBvKRemap<std::size_t> &remap,
                                          const Vector3_Order<double> &k)
{
    BandMatrix result{};
    for (std::size_t i = 0; i != 4; ++i)
    {
        for (std::size_t j = 0; j != 4; ++j)
        {
            for (std::size_t ir = 0; ir != Rs.size(); ++ir)
            {
                const auto *images = remap.find_R_bvk({i / 2, j / 2}, Rs[ir]);
                std::complex<double> phase{};
                if (images == nullptr)
                    phase = std::polar(1.0, TWO_PI * (k * Rs[ir]));
                else
                {
                    for (const auto &image : *images)
                        phase += std::polar(1.0, TWO_PI * (k * image));
                    phase /= static_cast<double>(images->size());
                }
                result[4 * i + j] += phase * cells[ir][4 * i + j];
            }
        }
    }
    return result;
}

static void test_band_remapping_preserves_adjoint_and_mesh()
{
    // The even mesh includes same-atom cells equal to their opposite modulo
    // the period. The skew lattice additionally exercises corner image ties.
    const Vector3_Order<int> period{4, 4, 4};
    const auto Rs = construct_R_grid(period);
    const auto ks = build_uniform_kmesh_frac(period);
    const std::map<std::size_t, Vector3<double>> coords{{0, {0.0, 0.0, 0.0}},
                                                        {1, {0.25, 0.25, 0.25}}};
    const std::array<Matrix3, 2> lattices{Matrix3{},
                                          Matrix3{0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0}};
    const std::array<Vector3_Order<double>, 3> off_mesh{Vector3_Order<double>{0.173, 0.271, -0.119},
                                                        Vector3_Order<double>{-0.31, 0.07, 0.43},
                                                        Vector3_Order<double>{0.125, 0.2, 0.375}};

    // Test the Hermitian real and imaginary operators separately: a general
    // frequency-domain self-energy H + iK is not itself Hermitian.
    for (const double shift : {0.0, 0.43})
    {
        std::vector<BandMatrix> mesh;
        for (std::size_t ik = 0; ik != ks.size(); ++ik)
            mesh.push_back(hermitian_mesh_matrix(ik, shift));
        std::vector<BandMatrix> cells(Rs.size());
        for (std::size_t ir = 0; ir != Rs.size(); ++ir)
        {
            for (std::size_t ik = 0; ik != ks.size(); ++ik)
            {
                const auto phase =
                    std::polar(1.0, -TWO_PI * (ks[ik] * Rs[ir])) / static_cast<double>(ks.size());
                for (std::size_t e = 0; e != 16; ++e) cells[ir][e] += phase * mesh[ik][e];
            }
        }
        for (const auto &lattice : lattices)
        {
            for (const int mode : {0, 1})
            {
                const AtomPairBvKRemap<std::size_t> remap(coords, Rs, period, lattice, mode);
                for (std::size_t ik = 0; ik != ks.size(); ++ik)
                {
                    const auto actual = interpolate_band_matrix(cells, Rs, remap, ks[ik]);
                    for (std::size_t e = 0; e != 16; ++e)
                        assert(std::abs(actual[e] - mesh[ik][e]) < 5e-12);
                }
                for (const auto &k : off_mesh)
                {
                    const auto actual = interpolate_band_matrix(cells, Rs, remap, k);
                    for (std::size_t i = 0; i != 4; ++i)
                        for (std::size_t j = 0; j != 4; ++j)
                            assert(std::abs(actual[4 * i + j] - std::conj(actual[4 * j + i])) <
                                   5e-12);
                }
            }
        }
    }
}

using AtomMatrix = std::array<std::complex<double>, 4>;

static AtomMatrix interpolate_displaced_atom(const std::vector<Vector3_Order<int>> &Rs,
                                             const std::vector<AtomMatrix> &cells,
                                             const AtomPairBvKRemap<atom_t> &remap,
                                             const Vector3_Order<double> &k)
{
    AtomMatrix result{};
    for (std::size_t ir = 0; ir != Rs.size(); ++ir)
    {
        const auto *images = remap.find_R_bvk({atom_t{1}, atom_t{1}}, Rs[ir]);
        std::complex<double> phase{};
        if (images == nullptr)
            phase = std::polar(1.0, TWO_PI * (k * Rs[ir]));
        else
        {
            for (const auto &image : *images) phase += std::polar(1.0, TWO_PI * (k * image));
            phase /= static_cast<double>(images->size());
        }
        for (std::size_t e = 0; e != result.size(); ++e) result[e] += phase * cells[ir][e];
    }
    return result;
}

static void test_band_remap_api_preserves_displaced_atom_adjoint_and_mesh()
{
    // Exercise the helper shared by EXX, G0W0 and the active QSGW band path,
    // including Cartesian atom input converted to fractional coordinates.
    constexpr double a = 5.1315533365962613;
    const Matrix3 lattice(0.0, a, a, a, 0.0, a, a, a, 0.0);
    PeriodicBoundaryData pbc;
    pbc.set_latvec({0.0, a, a, a, 0.0, a, a, a, 0.0});
    pbc.set_period(4, 4, 4);
    Atoms atoms;
    atoms.set({14, 14},
              {Vector3<double>{0.0, 0.0, 0.0}, Vector3<double>{0.5 * a, 0.5 * a, 0.5 * a}},
              lattice);

    std::vector<AtomMatrix> cells(pbc.Rlist.size());
    const Vector3_Order<int> corner{1, 1, 1};
    const auto positive = static_cast<std::size_t>(pbc.get_R_index(corner));
    const auto negative = static_cast<std::size_t>(pbc.get_R_index(-corner));
    cells[positive] = AtomMatrix{{{0.3, 0.1}, {0.7, 0.2}, {-0.4, 0.6}, {0.1, -0.2}}};
    for (std::size_t i = 0; i != 2; ++i)
        for (std::size_t j = 0; j != 2; ++j)
            cells[negative][2 * i + j] = std::conj(cells[positive][2 * j + i]);
    const auto boundary = static_cast<std::size_t>(pbc.get_R_index(Vector3_Order<int>{-2, 0, 0}));
    cells[boundary] = AtomMatrix{{{0.8, 0.0}, {0.25, 0.17}, {0.25, -0.17}, {-0.3, 0.0}}};

    const std::array<Vector3_Order<double>, 2> off_mesh{Vector3_Order<double>{0.173, 0.271, -0.119},
                                                        Vector3_Order<double>{-0.31, 0.07, 0.43}};
    const auto mesh = build_uniform_kmesh_frac(pbc.period);
    const AtomPairBvKRemap<atom_t> identity;
    for (const int mode : {0, 1})
    {
        const auto remap = api::build_band_bvk_remap(atoms, pbc, mode);
        for (const auto &k : mesh)
        {
            const auto expected = interpolate_displaced_atom(pbc.Rlist, cells, identity, k);
            const auto actual = interpolate_displaced_atom(pbc.Rlist, cells, remap, k);
            for (std::size_t e = 0; e != actual.size(); ++e)
                assert(std::abs(actual[e] - expected[e]) < 5e-12);
        }
        for (const auto &k : off_mesh)
        {
            const auto actual = interpolate_displaced_atom(pbc.Rlist, cells, remap, k);
            for (std::size_t i = 0; i != 2; ++i)
                for (std::size_t j = 0; j != 2; ++j)
                    assert(std::abs(actual[2 * i + j] - std::conj(actual[2 * j + i])) < 5e-12);
        }
    }
}

int main(int argc, char *argv[])
{
    test_is_gamma_point();
    test_get_R_index();
    test_periodic_boundary_data();
    test_kgrids_with_weighted_coulomb_mapping();
    test_full_scf_kgrids_keep_loaded_order();
    test_full_scf_kgrids_canonicalize_rounded_cartesian_coordinates();
    test_reduced_scf_kgrids();
    test_incomplete_time_reversal_reduced_scf_kgrids();
    test_irreducible_kgrids_from_symmetry_stars();
    test_atom_pair_bvk_remap();
    test_band_remapping_preserves_adjoint_and_mesh();
    test_band_remap_api_preserves_displaced_atom_adjoint_and_mesh();
    return 0;
}
