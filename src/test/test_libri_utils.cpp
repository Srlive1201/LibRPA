#include <array>
#include <cassert>
#include <complex>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include "../math/vector3_order.h"
#include "../utils/libri_utils.h"

#ifdef LIBRPA_USE_LIBRI
#include "RI/physics/GW.h"
#endif

void test_comm_map2_targets_from_vector()
{
    using namespace librpa_int;

    const std::vector<std::pair<int, int>> atpairs{{0, 0}, {0, 1}, {2, 3}};
    const std::vector<Vector3_Order<int>> cells{{0, 0, 0}, {1, -1, 2}};
    const auto targets = get_s0_s1_for_comm_map2<int, int>(atpairs, cells);

    assert(targets.first.size() == 2);
    assert(targets.first.count(0) == 1);
    assert(targets.first.count(2) == 1);

    assert(targets.second.size() == 6);
    assert(targets.second.count({0, {0, 0, 0}}) == 1);
    assert(targets.second.count({0, {1, -1, 2}}) == 1);
    assert(targets.second.count({1, {0, 0, 0}}) == 1);
    assert(targets.second.count({1, {1, -1, 2}}) == 1);
    assert(targets.second.count({3, {0, 0, 0}}) == 1);
    assert(targets.second.count({3, {1, -1, 2}}) == 1);
}

void test_comm_map2_targets_from_set()
{
    using namespace librpa_int;

    const std::set<std::pair<int, int>> atpairs{{0, 0}, {0, 1}, {2, 3}};
    const std::vector<Vector3_Order<int>> cells{{0, 0, 0}, {1, -1, 2}};
    const auto targets = get_s0_s1_for_comm_map2<int, int>(atpairs, cells);

    const std::set<int> expected_s0{0, 2};
    const std::set<libri_types<int, int>::TAC> expected_s1{
        {0, {0, 0, 0}}, {0, {1, -1, 2}},
        {1, {0, 0, 0}}, {1, {1, -1, 2}},
        {3, {0, 0, 0}}, {3, {1, -1, 2}},
    };

    assert(targets.first == expected_s0);
    assert(targets.second == expected_s1);
}

void test_comm_map2_targets_output_types()
{
    using namespace librpa_int;

    const std::vector<std::pair<int, int>> atpairs{{1, 2}};
    const std::vector<Vector3_Order<int>> cells{{3, 4, 5}};
    const auto targets = get_s0_s1_for_comm_map2<int, int, long, long>(atpairs, cells);

    const std::set<long> expected_s0{1L};
    const std::set<libri_types<long, long>::TAC> expected_s1{{2L, {3L, 4L, 5L}}};
    assert(targets.first == expected_s0);
    assert(targets.second == expected_s1);
}

void test_comm_map2_targets_empty_cells()
{
    using namespace librpa_int;

    const std::vector<std::pair<int, int>> atpairs{{0, 1}, {2, 3}};
    const std::vector<Vector3_Order<int>> cells;
    const auto targets = get_s0_s1_for_comm_map2<int, int>(atpairs, cells);

    const std::set<int> expected_s0{0, 2};
    assert(targets.first == expected_s0);
    assert(targets.second.empty());
}

#ifdef LIBRPA_USE_LIBRI
void test_libri_gw_scalar_complex_contraction()
{
    using complex_t = std::complex<double>;
    using tensor_t = RI::Tensor<complex_t>;
    using tensor_map_t =
        std::map<int, std::map<std::pair<int, std::array<int, 3>>, tensor_t>>;

    const std::array<int, 3> origin{0, 0, 0};
    const complex_t c{1.25, 0.0};
    const complex_t w{0.75, -0.2};
    const complex_t g{-0.4, 0.3};

    tensor_t c_tensor({1, 1, 1});
    tensor_t w_tensor({1, 1});
    tensor_t g_tensor({1, 1});
    c_tensor(0, 0, 0) = c;
    w_tensor(0, 0) = w;
    g_tensor(0, 0) = g;

    tensor_map_t cs{{0, {{{0, origin}, c_tensor}}}};
    tensor_map_t ws{{0, {{{0, origin}, w_tensor}}}};
    tensor_map_t gs{{0, {{{0, origin}, g_tensor}}}};

    RI::GW<int, int, 3, complex_t> gw;
    const std::array<std::array<double, 3>, 3> lattice{{
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
    }};
    gw.set_parallel(MPI_COMM_WORLD, {{0, {0.0, 0.0, 0.0}}},
                    lattice, {1, 1, 1});
    gw.set_symmetry(false, {});
    gw.set_Cs(cs, 0.0);
    gw.set_Ws(ws, 0.0);
    gw.set_Gs(gs, 0.0);
    gw.cal_Sigmas();

    const auto &sigma = gw.Sigmas.at(0).at({0, origin});
    assert(sigma.shape.size() == 2);
    assert(sigma.shape[0] == 1);
    assert(sigma.shape[1] == 1);
    const complex_t expected = 4.0 * c * c * w * g;
    assert(std::abs(sigma(0, 0) - expected) < 1e-12);
}

void test_libri_gw_scalar_translation_contraction()
{
    using complex_t = std::complex<double>;
    using tensor_t = RI::Tensor<complex_t>;
    using tensor_map_t =
        std::map<int, std::map<std::pair<int, std::array<int, 3>>, tensor_t>>;

    const std::array<int, 3> origin{0, 0, 0};
    const std::array<int, 3> r1{1, 0, 0};
    const complex_t c{0.8, 0.0};
    const complex_t w{0.6, -0.15};
    const complex_t g{-0.35, 0.22};

    tensor_t c_tensor({1, 1, 1});
    tensor_t w_tensor({1, 1});
    tensor_t g_tensor({1, 1});
    c_tensor(0, 0, 0) = c;
    w_tensor(0, 0) = w;
    g_tensor(0, 0) = g;

    tensor_map_t cs{{0, {{{0, origin}, c_tensor}}}};
    tensor_map_t ws{{0, {{{0, r1}, w_tensor}}}};
    tensor_map_t gs{{0, {{{0, r1}, g_tensor}}}};

    RI::GW<int, int, 3, complex_t> gw;
    const std::array<std::array<double, 3>, 3> lattice{{
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
    }};
    gw.set_parallel(MPI_COMM_WORLD, {{0, {0.0, 0.0, 0.0}}}, lattice, {3, 1, 1});
    gw.set_symmetry(false, {});
    gw.set_Cs(cs, 0.0);
    gw.set_Ws(ws, 0.0);
    gw.set_Gs(gs, 0.0);
    gw.cal_Sigmas();

    const auto &sigma = gw.Sigmas.at(0).at({0, r1});
    const complex_t expected = 4.0 * c * c * w * g;
    assert(std::abs(sigma(0, 0) - expected) < 1e-12);
}

void test_libri_gw_scalar_mesh_matches_k_minus_q_convolution()
{
    using complex_t = std::complex<double>;
    using tensor_t = RI::Tensor<complex_t>;
    using tensor_map_t =
        std::map<int, std::map<std::pair<int, std::array<int, 3>>, tensor_t>>;

    constexpr int mesh = 3;
    const std::array<int, 3> origin{0, 0, 0};
    const complex_t c{0.7, 0.0};
    const std::array<complex_t, mesh> wc_q{{
        {1.2, -0.3},
        {-0.4, 0.7},
        {0.9, 0.2},
    }};
    const std::array<complex_t, mesh> weighted_green_k{{
        {0.13, 0.04},
        {-0.08, 0.02},
        {0.05, -0.06},
    }};

    tensor_t c_tensor({1, 1, 1});
    c_tensor(0, 0, 0) = c;
    tensor_map_t cs{{0, {{{0, origin}, c_tensor}}}};
    tensor_map_t ws;
    tensor_map_t gs;
    for (int iR = 0; iR != mesh; ++iR)
    {
        complex_t wc_R{};
        complex_t green_R{};
        for (int iq = 0; iq != mesh; ++iq)
        {
            const double angle = -2.0 * std::acos(-1.0) * static_cast<double>(iq * iR) / mesh;
            const complex_t phase{std::cos(angle), std::sin(angle)};
            wc_R += wc_q[iq] * phase / static_cast<double>(mesh);
            green_R += weighted_green_k[iq] * phase;
        }
        tensor_t w_tensor({1, 1});
        tensor_t g_tensor({1, 1});
        w_tensor(0, 0) = wc_R;
        g_tensor(0, 0) = green_R;
        const std::array<int, 3> r{iR, 0, 0};
        ws[0][{0, r}] = w_tensor;
        gs[0][{0, r}] = g_tensor;
    }

    RI::GW<int, int, 3, complex_t> gw;
    const std::array<std::array<double, 3>, 3> lattice{{
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
    }};
    gw.set_parallel(MPI_COMM_WORLD, {{0, {0.0, 0.0, 0.0}}}, lattice, {mesh, 1, 1});
    gw.set_symmetry(false, {});
    gw.set_Cs(cs, 0.0);
    gw.set_Ws(ws, 0.0);
    gw.set_Gs(gs, 0.0);
    gw.cal_Sigmas();

    for (int ik = 0; ik != mesh; ++ik)
    {
        complex_t sigma_from_r{};
        for (const auto &[jr, sigma] : gw.Sigmas.at(0))
        {
            assert(jr.first == 0);
            const int iR = jr.second[0];
            const double angle = 2.0 * std::acos(-1.0) * static_cast<double>(ik * iR) / mesh;
            sigma_from_r += sigma(0, 0) * complex_t{std::cos(angle), std::sin(angle)};
        }

        complex_t direct{};
        for (int iq = 0; iq != mesh; ++iq)
            direct += wc_q[iq] * weighted_green_k[(ik - iq + mesh) % mesh];
        direct *= 4.0 * c * c;
        assert(std::abs(sigma_from_r - direct) < 1e-12);
    }
}

void test_libri_gw_two_atom_cross_block_translation()
{
    using complex_t = std::complex<double>;
    using tensor_t = RI::Tensor<complex_t>;
    using tensor_map_t =
        std::map<int, std::map<std::pair<int, std::array<int, 3>>, tensor_t>>;

    const std::array<int, 3> origin{0, 0, 0};
    const std::array<int, 3> r1{1, 0, 0};
    const complex_t c0{0.8, 0.0};
    const complex_t c1{1.1, 0.0};
    const complex_t w01{0.6, -0.15};
    const complex_t g01{-0.35, 0.22};

    tensor_t c0_tensor({1, 1, 1});
    tensor_t c1_tensor({1, 1, 1});
    tensor_t w_tensor({1, 1});
    tensor_t g_tensor({1, 1});
    c0_tensor(0, 0, 0) = c0;
    c1_tensor(0, 0, 0) = c1;
    w_tensor(0, 0) = w01;
    g_tensor(0, 0) = g01;

    tensor_map_t cs{
        {0, {{{0, origin}, c0_tensor}}},
        {1, {{{1, origin}, c1_tensor}}},
    };
    tensor_map_t ws{{0, {{{1, r1}, w_tensor}}}};
    tensor_map_t gs{{0, {{{1, r1}, g_tensor}}}};

    RI::GW<int, int, 3, complex_t> gw;
    const std::array<std::array<double, 3>, 3> lattice{{
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
    }};
    gw.set_parallel(MPI_COMM_WORLD,
                    {{0, {0.0, 0.0, 0.0}}, {1, {0.25, 0.0, 0.0}}},
                    lattice, {3, 1, 1});
    gw.set_symmetry(false, {});
    gw.set_Cs(cs, 0.0);
    gw.set_Ws(ws, 0.0);
    gw.set_Gs(gs, 0.0);
    gw.cal_Sigmas();

    const auto &sigma = gw.Sigmas.at(0).at({1, r1});
    const complex_t expected = 4.0 * c0 * c1 * w01 * g01;
    assert(std::abs(sigma(0, 0) - expected) < 1e-12);
}

void test_libri_gw_two_atom_cell_gauge_covariance()
{
    using complex_t = std::complex<double>;
    using tensor_t = RI::Tensor<complex_t>;
    using tensor_map_t =
        std::map<int, std::map<std::pair<int, std::array<int, 3>>, tensor_t>>;

    constexpr int period_x = 3;
    const std::array<int, 3> origin{0, 0, 0};
    const std::array<int, 3> r1{1, 0, 0};
    const std::array<int, 3> r2{2, 0, 0};

    auto scalar_tensor = [](const complex_t value, const bool coefficient) {
        if (coefficient)
        {
            tensor_t tensor({1, 1, 1});
            tensor(0, 0, 0) = value;
            return tensor;
        }
        else
        {
            tensor_t tensor({1, 1});
            tensor(0, 0) = value;
            return tensor;
        }
    };

    tensor_map_t cs{
        {0, {{{0, origin}, scalar_tensor({0.70, 0.00}, true)},
             {{1, r1}, scalar_tensor({0.25, -0.10}, true)}}},
        {1, {{{0, r2}, scalar_tensor({-0.15, 0.08}, true)},
             {{1, origin}, scalar_tensor({0.90, 0.00}, true)}}},
    };
    tensor_map_t ws{
        {0, {{{0, origin}, scalar_tensor({0.80, -0.12}, false)},
             {{1, r1}, scalar_tensor({0.35, 0.21}, false)}}},
        {1, {{{0, r2}, scalar_tensor({-0.18, 0.17}, false)},
             {{1, origin}, scalar_tensor({0.55, -0.09}, false)}}},
    };
    tensor_map_t gs{
        {0, {{{0, origin}, scalar_tensor({-0.31, 0.14}, false)},
             {{1, r1}, scalar_tensor({0.22, -0.19}, false)}}},
        {1, {{{0, r2}, scalar_tensor({0.11, 0.16}, false)},
             {{1, origin}, scalar_tensor({-0.27, -0.05}, false)}}},
    };

    const std::array<std::array<double, 3>, 3> lattice{{
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
    }};
    auto contract = [&](const tensor_map_t &coefficients, const tensor_map_t &screened,
                        const tensor_map_t &green) {
        RI::GW<int, int, 3, complex_t> gw;
        gw.set_parallel(MPI_COMM_WORLD,
                        {{0, {0.0, 0.0, 0.0}}, {1, {0.25, 0.0, 0.0}}},
                        lattice, {period_x, 1, 1});
        gw.set_symmetry(false, {});
        gw.set_Cs(coefficients, 0.0);
        gw.set_Ws(screened, 0.0);
        gw.set_Gs(green, 0.0);
        gw.cal_Sigmas();
        return gw.Sigmas;
    };

    const std::array<int, 2> atom_cell_shift{0, 1};
    auto gauge_transform = [&](const tensor_map_t &input) {
        tensor_map_t output;
        for (const auto &[atom_i, blocks] : input)
        {
            for (const auto &[jr, tensor] : blocks)
            {
                auto shifted_r = jr.second;
                shifted_r[0] = (shifted_r[0] + atom_cell_shift[jr.first] -
                                atom_cell_shift[atom_i]) % period_x;
                if (shifted_r[0] < 0) shifted_r[0] += period_x;
                output[atom_i][{jr.first, shifted_r}] = tensor;
            }
        }
        return output;
    };
    auto normalize_cell_labels = [&](const tensor_map_t &input) {
        tensor_map_t output;
        for (const auto &[atom_i, blocks] : input)
        {
            for (const auto &[jr, tensor] : blocks)
            {
                auto normalized_r = jr.second;
                normalized_r[0] %= period_x;
                if (normalized_r[0] < 0) normalized_r[0] += period_x;
                output[atom_i][{jr.first, normalized_r}] = tensor;
            }
        }
        return output;
    };

    const auto sigma = contract(cs, ws, gs);
    const auto shifted_sigma = normalize_cell_labels(
        contract(gauge_transform(cs), gauge_transform(ws), gauge_transform(gs)));
    const auto expected_shifted_sigma =
        normalize_cell_labels(gauge_transform(sigma));

    assert(shifted_sigma.size() == expected_shifted_sigma.size());
    for (const auto &[atom_i, expected_blocks] : expected_shifted_sigma)
    {
        const auto actual_i = shifted_sigma.find(atom_i);
        assert(actual_i != shifted_sigma.end());
        assert(actual_i->second.size() == expected_blocks.size());
        for (const auto &[jr, expected] : expected_blocks)
        {
            const auto actual = actual_i->second.find(jr);
            assert(actual != actual_i->second.end());
            assert(actual->second.shape.size() == expected.shape.size());
            for (std::size_t idim = 0; idim < expected.shape.size(); ++idim)
                assert(actual->second.shape[idim] == expected.shape[idim]);
            assert(std::abs(actual->second(0, 0) - expected(0, 0)) < 1e-12);
        }
    }
}
#endif

int main(int argc, char **argv)
{
    int mpi_initialized = 0;
    MPI_Initialized(&mpi_initialized);
    if (!mpi_initialized) MPI_Init(&argc, &argv);

    test_comm_map2_targets_from_vector();
    test_comm_map2_targets_from_set();
    test_comm_map2_targets_output_types();
    test_comm_map2_targets_empty_cells();

#ifdef LIBRPA_USE_LIBRI
    test_libri_gw_scalar_complex_contraction();
    test_libri_gw_scalar_translation_contraction();
    test_libri_gw_scalar_mesh_matches_k_minus_q_convolution();
    test_libri_gw_two_atom_cross_block_translation();
    test_libri_gw_two_atom_cell_gauge_covariance();
#endif

    int mpi_finalized = 0;
    MPI_Finalized(&mpi_finalized);
    if (!mpi_finalized) MPI_Finalize();

    return 0;
}
