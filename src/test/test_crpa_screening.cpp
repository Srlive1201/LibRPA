#include <omp.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <valarray>

#include "../core/crpa.h"
#include "../core/crpa_workflow.h"
#include "../core/dataset.h"
#include "../io/global_io.h"
#include "../math/utils_matrix_m_mpi.h"
#include "../mpi/global_mpi.h"

using namespace librpa_int;

namespace
{
using Complex = std::complex<double>;

void close(Complex actual, Complex expected, double tolerance)
{
    if (!std::isfinite(actual.real()) || !std::isfinite(actual.imag()) ||
        std::abs(actual - expected) > tolerance)
    {
        std::cerr << std::setprecision(17) << "screening actual=" << actual
                  << " expected=" << expected << ", tolerance=" << tolerance << '\n';
        throw std::runtime_error("cRPA screening numerical mismatch");
    }
}

void test_distributed_complex_contraction()
{
    Dataset ds{MPI_COMM_WORLD};
    constexpr int naux = 5, nleft = 2, nright = 3;
    ArrayDesc descriptor(ds.blacs_h);
    descriptor.init_square_blk(naux, naux, 0, 0);
    auto local_kernel = init_local_mat<Complex>(descriptor, MAJOR::COL);
    ComplexMatrix kernel(naux, naux), left(nleft * nleft, naux), right(nright * nright, naux);
    for (int i = 0; i < naux; ++i)
        for (int j = 0; j < naux; ++j)
        {
            kernel(i, j) = Complex(0.12 * (i + 1) * (j + 1) + (i == j ? 2.0 : 0.0), 0.17 * (i - j));
            const int row = descriptor.indx_g2l_r(i), column = descriptor.indx_g2l_c(j);
            if (row >= 0 && column >= 0) local_kernel(row, column) = kernel(i, j);
        }
    for (int pair = 0; pair < left.nr; ++pair)
        for (int i = 0; i < naux; ++i)
            left(pair, i) = Complex(0.13 * (pair + 1) + 0.07 * i, 0.09 * (pair - 2 * i));
    for (int pair = 0; pair < right.nr; ++pair)
        for (int i = 0; i < naux; ++i)
            right(pair, i) = Complex(0.11 * (pair + 2) - 0.06 * i, 0.08 * (2 * pair + i));
    const auto actual =
        contract_crpa_blacs(left, right, local_kernel, descriptor, nleft, nright, ds.comm_h);
    for (int a = 0; a < nleft; ++a)
        for (int b = 0; b < nleft; ++b)
            for (int c = 0; c < nright; ++c)
                for (int d = 0; d < nright; ++d)
                {
                    Complex expected = 0;
                    for (int i = 0; i < naux; ++i)
                        for (int j = 0; j < naux; ++j)
                            expected += std::conj(left(b * nleft + a, i)) * kernel(i, j) *
                                        right(c * nright + d, j);
                    close(actual(a * nleft + b, c * nright + d), expected, 2e-12);
                }
}

void test_native_screening_pipeline(bool negative_mode, bool spin_polarized = false)
{
    Dataset ds{MPI_COMM_WORLD};
    const std::size_t naux = negative_mode ? 2 : 1;
    ds.basis_wfc.set(std::vector<std::size_t>{2});
    ds.basis_aux.set(std::vector<std::size_t>{naux});
    const int nspins = spin_polarized ? 2 : 1;
    ds.mf.set(nspins, 4, 2, 2, 1);
    ds.mf.get_efermi() = 0;
    ds.pbc.set_latvec({1, 0, 0, 0, 1, 0, 0, 0, 1});
    const double pi = std::acos(-1.0);
    ds.pbc.set_kgrids_kvec(4, 1, 1, {0, 0, 0, pi / 2, 0, 0, pi, 0, 0, -pi / 2, 0, 0});
    ds.scfk_blacs_ctxt.init(KPointBlacsProcessShape(1, ds.comm_h.nprocs, true), MPI_COMM_WORLD, 4);
    ds.desc_wfc_kb_full = ds.scfk_blacs_ctxt.create_array_desc(2, 2, 2, 2);
    ds.desc_abf.reset_handler(ds.blacs_h);
    ds.desc_abf.init_1b1p(naux, naux, 0, 0);
    ds.tfg.reset(16);
    ds.tfg.generate_minimax(0.4, 4.0);

    SiteOrbitalGroup site;
    site.label = "mixed-complex-orbital";
    site.atom_index = 0;
    site.orb_start = 0;
    site.n_orbitals = 1;
    CorrelatedSubspace frame({site}, ds.pbc.kfrac_list, ds.pbc.Rlist, 2, 2, 1);
    std::unique_ptr<CorrelatedSubspace> down_frame;
    if (spin_polarized)
        down_frame = std::make_unique<CorrelatedSubspace>(std::vector<SiteOrbitalGroup>{site},
                                                          ds.pbc.kfrac_list, ds.pbc.Rlist, 2, 2, 1);
    constexpr double down_angle = 0.31, down_gap = 1.6;
    for (int k = 0; k < 4; ++k)
    {
        ComplexMatrix eigenvectors(2, 2), overlap(2, 2), trial(2, 1);
        eigenvectors.set_as_identity_matrix();
        overlap.set_as_identity_matrix();
        trial(0, 0) = 1.0 / std::sqrt(2.0);
        trial(1, 0) = Complex(0.0, 1.0 / std::sqrt(2.0));
        ds.mf.get_eigenvectors()[0][0][k] = eigenvectors;
        ds.mf.get_eigenvals()[0](k, 0) = -0.5;
        ds.mf.get_eigenvals()[0](k, 1) = 0.5;
        ds.mf.get_weight()[0](k, 0) = spin_polarized ? 0.25 : 0.5;
        ds.mf.get_weight()[0](k, 1) = 0;
        frame.set_S_k(k, overlap);
        frame.set_W_k(k, trial);
        frame.build_spin_k(0, k, eigenvectors);
        if (spin_polarized)
        {
            ds.mf.get_eigenvectors()[1][0][k] = eigenvectors;
            ds.mf.get_eigenvals()[1](k, 0) = -down_gap / 2;
            ds.mf.get_eigenvals()[1](k, 1) = down_gap / 2;
            ds.mf.get_weight()[1](k, 0) = 0.25;
            ds.mf.get_weight()[1](k, 1) = 0;
            trial(0, 0) = std::cos(down_angle);
            trial(1, 0) = Complex(0, -std::sin(down_angle));
            down_frame->set_S_k(k, overlap);
            down_frame->set_W_k(k, trial);
            down_frame->build_spin_k(0, k, eigenvectors);
        }
        if (ds.comm_h.is_root())
        {
            auto v = std::make_shared<ComplexMatrix>(naux, naux);
            if (negative_mode)
            {
                // Eigenvalues 2 and -0.4 in a rotated auxiliary basis.
                (*v)(0, 0) = (*v)(1, 1) = 0.8;
                (*v)(0, 1) = (*v)(1, 0) = 1.2;
            }
            else
                (*v)(0, 0) = 2.0;
            ds.vq[0][0][ds.pbc.klist[k]] = v;
        }
    }
    frame.compute_T_R();
    if (down_frame) down_frame->compute_T_R();
    // Model normal distributed API input: only block owners set the tag.
    if (ds.comm_h.is_root())
    {
        ds.cs_data.use_libri = true;
        auto data = std::make_shared<std::valarray<double>>(0.0, 4 * naux);
        const double retained[] = {0.15, 0.2, 0.2, 0.4};
        const double discarded[] = {0.3, 0.05, 0.05, 0.1};
        for (int i = 0; i < 4; ++i)
        {
            (*data)[i] =
                negative_mode ? (retained[i] + discarded[i]) / std::sqrt(2.0) : retained[i];
            if (negative_mode) (*data)[4 + i] = (retained[i] - discarded[i]) / std::sqrt(2.0);
        }
        ds.cs_data.data_libri[0][{0, {0, 0, 0}}] = RI::Tensor<double>({naux, 2, 2}, data);
    }
    const auto assigned =
        dispatch_upper_triangular_tasks(1, ds.blacs_h.myid, ds.blacs_h.nprows, ds.blacs_h.npcols,
                                        ds.blacs_h.myprow, ds.blacs_h.mypcol);
    const std::vector<atpair_t> pairs(assigned.begin(), assigned.end());
    Chi0 chi(ds.mf, ds.basis_wfc, ds.basis_aux, ds.pbc, ds.symmetry_context, ds.tfg,
             ds.scfk_blacs_ctxt, ds.desc_wfc_kb_full, false, false);
    chi.gf_threshold = chi.libri_threshold_C = chi.libri_threshold_G = 0;
    std::map<Vector3_Order<double>, ComplexMatrix> no_shrink;
    CrpaContext context{chi,       ds.cs_data, ds.basis_aux, pairs,
                        no_shrink, ds.vq,      ds.blacs_h,   ds.desc_abf};
    std::vector<CorrelatedSubspace*> frames{&frame};
    if (down_frame) frames.push_back(down_frame.get());
    const auto original_mf = ds.mf;
    Chi0::BandSelection selection(nspins, std::vector<std::vector<double>>(4, {1, 0}));
    for (const bool all_bands : {false, true})
    {
        for (auto& spin : selection)
            for (auto& row : spin) row[1] = all_bands;
        const auto result = compute_crpa_onsite(context, selection, frames);
        if (result.tensors.size() != nspins * nspins || result.frequencies.size() != 16)
            throw std::runtime_error("native cRPA pipeline returned incomplete tensors");
        const double density_vertex[] = {0.55, 2 * (0.15 * std::pow(std::cos(down_angle), 2) +
                                                    0.4 * std::pow(std::sin(down_angle), 2))};
        // Match the reductions for the tight shared-response comparison.
        // Keep result above threaded for the independent analytic U/W checks.
        const int requested_threads = omp_get_max_threads();
        omp_set_num_threads(1);
        const auto repeated = compute_crpa_onsite(context, selection, frames);
        const auto full = build_crpa_total_response(context);
        auto resident_full = full;
        chi.swap_chi0_q(resident_full);
        auto incoming = selection;
        for (auto& row : incoming[0]) row = {0, 1};
        chi.set_band_selection(incoming);
        if (ds.comm_h.is_root())
            std::cout << "Shared-P0 U-only begin: "
                      << (spin_polarized ? "unequal collinear spins" : "scalar")
                      << " retained Coulomb pipeline" << std::endl;
        const auto shared =
            evaluate_crpa_onsite_from_response(context, selection, frames, chi.get_chi0_q());
        omp_set_num_threads(requested_threads);
        if (ds.comm_h.is_root()) std::cout << "Shared-P0 U-only end\n";
        if (shared.tensors.size() != result.tensors.size() ||
            shared.frequencies != result.frequencies || chi.get_band_selection() != incoming)
            throw std::runtime_error("shared cRPA lost spin frames, nodes or incoming selector");
        const auto& resident = chi.get_chi0_q();
        if (resident.size() != full.size()) throw std::runtime_error("shared cRPA consumed P0");
        for (const auto& [frequency, qs] : full)
        {
            if (resident.at(frequency).size() != qs.size())
                throw std::runtime_error("shared cRPA changed P0 q ownership");
            for (const auto& [q, atoms] : qs)
                for (const auto& [left, rights] : atoms)
                    for (const auto& [right, block] : rights)
                    {
                        const auto& preserved = resident.at(frequency).at(q).at(left).at(right);
                        for (int i = 0; i < block.size; ++i) close(preserved.c[i], block.c[i], 0);
                    }
        }
        double maximum = 0;
        for (std::size_t i = 0; i < result.tensors.size(); ++i)
        {
            const auto& tensor = result.tensors[i];
            const auto& borrowed = shared.tensors[i];
            const double bare =
                2 * density_vertex[tensor.spin_left] * density_vertex[tensor.spin_right];
            close(tensor.bare(0, 0), bare, 1e-12);
            close(borrowed.bare(0, 0), bare, 1e-12);
            if (!borrowed.w.empty() || borrowed.u.size() != result.frequencies.size() ||
                borrowed.spin_left != tensor.spin_left || borrowed.spin_right != tensor.spin_right)
                throw std::runtime_error("shared cRPA returned W or collapsed physical spin pairs");
            for (std::size_t iw = 0; iw < result.frequencies.size(); ++iw)
            {
                const double omega = result.frequencies[iw];
                const double polarizability =
                    spin_polarized ? 0.64 * (1 / (1 + omega * omega) +
                                             down_gap / (down_gap * down_gap + omega * omega))
                                   : 1.28 / (1 + omega * omega);
                const double screened = bare / (1 + polarizability);
                close(tensor.w[iw](0, 0), screened, 2e-5);
                close(tensor.u[iw](0, 0), all_bands ? bare : screened, 2e-5);
                close(borrowed.u[iw](0, 0), all_bands ? bare : screened, 2e-5);
                const auto paired = repeated.tensors.at(i).u.at(iw)(0, 0);
                close(borrowed.u[iw](0, 0), paired, 1e-12);
                maximum = std::max(maximum, std::abs(borrowed.u[iw](0, 0) - paired));
            }
        }
        for (int spin = 0; spin < nspins; ++spin)
            for (int k = 0; k < 4; ++k)
            {
                for (int band = 0; band < 2; ++band)
                {
                    close(ds.mf.get_eigenvals()[spin](k, band),
                          original_mf.get_eigenvals()[spin](k, band), 0);
                    close(ds.mf.get_weight()[spin](k, band),
                          original_mf.get_weight()[spin](k, band), 0);
                }
                const auto& actual = *ds.mf.find_wfc(spin, 0, k);
                const auto& original = *original_mf.find_wfc(spin, 0, k);
                for (int i = 0; i < actual.size; ++i) close(actual.c[i], original.c[i], 0);
            }
        if (ds.comm_h.is_root())
            std::cout << "Native "
                      << (spin_polarized ? "four unequal collinear spin pairs" : "scalar")
                      << " shared-P0 v/U: max=" << maximum << "; P0, selector, E/C/f unchanged\n";
    }
}
}  // namespace

int main(int argc, char** argv)
{
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    global::init_global_mpi(MPI_COMM_WORLD);
    global::init_global_io();
    test_distributed_complex_contraction();
    test_native_screening_pipeline(false);
    test_native_screening_pipeline(true);
    test_native_screening_pipeline(false, true);
    global::finalize_global_io();
    global::finalize_global_mpi();
    MPI_Finalize();
    std::cout << "test_crpa_screening: passed\n";
}
