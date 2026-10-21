#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <vector>

#include "../core/epsilon.h"
#include "../core/qsgw/headwing_velocity.h"
#include "../mpi/collective_stage.h"
#include "../utils/constants.h"
#include "dataset_helper.h"
#include "instance_manager.h"
#include "librpa_input.h"
#include "librpa_optical.h"

namespace
{
using namespace librpa_int;

qsgw::LocalizedMatrix decode_operator(const LibrpaLocalizedOperator &input, int dimension)
{
    if (dimension <= 0 || input.dimension != dimension || input.n_blocks <= 0 || !input.blocks)
        throw LIBRPA_INVALID_ARGUMENT("Native operator dimensions or blocks are incomplete");
    qsgw::LocalizedMatrix result;
    for (int block = 0; block < input.n_blocks; ++block)
    {
        const auto &source = input.blocks[block];
        if (!source.values_ri || (source.row_major != 0 && source.row_major != 1))
            throw LIBRPA_INVALID_ARGUMENT("Native operator storage is invalid");
        Matz matrix(dimension, dimension, source.row_major ? MAJOR::ROW : MAJOR::COL);
        for (std::size_t i = 0; i < matrix.size(); ++i)
            matrix.ptr()[i] = {source.values_ri[2 * i], source.values_ri[2 * i + 1]};
        const Vector3_Order<int> translation{source.translation[0], source.translation[1],
                                             source.translation[2]};
        if (!result.emplace(translation, std::move(matrix)).second)
            throw LIBRPA_INVALID_ARGUMENT("Native operator repeats a lattice translation");
    }
    return result;
}

void append_matrix(std::vector<double> &values, const Matz &matrix)
{
    for (int row = 0; row < matrix.nr(); ++row)
        for (int column = 0; column < matrix.nc(); ++column)
        {
            const auto value = matrix(row, column);
            values.push_back(value.real());
            values.push_back(value.imag());
        }
}
}  // namespace

void librpa_prepare_native_velocity(LibrpaHandler *h, const LibrpaOptions *p_opts,
                                    const LibrpaNativeVelocityInput *root_input,
                                    int prepare_optical)
{
    using namespace librpa_int;
    using namespace librpa_int::qsgw;
    const auto ds = api::get_dataset_instance(h);
    if (!p_opts) throw LIBRPA_INVALID_ARGUMENT("native velocity requires numerical options");
    const auto &opts = *p_opts;
    if (opts.option_dielect_func == 6)
        throw LIBRPA_INVALID_ARGUMENT("Native velocity API options must use numerical values 3/4");
    const auto &mf = ds->mf;
    const auto &comm = ds->comm_h;
    const int n = mf.get_n_states(), nk = mf.get_n_kpoints(), ns = mf.get_n_spins();
    const int nao = mf.get_n_aos(), nspinor = mf.get_n_spinor();
    const std::size_t count = static_cast<std::size_t>(ns) * nk * 3 * n * n;
    std::vector<cplxdb> velocity;
    HeadwingAoOperators operators;
    std::vector<LocalizedMatrix> hamiltonians;
    std::array<double, 2> coulomb_normalization{};

    collective_all_stage(
        comm, "G0W0 native velocity input",
        [&]
        {
            const auto period = ds->pbc.period;
            const long long nk_full = static_cast<long long>(period.x) * period.y * period.z;
            const bool collinear = nspinor == 1 && (ns == 1 || ns == 2);
            const bool spinor = nspinor == 2 && ns == 1;
            if ((!collinear && !spinor) || n <= 0 || nao * nspinor != n)
                throw LIBRPA_RUNTIME_ERROR(
                    "G0W0 native velocity requires one/two collinear spins or one two-component "
                    "spinor channel and complete AO-spinor/band frames");
            if (nk <= 0 || static_cast<int>(ds->pbc.kfrac_list.size()) != nk ||
                static_cast<int>(ds->pbc.weight_k.size()) != nk)
                throw LIBRPA_RUNTIME_ERROR(
                    "G0W0 native velocity requires matching SCF k points and weights");
            const bool reduced_grid = nk != nk_full;
            if (reduced_grid)
                throw LIBRPA_RUNTIME_ERROR("Native velocity requires a complete full-BZ grid");
            for (int ik = 0; ik != nk; ++ik)
            {
                const auto weight = ds->pbc.weight_k[ik];
                const double expected = 1.0 / nk;
                if (!std::isfinite(weight) || std::abs(weight - expected) > 1e-8)
                    throw LIBRPA_RUNTIME_ERROR(
                        "G0W0 native velocity k weights disagree with the full-mesh star counts");
            }
            if (count > static_cast<std::size_t>(std::numeric_limits<int>::max()))
                throw LIBRPA_RUNTIME_ERROR("Native velocity exceeds the MPI buffer-count limit");
            velocity.resize(count);
            if (!comm.is_root()) return;
            if (!root_input || root_input->n_spins != ns || !root_input->hamiltonian ||
                root_input->overlap.dimension != n)
                throw LIBRPA_INVALID_ARGUMENT(
                    "Native operators do not match the complete SCF frame");
            operators.dimension = n;
            operators.overlap = decode_operator(root_input->overlap, n);
            for (int axis = 0; axis < 3; ++axis)
                operators.position[axis] = decode_operator(root_input->position[axis], n);
            for (int spin = 0; spin < ns; ++spin)
                hamiltonians.push_back(decode_operator(root_input->hamiltonian[spin], n));
            coulomb_normalization = {root_input->inplane_area_bohr2,
                                     root_input->auxiliary_monopole_norm_squared};
        });

    std::array<double, 4> residuals{};
    for (int spin = 0; spin < ns; ++spin)
        for (int ik = 0; ik < nk; ++ik)
        {
            Matz ket, component;
            collective_all_stage(comm, "G0W0 native eigenframe shape",
                                 [&]
                                 {
                                     ket = Matz(n, n, MAJOR::ROW);
                                     component = Matz(nao, n, MAJOR::ROW);
                                 });
            for (int spinor = 0; spinor < nspinor; ++spinor)
            {
                const auto *wfc = mf.find_wfc(spin, spinor, ik);
                collective_all_stage(
                    comm, "G0W0 native eigenframe component shape",
                    [&]
                    {
                        if (wfc != nullptr && wfc->size != 0 && (wfc->nr != n || wfc->nc != nao))
                            throw LIBRPA_RUNTIME_ERROR(
                                "G0W0 native velocity found a partial SCF spinor block");
                    });
                const int owner = find_mpi_owner_rank(wfc != nullptr && wfc->size != 0, comm.comm);
                if (owner < 0)
                    throw LIBRPA_RUNTIME_ERROR(
                        "G0W0 native velocity is missing a SCF spinor block");
                collective_all_stage(comm, "G0W0 native eigenframe transfer",
                                     [&]
                                     {
                                         if (comm.myid != owner) return;
                                         // MeanField has band rows, spatial AO columns. Native
                                         // H/S/r use AO-major, spinor-minor ket rows.
                                         for (int ao = 0; ao < nao; ++ao)
                                             for (int band = 0; band < n; ++band)
                                                 component(ao, band) = (*wfc)(band, ao);
                                     });
                comm.bcast(component.ptr(), nao * n, owner);
                for (int ao = 0; ao < nao; ++ao)
                    for (int band = 0; band < n; ++band)
                        ket(ao * nspinor + spinor, band) = component(ao, band);
            }
            collective_all_stage(
                comm, "G0W0 native velocity frame",
                [&]
                {
                    if (!comm.is_root()) return;
                    std::vector<double> eigenvalues(n);
                    for (int band = 0; band < n; ++band)
                        eigenvalues[band] = mf.get_eigenvals().at(spin)(ik, band);
                    const auto result = compute_localized_velocity_in_frame(
                        hamiltonians.at(spin), operators, ds->pbc.kfrac_list.at(ik), ds->pbc.latvec,
                        eigenvalues, ket);
                    residuals[0] = std::max(residuals[0], result.metric_residual);
                    residuals[1] = std::max(residuals[1], result.hamiltonian_frame_residual_ha);
                    residuals[2] =
                        std::max(residuals[2], result.hamiltonian_hermiticity_residual_ha);
                    residuals[3] =
                        std::max(residuals[3], result.velocity_upper_diagonal_imaginary_ha_bohr);
                    if (result.metric_residual > 1e-6 ||
                        result.hamiltonian_frame_residual_ha > 1e-4 / HA2EV ||
                        result.hamiltonian_hermiticity_residual_ha > 1e-4 / HA2EV ||
                        result.velocity_upper_diagonal_imaginary_ha_bohr > 1e-8)
                    {
                        std::ostringstream message;
                        message << std::scientific << "spin=" << spin << " k=" << ik
                                << " metric=" << result.metric_residual
                                << " frame_Ha=" << result.hamiltonian_frame_residual_ha
                                << " hermiticity_Ha=" << result.hamiltonian_hermiticity_residual_ha
                                << " diagonal_imag_HaBohr="
                                << result.velocity_upper_diagonal_imaginary_ha_bohr;
                        throw LIBRPA_RUNTIME_ERROR(message.str());
                    }
                    for (int alpha = 0; alpha < 3; ++alpha)
                        for (int i = 0; i < n; ++i)
                            for (int j = 0; j < n; ++j)
                            {
                                const std::size_t offset =
                                    ((((static_cast<std::size_t>(spin) * nk + ik) * 3 + alpha) * n +
                                      i) *
                                         n +
                                     j);
                                velocity[offset] = result.velocity_ha_bohr[alpha](i, j);
                            }
                });
        }
    comm.bcast(velocity.data(), static_cast<int>(velocity.size()), 0);
    collective_all_stage(comm, "G0W0 native velocity setter",
                         [&]
                         {
                             librpa_set_velocity_matrix_packed(
                                 h, ns, nk, n, reinterpret_cast<const double *>(velocity.data()));
                         });
    if (comm.is_root())
        global::lib_printf(
            "G0W0 native velocity: metric %.8e frame_Ha %.8e hermiticity_Ha %.8e "
            "diagonal_imag_HaBohr %.8e; SCF E/C/occupations retained\n",
            residuals[0], residuals[1], residuals[2], residuals[3]);

    if (!prepare_optical) return;

    // Preserve the numerical grid/head/wing lifecycle behind the public API.
    initialize_ds_global_ddla(*ds, opts);
    initialize_ds_tfgrids(*ds, opts);
    initialize_ds_headwing(*ds, opts, opts.option_dielect_func == 3);
    comm.bcast(coulomb_normalization.data(), 2, 0);
    collective_all_stage(
        comm, "Native head/wing Coulomb normalization",
        [&]
        {
            if (opts.use_2d_dielectric == LIBRPA_SWITCH_ON)
            {
                const auto normalization =
                    strict_2d_coulomb_head_normalization(ds->pbc, coulomb_normalization[1]);
                const double area = coulomb_normalization[0];
                if (!std::isfinite(area) || area <= 0.0 ||
                    std::abs(normalization.inplane_area_bohr2 - area) >
                        1.0e-10 * std::max(normalization.inplane_area_bohr2, area))
                    throw LIBRPA_INVALID_ARGUMENT(
                        "Coulomb normalization area disagrees with the lattice");
                ds->p_headwing->configure_strict_2d_coulomb_head(true, coulomb_normalization[1]);
            }
            else
                ds->p_headwing->configure_strict_2d_coulomb_head(false);
        });
}

struct LibrpaOpticalResult
{
    LibrpaOpticalView view{};
    std::vector<double> velocity, frequencies, head, wing;
};

LibrpaOpticalResult *librpa_get_optical_result(const LibrpaHandler *h)
{
    const auto ds = librpa_int::api::get_dataset_instance(h);
    if (!ds->p_headwing) throw LIBRPA_RUNTIME_ERROR("Analytic head/wing is not initialized");
    auto result = std::make_unique<LibrpaOpticalResult>();
    auto &view = result->view;
    view.n_spins = ds->mf.get_n_spins();
    view.n_kpoints = ds->mf.get_n_kpoints();
    view.n_states = ds->mf.get_n_states();
    for (const auto &spin : ds->velocity_matrix)
        for (const auto &k : spin)
            for (const auto &axis : k)
                for (int row = 0; row < axis.nr; ++row)
                    for (int column = 0; column < axis.nc; ++column)
                    {
                        result->velocity.push_back(axis(row, column).real());
                        result->velocity.push_back(axis(row, column).imag());
                    }
    result->frequencies = ds->p_headwing->get_frequency_nodes();
    view.n_frequencies = result->frequencies.size();
    const auto &heads = ds->p_headwing->get_head_tensors();
    const auto &wings = ds->p_headwing->get_wing_auxiliary();
    if (!heads.empty())
    {
        view.head_rows = heads.front().nr();
        view.head_columns = heads.front().nc();
        for (const auto &matrix : heads) append_matrix(result->head, matrix);
    }
    if (!wings.empty())
    {
        view.wing_rows = wings.front().nr();
        view.wing_columns = wings.front().nc();
        for (const auto &matrix : wings) append_matrix(result->wing, matrix);
    }
    view.velocity_ri = result->velocity.data();
    view.frequencies_ha = result->frequencies.data();
    view.head_ri = result->head.data();
    view.wing_ri = result->wing.data();
    return result.release();
}

LibrpaOpticalView librpa_get_optical_view(const LibrpaOpticalResult *result)
{
    if (!result) throw LIBRPA_INVALID_ARGUMENT("Optical result must not be null");
    return result->view;
}

void librpa_delete_optical_result(LibrpaOpticalResult *result) { delete result; }
