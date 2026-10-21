#include "optical_task.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <vector>

#include "../mpi/collective_stage.h"
#include "../utils/constants.h"
#include "dataset_helper.h"
#include "epsilon.h"
#include "qsgw/headwing_velocity.h"

namespace librpa_int
{
namespace
{

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

velocity_matrix_t decode_velocity(const double *values, const int ns, const int nk, const int n)
{
    if (!values || ns <= 0 || nk <= 0 || n <= 0)
        throw LIBRPA_INVALID_ARGUMENT("Optical velocity dimensions or values are incomplete");
    velocity_matrix_t velocity;
    initialize_velocity_matrix(velocity, ns, nk, n);
    std::size_t offset = 0;
    for (int spin = 0; spin != ns; ++spin)
        for (int kpoint = 0; kpoint != nk; ++kpoint)
            for (int axis = 0; axis != 3; ++axis)
                for (int row = 0; row != n; ++row)
                    for (int column = 0; column != n; ++column, ++offset)
                        velocity[spin][kpoint][axis](row, column) = {values[2 * offset],
                                                                     values[2 * offset + 1]};
    return velocity;
}

MeanField decode_optical_meanfield(const LibrpaOpticalMeanFieldInput &input,
                                   const bool require_spectrum)
{
    if (input.n_spins <= 0 || input.n_kpoints <= 0 || input.n_states <= 0 || input.n_aos <= 0 ||
        (input.n_spinor != 1 && input.n_spinor != 2) || input.n_wavefunction_blocks < 0 ||
        (input.n_wavefunction_blocks && !input.wavefunction_blocks) ||
        ((input.eigenvalues_ha == nullptr) != (input.occupation_weights == nullptr)) ||
        (require_spectrum && !input.eigenvalues_ha))
        throw LIBRPA_INVALID_ARGUMENT("Optical mean-field dimensions or arrays are incomplete");
    MeanField mf(input.n_spins, input.n_kpoints, input.n_states, input.n_aos, input.n_spinor);
    mf.get_efermi() = input.fermi_energy_ha;
    if (input.eigenvalues_ha)
        for (int spin = 0; spin != input.n_spins; ++spin)
            for (int kpoint = 0; kpoint != input.n_kpoints; ++kpoint)
                for (int band = 0; band != input.n_states; ++band)
                {
                    const auto index = (static_cast<std::size_t>(spin) * input.n_kpoints + kpoint) *
                                           input.n_states +
                                       band;
                    mf.get_eigenvals()[spin](kpoint, band) = input.eigenvalues_ha[index];
                    mf.get_weight()[spin](kpoint, band) = input.occupation_weights[index];
                }
    for (int i = 0; i != input.n_wavefunction_blocks; ++i)
    {
        const auto &block = input.wavefunction_blocks[i];
        if (block.spin < 0 || block.spin >= input.n_spins || block.spinor < 0 ||
            block.spinor >= input.n_spinor || block.kpoint < 0 || block.kpoint >= input.n_kpoints ||
            block.rows < 0 || block.rows > input.n_states || block.columns < 0 ||
            block.columns > input.n_aos || (block.rows && block.columns && !block.values_ri))
            throw LIBRPA_INVALID_ARGUMENT("Optical wavefunction block is inconsistent");
        if (!block.rows || !block.columns) continue;
        ComplexMatrix values(block.rows, block.columns);
        for (std::size_t index = 0; index != static_cast<std::size_t>(values.size); ++index)
            values.c[index] = {block.values_ri[2 * index], block.values_ri[2 * index + 1]};
        if (!mf.get_eigenvectors()[block.spin][block.spinor]
                 .emplace(block.kpoint, std::move(values))
                 .second)
            throw LIBRPA_INVALID_ARGUMENT("Optical wavefunction block is repeated");
    }
    return mf;
}

}  // namespace

void prepare_native_velocity(Dataset &dataset, const LibrpaOptions *p_opts,
                             const LibrpaNativeVelocityInput *root_input, int prepare_optical)
{
    using namespace librpa_int;
    using namespace librpa_int::qsgw;
    auto *ds = &dataset;
    const auto &comm = ds->comm_h;
    collective_all_stage(
        comm, "Native velocity options",
        [&]
        {
            if (!p_opts)
                throw LIBRPA_INVALID_ARGUMENT("Native velocity requires numerical options");
            if (p_opts->option_dielect_func != 3 && p_opts->option_dielect_func != 4)
                throw LIBRPA_INVALID_ARGUMENT(
                    "Native velocity API options must use numerical values 3/4");
            if (prepare_optical && p_opts->replace_w_head != LIBRPA_SWITCH_ON)
                throw LIBRPA_INVALID_ARGUMENT(
                    "Native optical preparation requires replace_w_head enabled");
        });
    const auto &opts = *p_opts;
    const auto &mf = ds->mf;
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
            {
                const bool use_symmetry = (opts.use_symmetry_rpa == LIBRPA_SWITCH_ON) ||
                                          (opts.use_symmetry_gw == LIBRPA_SWITCH_ON);
                if (!spinor || !use_symmetry || (opts.use_kpara_scf_eigvec == LIBRPA_SWITCH_ON))
                    throw LIBRPA_RUNTIME_ERROR(
                        "IBZ native velocity requires SOC magnetic symmetry and undistributed "
                        "AO/band WFC blocks");
                initialize_symmetry_context(*ds, true);
                require_symmetry_shell_layouts(*ds, "native SOC head/wing");
                const auto &ctx = ds->symmetry_context;
                if (ctx.kstars.size() != static_cast<std::size_t>(nk) ||
                    ctx.count_kstar_members() != static_cast<std::size_t>(nk_full) ||
                    ctx.spin_operations.empty())
                    throw LIBRPA_RUNTIME_ERROR(
                        "IBZ native SOC velocity requires magnetic k-stars covering the full mesh");
            }
            for (int ik = 0; ik != nk; ++ik)
            {
                const auto weight = ds->pbc.weight_k[ik];
                const double expected =
                    reduced_grid
                        ? static_cast<double>(find_symmetry_kstar_for_ibz_kpoint(
                                                  ds->symmetry_context, ds->pbc.kfrac_list[ik])
                                                  .members.size()) /
                              nk_full
                        : 1.0 / nk;
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
                             ds->velocity_matrix = decode_velocity(
                                 reinterpret_cast<const double *>(velocity.data()), ns, nk, n);
                             ds->p_headwing.reset();
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

void prepare_supplied_velocity(Dataset &ds, const LibrpaOptions *options,
                               const LibrpaSuppliedVelocityInput *input, const int need_wing)
{
    const auto &comm = ds.comm_h;
    MeanField supplied_mf, full_source_mf;
    velocity_matrix_t velocity, full_source_velocity;
    std::vector<Vector3_Order<double>> full_source_kfrac;
    std::vector<double> frequencies;
    collective_all_stage(
        comm, "Supplied optical input",
        [&]
        {
            if (!options || !input || options->replace_w_head != LIBRPA_SWITCH_ON ||
                (options->option_dielect_func != 3 && options->option_dielect_func != 4 &&
                 options->option_dielect_func != 53))
                throw LIBRPA_INVALID_ARGUMENT(
                    "Supplied optical input requires numerical options 3/4/53 and head "
                    "replacement");
            if ((input->coefficients_are_kblacs_2d != 0 &&
                 input->coefficients_are_kblacs_2d != 1) ||
                (input->normalize_ibz_occupations != 0 && input->normalize_ibz_occupations != 1) ||
                input->n_frequencies < 0 || (input->n_frequencies && !input->frequencies_ha))
                throw LIBRPA_INVALID_ARGUMENT("Supplied optical conventions are incomplete");
            if (input->meanfield) supplied_mf = decode_optical_meanfield(*input->meanfield, true);
            const auto &mf = input->meanfield ? supplied_mf : ds.mf;
            if (input->n_spins != mf.get_n_spins() || input->n_kpoints != mf.get_n_kpoints() ||
                input->n_states != mf.get_n_states() || mf.get_n_spins() != ds.mf.get_n_spins() ||
                mf.get_n_aos() != ds.mf.get_n_aos() || mf.get_n_spinor() != ds.mf.get_n_spinor() ||
                mf.get_n_kpoints() != static_cast<int>(ds.pbc.kfrac_list.size()))
                throw LIBRPA_INVALID_ARGUMENT(
                    "Supplied optical frame disagrees with the active AO/spin/k grid");
            velocity = decode_velocity(input->velocity_ri, input->n_spins, input->n_kpoints,
                                       input->n_states);
            if (input->n_frequencies)
                frequencies.assign(input->frequencies_ha,
                                   input->frequencies_ha + input->n_frequencies);
            if (input->full_source_meanfield)
            {
                full_source_mf = decode_optical_meanfield(*input->full_source_meanfield, false);
                if (!input->full_source_kfrac || full_source_mf.get_n_spins() != mf.get_n_spins() ||
                    full_source_mf.get_n_spinor() != mf.get_n_spinor() ||
                    full_source_mf.get_n_states() != mf.get_n_states() ||
                    full_source_mf.get_n_aos() != mf.get_n_aos())
                    throw LIBRPA_INVALID_ARGUMENT("Full-source optical frame is inconsistent");
                full_source_velocity =
                    decode_velocity(input->full_source_velocity_ri, mf.get_n_spins(),
                                    full_source_mf.get_n_kpoints(), mf.get_n_states());
                for (int ik = 0; ik != full_source_mf.get_n_kpoints(); ++ik)
                    full_source_kfrac.emplace_back(input->full_source_kfrac[3 * ik],
                                                   input->full_source_kfrac[3 * ik + 1],
                                                   input->full_source_kfrac[3 * ik + 2]);
            }
        });
    const auto &opts = *options;
    const auto &mf = input->meanfield ? supplied_mf : ds.mf;
    const std::array<int, 4> branches{input->meanfield != nullptr,
                                      input->full_source_meanfield != nullptr, input->n_frequencies,
                                      need_wing != 0};
    std::array<int, 4> minimum{}, maximum{};
    comm.allreduce(branches.data(), minimum.data(), 4, MPI_MIN);
    comm.allreduce(branches.data(), maximum.data(), 4, MPI_MAX);
    if (minimum != maximum)
        throw LIBRPA_INVALID_ARGUMENT("Supplied optical layout differs between MPI ranks");
    if (input->meanfield)
        for (int spin = 0; spin != mf.get_n_spins(); ++spin)
            for (int component = 0; component != mf.get_n_spinor(); ++component)
                for (int ik = 0; ik != mf.get_n_kpoints(); ++ik)
                {
                    const auto *block = mf.find_wfc(spin, component, ik);
                    const int have_local = block && block->nr > 0;
                    int have_any = 0;
                    comm.allreduce(&have_local, &have_any, 1, MPI_MAX);
                    if (!have_any)
                        throw LIBRPA_RUNTIME_ERROR("Supplied optical frame is missing a WFC block");
                }
    std::unique_ptr<diele_func> optical;
    collective_all_stage(
        comm, "Supplied optical model",
        [&]
        {
            if (frequencies.empty())
            {
                initialize_ds_tfgrids(ds, opts);
                frequencies = ds.tfg.get_freq_nodes();
            }
            if (frequencies.empty())
                throw LIBRPA_RUNTIME_ERROR("Head/wing frequency grid is empty");
            const bool use_symmetry = opts.use_symmetry_rpa == LIBRPA_SWITCH_ON ||
                                      opts.use_symmetry_gw == LIBRPA_SWITCH_ON;
            if (use_symmetry)
            {
                initialize_symmetry_context(ds, true);
                require_symmetry_shell_layouts(ds, "analytic head/wing");
            }
            const auto &basis_aux =
                opts.use_shrink_abfs == LIBRPA_SWITCH_ON ? ds.basis_aux_shrink : ds.basis_aux;
            if (!basis_aux.initialized())
                throw LIBRPA_RUNTIME_ERROR("Head/wing auxiliary basis is not initialized");
            // The optical object keeps a copy of its eigenframe and a reference to the
            // Dataset velocity; no substitution of the SCF mean field is needed.
            ds.p_headwing.reset();
            ds.velocity_matrix = std::move(velocity);
            const auto &descriptor =
                (input->meanfield ? input->coefficients_are_kblacs_2d : ds.eigvecs_kpara_2d_ready())
                    ? ds.desc_wfc_kb
                    : ds.desc_wfc_kb_full;
            optical = std::make_unique<diele_func>(
                mf, ds.velocity_matrix, ds.pbc.kfrac_list, ds.basis_wfc, basis_aux, frequencies,
                mf.get_n_aos(), mf.get_n_states(), mf.get_n_spins(), basis_aux.nb_total, ds.pbc,
                comm, ds.blacs_h, opts.use_kpara_scf_eigvec == LIBRPA_SWITCH_ON,
                &ds.scfk_blacs_ctxt, &descriptor);
            optical->use_soc = mf.get_n_spinor() > 1;
            optical->debug = global::should_output(LIBRPA_VERBOSE_DEBUG);
            if (opts.use_2d_dielectric == LIBRPA_SWITCH_ON)
            {
                const auto normalization = strict_2d_coulomb_head_normalization(
                    ds.pbc, input->auxiliary_monopole_norm_squared);
                const double area = input->inplane_area_bohr2;
                if (!std::isfinite(area) || area <= 0.0 ||
                    std::abs(normalization.inplane_area_bohr2 - area) >
                        1e-10 * std::max(normalization.inplane_area_bohr2, area))
                    throw LIBRPA_INVALID_ARGUMENT(
                        "Coulomb normalization area disagrees with the lattice");
                optical->configure_strict_2d_coulomb_head(true,
                                                          input->auxiliary_monopole_norm_squared);
            }
            else
                optical->configure_strict_2d_coulomb_head(false);
            if (use_symmetry && ds.symmetry_context.available &&
                !ds.symmetry_context.kstars.empty())
            {
                optical->set_symmetry_context(ds.symmetry_context);
                optical->use_symmetry = true;
                optical->atom_nw = ds.basis_wfc.get_atom_nb_map();
                optical->coord_frac = ds.symmetry_context.input_coord_frac;
            }
            if (full_source_mf.initialized())
            {
                if (!use_symmetry || mf.get_n_kpoints() == ds.pbc.get_n_cells_bvk())
                    throw LIBRPA_INVALID_ARGUMENT(
                        "Full-source optical cache requires an active IBZ grid");
                auto member_indices = map_symmetry_kstar_members_to_source_kpoints(
                    ds.symmetry_context, ds.pbc.kfrac_list, full_source_kfrac, 1e-5);
                optical->set_direct_full_bz_headwing_inputs(std::move(full_source_velocity),
                                                            std::move(full_source_mf),
                                                            std::move(member_indices));
            }
            if (input->normalize_ibz_occupations)
                optical->normalize_ibz_star_weighted_occupations();
        });
    optical->init(opts.sqrt_coulomb_threshold, ds.vq);
    optical->cal_head();
    optical->test_head();
    if (need_wing)
    {
        const auto &coefficients =
            opts.use_shrink_abfs == LIBRPA_SWITCH_ON ? ds.cs_data_shrink : ds.cs_data;
        optical->cal_wing(coefficients, opts.sqrt_coulomb_threshold, ds.vq);
        optical->test_wing();
    }
    ds.epsmacs_imagfreq = optical->get_head_vec();
    ds.omegas_imagfreq = frequencies;
    ds.p_headwing = std::move(optical);
}

}  // namespace librpa_int
