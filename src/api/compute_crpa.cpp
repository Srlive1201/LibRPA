#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <string>

#include "../core/crpa.h"
#include "../core/crpa_continuation.h"
#include "../core/crpa_model.h"
#include "../core/crpa_workflow.h"
#include "../utils/constants.h"
#include "crpa_helper.h"
#include "dataset_helper.h"
#include "instance_manager.h"
#include "librpa_crpa.h"

struct LibrpaCrpaResult
{
    librpa_int::CrpaResult values;
    std::vector<std::vector<std::vector<int>>> response_bands;
    std::vector<std::vector<std::vector<int>>> orbital_bands;
    librpa_int::Chi0::BandSelection response_weights;
    librpa_int::CrpaModelSnapshot model;
    std::vector<double> lattice, positions, k_fractional, k_weights;
    std::vector<int> atom_types, ao_counts;
    std::vector<int> shell_offsets{0}, shell_l, basis_convention;
};

struct LibrpaCrpaContinuedResult
{
    std::vector<librpa_int::CrpaContinuedTensor> tensors;
};

namespace
{
void collective_error(const librpa_int::MpiCommHandler& comm, const std::string& error)
{
    const int failed = !error.empty();
    int any_failed = 0;
    comm.allreduce(&failed, &any_failed, 1, MPI_MAX);
    if (any_failed)
        throw LIBRPA_RUNTIME_ERROR("cRPA input: " +
                                   (error.empty() ? "another MPI rank failed validation" : error));
}
}  // namespace

void librpa_init_crpa_input(LibrpaCrpaInput* input)
{
    if (!input) throw LIBRPA_INVALID_ARGUMENT("null cRPA input");
    *input = {};
    input->parent_orbitals = "d";
    input->output_orbitals = "d";
    input->gram_abs = input->gram_rel = input->s_rel = 1.0e-12;
    input->gram_cond_max = 1.0e12;
    input->s_abs = input->residual_tol = input->t_roundtrip_tol = 1.0e-10;
    input->s_cond_max = 1.0e10;
    input->tracking_min_capture = 0.99;
}

void librpa_get_crpa_kgrid(const LibrpaHandler* h, int n_kpoints, double* fractional_k)
{
    const auto ds = librpa_int::api::get_dataset_instance(h);
    if (!fractional_k || n_kpoints <= 0 ||
        ds->pbc.kfrac_list.size() != static_cast<std::size_t>(n_kpoints))
        throw LIBRPA_INVALID_ARGUMENT("cRPA k-grid query has inconsistent dimensions");
    for (int k = 0; k < n_kpoints; ++k)
    {
        fractional_k[3 * k] = ds->pbc.kfrac_list[k].x;
        fractional_k[3 * k + 1] = ds->pbc.kfrac_list[k].y;
        fractional_k[3 * k + 2] = ds->pbc.kfrac_list[k].z;
    }
}

librpa_int::CrpaPreparedModel librpa_int::api::prepare_crpa_api_model(
    Dataset& dataset, const LibrpaOptions& numerical_options, const LibrpaCrpaInput* input,
    CrpaApiMode mode)
{
    using namespace librpa_int;
    if (numerical_options.use_symmetry_rpa == LIBRPA_SWITCH_ON ||
        numerical_options.use_symmetry_gw == LIBRPA_SWITCH_ON ||
        numerical_options.use_symmetry_exx == LIBRPA_SWITCH_ON)
        throw LIBRPA_INVALID_ARGUMENT(
            "cRPA requires full-grid inputs with symmetry reduction disabled");
    auto* ds = &dataset;
    const auto* options = &numerical_options;
    CrpaWindowOptions window;
    std::vector<std::string> species;
    std::vector<ComplexMatrix> overlap;
    std::string error;
    try
    {
        if (!input) throw LIBRPA_INVALID_ARGUMENT("missing cRPA input or options");
        if (options->use_kpara_scf_eigvec == LIBRPA_SWITCH_ON ||
            (options->parallel_routing != LIBRPA_ROUTING_LIBRI &&
             options->parallel_routing != LIBRPA_ROUTING_AUTO))
            throw LIBRPA_INVALID_ARGUMENT(
                "cRPA requires full-grid LIBRI and replicated eigenvectors");
        if (mode == CrpaApiMode::SharedGw && options->read_sigc_mat_rf == LIBRPA_SWITCH_ON)
            throw LIBRPA_INVALID_ARGUMENT(
                "shared GW/cRPA requires a resident response and no SigC file-read shortcut");
        const bool optical = options->replace_w_head == LIBRPA_SWITCH_ON;
        if ((optical && (options->option_dielect_func != 3 ||
                         options->use_2d_dielectric == LIBRPA_SWITCH_ON)) ||
            (!optical && options->option_dielect_func != 0) ||
            options->tfgrids_type != LIBRPA_TFGRID_MINIMAX || options->n_bands_chi0 >= 0 ||
            options->n_bands_exclude_chi0 != 0 ||
            (mode == CrpaApiMode::Standalone && (options->use_fullcoul_eps != LIBRPA_SWITCH_ON ||
                                                 options->use_fullcoul_wc != LIBRPA_SWITCH_ON)))
            throw LIBRPA_INVALID_ARGUMENT(
                "cRPA requires native minimax, all loaded bands, full/full Coulomb and either "
                "no head correction or 3D analytic option 3");
        const bool shrink = options->use_shrink_abfs == LIBRPA_SWITCH_ON;
        if (shrink && (options->use_shrink_chi != LIBRPA_SWITCH_ON ||
                       !ds->basis_aux_shrink.initialized() || ds->sinvS.empty()))
            throw LIBRPA_INVALID_ARGUMENT(
                "cRPA compressed basis requires use_shrink_chi and loaded auxiliary transforms");
        if (!shrink && options->use_shrink_chi == LIBRPA_SWITCH_ON)
            throw LIBRPA_INVALID_ARGUMENT("use_shrink_chi requires use_shrink_abfs for cRPA");
        if (!std::isfinite(options->sqrt_coulomb_threshold) ||
            options->sqrt_coulomb_threshold < 0.0)
            throw LIBRPA_INVALID_ARGUMENT("cRPA Coulomb threshold must be finite and nonnegative");
        if (input->n_species <= 0 || !input->species_labels || !input->correlated_species ||
            !input->parent_orbitals || !input->output_orbitals || !input->overlap_k_ri ||
            input->n_kpoints != ds->mf.get_n_kpoints() || input->n_aos != ds->mf.get_n_aos())
            throw LIBRPA_INVALID_ARGUMENT("incomplete cRPA species/window/overlap input");
        for (int type = 0; type < input->n_species; ++type)
        {
            if (!input->species_labels[type]) throw LIBRPA_INVALID_ARGUMENT("null species label");
            species.emplace_back(input->species_labels[type]);
        }
        window.correlated_species = input->correlated_species;
        window.ligand_species = input->ligand_species ? input->ligand_species : "";
        window.parent_orbitals = input->parent_orbitals;
        window.output_orbitals = input->output_orbitals;
        if (input->use_orbital_weights != LIBRPA_SWITCH_OFF &&
            input->use_orbital_weights != LIBRPA_SWITCH_ON)
            throw LIBRPA_INVALID_ARGUMENT("invalid cRPA orbital-weight switch");
        window.use_orbital_weights = input->use_orbital_weights == LIBRPA_SWITCH_ON;
        const auto copy_array = [](int count, const auto* source, auto& destination)
        {
            if (count < 0 || (count > 0 && !source))
                throw LIBRPA_INVALID_ARGUMENT("invalid cRPA selection array");
            if (count > 0) destination.assign(source, source + count);
        };
        copy_array(input->n_response_edges, input->response_windows_ha, window.response_windows_ha);
        copy_array(input->n_orbital_edges, input->orbital_windows_ha, window.orbital_windows_ha);
        copy_array(input->n_response_bands, input->response_bands, window.response_bands);
        copy_array(input->n_orbital_bands, input->orbital_bands, window.orbital_bands);
        window.gram_abs = input->gram_abs;
        window.gram_rel = input->gram_rel;
        window.gram_cond_max = input->gram_cond_max;
        window.s_abs = input->s_abs;
        window.s_rel = input->s_rel;
        window.s_cond_max = input->s_cond_max;
        window.residual_tol = input->residual_tol;
        window.t_roundtrip_tol = input->t_roundtrip_tol;
        window.tracking_min_capture = input->tracking_min_capture;
        if (input->n_previous_parent < 0 ||
            (input->n_previous_parent > 0 && !input->previous_parent_phi_ri))
            throw LIBRPA_INVALID_ARGUMENT("invalid cRPA previous-parent frame");
        if (input->n_previous_parent > 0)
        {
            const int rows = input->n_aos * ds->mf.get_n_spinor();
            const int columns = input->n_previous_parent;
            std::size_t offset = 0;
            window.previous_parent_spin_k.resize(ds->mf.get_n_spins());
            for (auto& spin : window.previous_parent_spin_k)
                for (int k = 0; k < input->n_kpoints; ++k)
                {
                    ComplexMatrix phi(rows, columns);
                    for (int i = 0; i < phi.size; ++i, offset += 2)
                        phi.c[i] = {input->previous_parent_phi_ri[offset],
                                    input->previous_parent_phi_ri[offset + 1]};
                    spin.push_back(std::move(phi));
                }
        }
        const std::size_t block_size = static_cast<std::size_t>(input->n_aos) * input->n_aos;
        for (int k = 0; k < input->n_kpoints; ++k)
        {
            overlap.emplace_back(input->n_aos, input->n_aos);
            for (std::size_t i = 0; i < block_size; ++i)
            {
                const auto offset = 2 * (static_cast<std::size_t>(k) * block_size + i);
                overlap.back().c[i] = {input->overlap_k_ri[offset],
                                       input->overlap_k_ri[offset + 1]};
            }
        }
    }
    catch (const std::exception& failure)
    {
        error = failure.what();
    }
    collective_error(ds->comm_h, error);

    CrpaPreparedModel model;
    error.clear();
    try
    {
        std::vector<int> atom_types(ds->basis_wfc.n_atoms);
        for (std::size_t atom = 0; atom < atom_types.size(); ++atom)
            atom_types[atom] = ds->atoms.types.at(atom);
        model = prepare_crpa_model(ds->basis_wfc, ds->basis_convention, atom_types, species,
                                   overlap, ds->pbc, ds->mf, window);
        model.retain_raw_model = input->retain_raw_model == LIBRPA_SWITCH_ON;
    }
    catch (const std::exception& failure)
    {
        error = failure.what();
    }
    collective_error(ds->comm_h, error);

    return model;
}

namespace
{
librpa_int::CrpaContext crpa_context(librpa_int::Dataset& ds, const LibrpaOptions& opts)
{
    const bool shrink = opts.use_shrink_abfs == LIBRPA_SWITCH_ON;
    // cRPA always screens and projects with full v. The independent GW
    // eps/Wc/EXX Coulomb flags do not enter this numerical context.
    return {*ds.p_chi0,   ds.cs_data,
            ds.basis_aux, ds.atpairs_local,
            ds.sinvS,     ds.vq,
            ds.blacs_h,   shrink ? ds.desc_abf_shrink : ds.desc_abf,
            shrink,       opts.sqrt_coulomb_threshold};
}

void initialize_crpa_response(librpa_int::Dataset& ds, const LibrpaOptions& opts)
{
    using namespace librpa_int;
    initialize_ds_global_ddla(ds, opts);
    initialize_ds_tfgrids(ds, opts);
    initialize_ds_atpairs_local(ds, LIBRPA_ROUTING_LIBRI);
    ds.redistribute_coulomb_blacs2ap();
    initialize_ds_chi0(ds, opts);
}

LibrpaCrpaResult* own_crpa_result(librpa_int::CrpaResult values,
                                  librpa_int::CrpaPreparedModel& model,
                                  const librpa_int::Dataset& dataset)
{
    std::unique_ptr<LibrpaCrpaResult> result;
    std::string error;
    try
    {
        result = std::make_unique<LibrpaCrpaResult>();
        result->values = std::move(values);
        result->model = librpa_int::snapshot_crpa_model(model, dataset.mf);
        for (const auto& row : dataset.pbc.latvec_array)
            for (double value : row) result->lattice.push_back(value);
        for (int atom = 0; atom < dataset.basis_wfc.n_atoms; ++atom)
        {
            const auto& position = dataset.atoms.coords.at(atom);
            result->positions.insert(result->positions.end(), {position.x, position.y, position.z});
            result->atom_types.push_back(dataset.atoms.types.at(atom));
            result->ao_counts.push_back(dataset.basis_wfc[atom]);
            const auto& shells = dataset.basis_wfc.get_l_shells(atom);
            result->shell_l.insert(result->shell_l.end(), shells.begin(), shells.end());
            result->shell_offsets.push_back(result->shell_l.size());
        }
        const auto& convention = dataset.basis_convention;
        result->basis_convention = {convention.bloch_phase, convention.bloch_ratom,
                                    convention.order, convention.coeff_m_negative,
                                    convention.coeff_m_positive};
        for (const auto& k : dataset.pbc.kfrac_list)
            result->k_fractional.insert(result->k_fractional.end(), {k.x, k.y, k.z});
        result->k_weights = dataset.pbc.weight_k;
        if (result->k_weights.empty())
            result->k_weights.assign(dataset.mf.get_n_kpoints(), 1.0 / dataset.mf.get_n_kpoints());
        result->response_bands = std::move(model.windows.response_bands);
        result->orbital_bands = std::move(model.windows.orbital_bands);
        result->response_weights = std::move(model.selection);
    }
    catch (const std::exception& failure)
    {
        error = "model snapshot: " + std::string(failure.what());
    }
    collective_error(dataset.comm_h, error);
    return result.release();
}
}  // namespace

LibrpaCrpaResult* librpa_int::api::evaluate_crpa_api_from_response(Dataset& ds,
                                                                   const LibrpaOptions& opts,
                                                                   CrpaPreparedModel& model,
                                                                   const Chi0QMap& full_response)
{
    auto context = crpa_context(ds, opts);
    std::unique_ptr<diele_func> constrained_headwing;
    if (opts.replace_w_head == LIBRPA_SWITCH_ON)
    {
        auto optical_options = opts;
        optical_options.use_symmetry_rpa = optical_options.use_symmetry_gw = LIBRPA_SWITCH_OFF;
        constrained_headwing = make_ds_headwing(ds, optical_options, true, &model.selection);
        constrained_headwing->normalize_gamma_cell_weights = true;
        context.constrained_headwing = constrained_headwing.get();
    }
    std::vector<CorrelatedSubspace*> frames;
    for (const auto& frame : model.frames) frames.push_back(frame.get());
    return own_crpa_result(
        evaluate_crpa_onsite_from_response(context, model.selection, frames, full_response), model,
        ds);
}

LibrpaCrpaResult* librpa_prepare_crpa_model(const LibrpaHandler* h, const LibrpaOptions* options,
                                            const LibrpaCrpaInput* input)
{
    if (!options) throw LIBRPA_INVALID_ARGUMENT("missing cRPA model options");
    auto ds = librpa_int::api::get_dataset_instance(h);
    auto model = librpa_int::api::prepare_crpa_api_model(
        *ds, *options, input, librpa_int::api::CrpaApiMode::ResponseOnly);
    return own_crpa_result({}, model, *ds);
}

LibrpaCrpaResult* librpa_compute_crpa_u_window(const LibrpaHandler* h, const LibrpaOptions* options,
                                               const LibrpaCrpaInput* input)
{
    using namespace librpa_int;
    if (!options) throw LIBRPA_INVALID_ARGUMENT("missing cRPA options");
    auto ds = api::get_dataset_instance(h);
    auto model = api::prepare_crpa_api_model(*ds, *options, input, api::CrpaApiMode::ResponseOnly);
    initialize_crpa_response(*ds, *options);
    auto context = crpa_context(*ds, *options);
    auto full_response = build_crpa_total_response(context);
    return api::evaluate_crpa_api_from_response(*ds, *options, model, full_response);
}

LibrpaCrpaResult* librpa_compute_crpa_window(const LibrpaHandler* h, const LibrpaOptions* options,
                                             const LibrpaCrpaInput* input)
{
    using namespace librpa_int;
    if (!options) throw LIBRPA_INVALID_ARGUMENT("missing cRPA options");
    auto ds = api::get_dataset_instance(h);
    auto model = api::prepare_crpa_api_model(*ds, *options, input);
    const auto& opts = *options;
    initialize_crpa_response(*ds, opts);
    auto context = crpa_context(*ds, opts);
    std::unique_ptr<diele_func> full_headwing, constrained_headwing;
    if (opts.replace_w_head == LIBRPA_SWITCH_ON)
    {
        // Both optical sums retain every original full-BZ state and velocity.
        auto optical_options = opts;
        optical_options.use_symmetry_rpa = optical_options.use_symmetry_gw = LIBRPA_SWITCH_OFF;
        full_headwing = make_ds_headwing(*ds, optical_options, true);
        constrained_headwing = make_ds_headwing(*ds, optical_options, true, &model.selection);
        full_headwing->normalize_gamma_cell_weights = true;
        constrained_headwing->normalize_gamma_cell_weights = true;
        context.full_headwing = full_headwing.get();
        context.constrained_headwing = constrained_headwing.get();
    }
    std::vector<CorrelatedSubspace*> frame_views;
    for (const auto& frame : model.frames) frame_views.push_back(frame.get());
    return own_crpa_result(compute_crpa_onsite(context, model.selection, frame_views), model, *ds);
}

void librpa_delete_crpa_result(LibrpaCrpaResult* result) { delete result; }

void librpa_crpa_model_info(const LibrpaCrpaResult* result, LibrpaCrpaModelInfo* info)
{
    if (!result || !info) throw LIBRPA_INVALID_ARGUMENT("missing cRPA model info");
    const auto& model = result->model;
    const auto& state = model.spin_k.at(0).at(0);
    *info = {model.n_spins,
             model.n_spinor,
             static_cast<int>(model.spin_k.at(0).size()),
             model.n_bands,
             model.n_aos,
             static_cast<int>(result->atom_types.size()),
             state.parent_projector.nr,
             state.output_projector.nr,
             static_cast<int>(model.parent_sites.size()),
             static_cast<int>(model.output_sites.size()),
             model.fermi,
             result->lattice.data(),
             result->positions.data(),
             result->atom_types.data(),
             result->ao_counts.data(),
             result->k_fractional.data(),
             result->k_weights.data(),
             model.output_columns.data(),
             static_cast<int>(result->shell_l.size()),
             result->shell_offsets.data(),
             result->shell_l.data(),
             result->basis_convention.data()};
}

void librpa_crpa_model_k(const LibrpaCrpaResult* result, int spin, int kpoint,
                         LibrpaCrpaModelK* view)
{
    if (!result || !view || spin < 0 || kpoint < 0)
        throw LIBRPA_INVALID_ARGUMENT("missing cRPA model k view");
    const auto& model = result->model;
    const auto& state = model.spin_k.at(spin).at(kpoint);
    const auto packed = [](const librpa_int::ComplexMatrix& matrix)
    { return reinterpret_cast<const double*>(matrix.c); };
    *view = {state.energies.data(),
             state.occupations.data(),
             packed(state.coefficients),
             model.overlap.empty() ? nullptr : packed(model.overlap.at(kpoint)),
             packed(state.parent_phi),
             packed(state.output_phi),
             packed(state.parent_projector),
             packed(state.output_projector),
             packed(state.parent_h),
             packed(state.output_h),
             packed(state.parent_density),
             packed(state.output_density),
             state.parent_capture_min};
}

void librpa_crpa_model_local(const LibrpaCrpaResult* result, int parent, int spin, int site,
                             double* h_ha_ri, double* density_ri)
{
    using namespace librpa_int;
    if (!result || !h_ha_ri || !density_ri || spin < 0 || site < 0)
        throw LIBRPA_INVALID_ARGUMENT("missing local cRPA model buffers");
    const auto& model = result->model;
    const auto& group = (parent ? model.parent_sites : model.output_sites).at(site);
    ComplexMatrix h(group.n_orbitals, group.n_orbitals), density(h);
    const auto& states = model.spin_k.at(spin);
    for (int k = 0; k < states.size(); ++k)
    {
        const auto& hk = parent ? states[k].parent_h : states[k].output_h;
        const auto& nk = parent ? states[k].parent_density : states[k].output_density;
        for (int a = 0; a < group.n_orbitals; ++a)
            for (int b = 0; b < group.n_orbitals; ++b)
            {
                h(a, b) += result->k_weights.at(k) * hk(group.orb_start + a, group.orb_start + b);
                density(a, b) +=
                    result->k_weights.at(k) * nk(group.orb_start + a, group.orb_start + b);
            }
    }
    std::copy_n(reinterpret_cast<const double*>(h.c), 2 * h.size, h_ha_ri);
    std::copy_n(reinterpret_cast<const double*>(density.c), 2 * density.size, density_ri);
}

void librpa_crpa_model_hopping(const LibrpaCrpaResult* result, int parent, int spin,
                               const int* lattice_R, double* h_ha_ri)
{
    using namespace librpa_int;
    if (!result || !lattice_R || !h_ha_ri || spin < 0)
        throw LIBRPA_INVALID_ARGUMENT("missing cRPA hopping buffers");
    const auto& states = result->model.spin_k.at(spin);
    const int n = parent ? states.front().parent_h.nr : states.front().output_h.nr;
    ComplexMatrix h(n, n);
    for (int k = 0; k < states.size(); ++k)
    {
        double dot = 0;
        for (int axis = 0; axis < 3; ++axis)
            dot += result->k_fractional[3 * k + axis] * lattice_R[axis];
        h += (result->k_weights.at(k) * std::exp(std::complex<double>(0, -2 * PI * dot))) *
             (parent ? states[k].parent_h : states[k].output_h);
    }
    std::copy_n(reinterpret_cast<const double*>(h.c), 2 * h.size, h_ha_ri);
}

namespace
{
librpa_int::ComplexMatrix unpack_complex(int rows, int columns, const double* values)
{
    if (rows <= 0 || columns <= 0 || !values)
        throw LIBRPA_INVALID_ARGUMENT("missing complex numerical array");
    librpa_int::ComplexMatrix matrix(rows, columns);
    for (int i = 0; i < matrix.size; ++i) matrix.c[i] = {values[2 * i], values[2 * i + 1]};
    return matrix;
}
}  // namespace

void librpa_crpa_project_one_body(int n_orbitals, int n_bands, const double* projector_ri,
                                  const double* diagonal, double* operator_ri)
{
    if (!diagonal || !operator_ri)
        throw LIBRPA_INVALID_ARGUMENT("missing projected-operator buffers");
    const auto p = unpack_complex(n_orbitals, n_bands, projector_ri);
    const auto result = librpa_int::project_crpa_one_body(p, {diagonal, diagonal + n_bands});
    std::copy_n(reinterpret_cast<const double*>(result.c), 2 * result.size, operator_ri);
}

void librpa_crpa_project_orbitals(int n_aos, int n_spinor, int n_bands, int n_orbitals,
                                  const double* coefficients_ri, const double* overlap_ri,
                                  const double* phi_ri, double* projector_ri)
{
    using namespace librpa_int;
    if (!projector_ri || (n_spinor != 1 && n_spinor != 2))
        throw LIBRPA_INVALID_ARGUMENT("missing orbital-projector buffers or spinor convention");
    const auto c = unpack_complex(n_aos * n_spinor, n_bands, coefficients_ri);
    const auto s = crpa_spinor_overlap(unpack_complex(n_aos, n_aos, overlap_ri), n_spinor);
    const auto phi = unpack_complex(n_aos * n_spinor, n_orbitals, phi_ri);
    const auto gram = transpose(phi, true) * s * phi;
    for (int i = 0; i < n_orbitals; ++i)
        for (int j = 0; j < n_orbitals; ++j)
            if (!std::isfinite(std::abs(gram(i, j))) ||
                std::abs(gram(i, j) - (i == j ? 1.0 : 0.0)) > 1e-10)
                throw LIBRPA_INVALID_ARGUMENT("explicit orbital frame is not metric-orthonormal");
    const auto p = transpose(phi, true) * s * c;
    std::copy_n(reinterpret_cast<const double*>(p.c), 2 * p.size, projector_ri);
}

void librpa_crpa_rotate_projector(int n_bands, int n_old, int n_new, const double* projector_ri,
                                  const double* rotation_ri, double* rotated_projector_ri)
{
    using namespace librpa_int;
    if (!rotated_projector_ri) throw LIBRPA_INVALID_ARGUMENT("missing rotated-projector buffer");
    const auto p = unpack_complex(n_old, n_bands, projector_ri);
    const auto rotation = unpack_complex(n_old, n_new, rotation_ri);
    const auto gram = transpose(rotation, true) * rotation;
    if (n_new > n_old) throw LIBRPA_INVALID_ARGUMENT("projector rotation must be an isometry");
    for (int i = 0; i < n_new; ++i)
        for (int j = 0; j < n_new; ++j)
            if (!std::isfinite(std::abs(gram(i, j))) ||
                std::abs(gram(i, j) - (i == j ? 1.0 : 0.0)) > 1e-10)
                throw LIBRPA_INVALID_ARGUMENT("projector rotation is not an isometry");
    const auto rotated = transpose(rotation, true) * p;
    std::copy_n(reinterpret_cast<const double*>(rotated.c), 2 * rotated.size, rotated_projector_ri);
}

void librpa_crpa_rotate_tensor(int n_old, int n_new, const double* tensor_ri,
                               const double* left_rotation_ri, const double* right_rotation_ri,
                               double* rotated_ri)
{
    if (!rotated_ri) throw LIBRPA_INVALID_ARGUMENT("missing rotated-interaction buffer");
    const auto left = unpack_complex(n_old, n_new, left_rotation_ri);
    const auto right = unpack_complex(n_old, n_new, right_rotation_ri);
    const auto input = unpack_complex(n_old * n_old, n_old * n_old, tensor_ri);
    const auto output = librpa_int::rotate_crpa_tensor(input, left, right);
    std::copy_n(reinterpret_cast<const double*>(output.c), 2 * output.size, rotated_ri);
}

void librpa_crpa_tensor_to_fermion(int n_orbitals, const double* density_ri, double* fermion_ri)
{
    if (!fermion_ri) throw LIBRPA_INVALID_ARGUMENT("missing fermionic-interaction buffer");
    const auto input = unpack_complex(n_orbitals * n_orbitals, n_orbitals * n_orbitals, density_ri);
    for (int i = 0; i < n_orbitals; ++i)
        for (int j = 0; j < n_orbitals; ++j)
            for (int k = 0; k < n_orbitals; ++k)
                for (int l = 0; l < n_orbitals; ++l)
                {
                    const auto value = input(i * n_orbitals + k, j * n_orbitals + l);
                    const int index = ((i * n_orbitals + j) * n_orbitals + k) * n_orbitals + l;
                    fermion_ri[2 * index] = value.real();
                    fermion_ri[2 * index + 1] = value.imag();
                }
}

LibrpaCrpaContinuedResult* librpa_continue_crpa(const LibrpaCrpaResult* result, int n_queries,
                                                const double* queries_ha_ri, int n_parameters,
                                                double tolerance_ha)
{
    if (!result || n_queries <= 0 || !queries_ha_ri)
        throw LIBRPA_INVALID_ARGUMENT("missing cRPA continuation queries");
    std::vector<std::complex<double>> queries;
    for (int i = 0; i < n_queries; ++i)
        queries.emplace_back(queries_ha_ri[2 * i], queries_ha_ri[2 * i + 1]);
    auto output = std::make_unique<LibrpaCrpaContinuedResult>();
    output->tensors =
        librpa_int::continue_crpa_interactions(result->values, queries, n_parameters, tolerance_ha);
    return output.release();
}

void librpa_delete_crpa_continued_result(LibrpaCrpaContinuedResult* result) { delete result; }
int librpa_crpa_continued_size(const LibrpaCrpaContinuedResult* result)
{
    if (!result) throw LIBRPA_INVALID_ARGUMENT("missing continued cRPA result");
    return result->tensors.size();
}
void librpa_crpa_continued_tensor(const LibrpaCrpaContinuedResult* result, int tensor, int query,
                                  LibrpaCrpaContinuedTensor* view)
{
    if (!result || !view || tensor < 0 || query < 0)
        throw LIBRPA_INVALID_ARGUMENT("missing continued cRPA tensor view");
    const auto& record = result->tensors.at(tensor);
    *view = {record.atom_index,
             record.spin_left,
             record.spin_right,
             record.n_orbitals,
             reinterpret_cast<const double*>(record.values.at(query).c),
             record.status.at(query).data(),
             record.spread.at(query).data(),
             record.max_source_residual,
             record.max_holdout_residual};
}

void librpa_crpa_model_site(const LibrpaCrpaResult* result, int parent, int site, int* atom,
                            int* column_start, int* n_orbitals)
{
    if (!result || !atom || !column_start || !n_orbitals || site < 0)
        throw LIBRPA_INVALID_ARGUMENT("missing cRPA model site");
    const auto& group = (parent ? result->model.parent_sites : result->model.output_sites).at(site);
    *atom = group.atom_index;
    *column_start = group.orb_start;
    *n_orbitals = group.n_orbitals;
}

const char* librpa_crpa_model_label(const LibrpaCrpaResult* result, int parent, int site,
                                    int orbital)
{
    if (!result || site < 0 || orbital < 0)
        throw LIBRPA_INVALID_ARGUMENT("missing cRPA model label");
    const auto& group = (parent ? result->model.parent_sites : result->model.output_sites).at(site);
    return group.orbital_labels.at(orbital).c_str();
}

int librpa_crpa_result_size(const LibrpaCrpaResult* result)
{
    if (!result) throw LIBRPA_INVALID_ARGUMENT("null cRPA result");
    const auto count = result->values.tensors.size() * result->values.frequencies.size();
    if (count > static_cast<std::size_t>(std::numeric_limits<int>::max()))
        throw LIBRPA_RUNTIME_ERROR("cRPA tensor count exceeds the public integer range");
    return static_cast<int>(count);
}

void librpa_crpa_result_tensor(const LibrpaCrpaResult* result, int index, LibrpaCrpaTensor* tensor)
{
    if (!tensor || index < 0 || index >= librpa_crpa_result_size(result))
        throw LIBRPA_INVALID_ARGUMENT("invalid cRPA tensor query");
    const auto nf = result->values.frequencies.size();
    const auto& record = result->values.tensors[index / nf];
    const auto& site = result->values.sites.at(record.site_index);
    tensor->atom_index = site.atom_index;
    tensor->spin_left = record.spin_left;
    tensor->spin_right = record.spin_right;
    tensor->n_orbitals = site.n_orbitals;
    tensor->frequency_ha = result->values.frequencies[index % nf];
    tensor->bare_ri = reinterpret_cast<const double*>(record.bare.c);
    tensor->partially_screened_ri = reinterpret_cast<const double*>(record.u.at(index % nf).c);
    tensor->fully_screened_ri =
        record.w.empty() ? nullptr : reinterpret_cast<const double*>(record.w.at(index % nf).c);
}

void librpa_crpa_result_bands(const LibrpaCrpaResult* result, int response_window, int spin,
                              int kpoint, int* count, const int** bands)
{
    if (!result || !count || !bands || spin < 0 || kpoint < 0)
        throw LIBRPA_INVALID_ARGUMENT("invalid cRPA window query");
    const auto& row =
        (response_window ? result->response_bands : result->orbital_bands).at(spin).at(kpoint);
    *count = static_cast<int>(row.size());
    *bands = row.data();
}

void librpa_crpa_result_weights(const LibrpaCrpaResult* result, int spin, int kpoint, int* count,
                                const double** weights)
{
    if (!result || !count || !weights || spin < 0 || kpoint < 0)
        throw LIBRPA_INVALID_ARGUMENT("invalid cRPA endpoint-weight query");
    const auto& row = result->response_weights.at(spin).at(kpoint);
    *count = static_cast<int>(row.size());
    *weights = row.data();
}

const char* librpa_crpa_result_orbital_label(const LibrpaCrpaResult* result, int atom, int orbital)
{
    if (!result || orbital < 0) throw LIBRPA_INVALID_ARGUMENT("invalid cRPA orbital query");
    for (const auto& site : result->values.sites)
        if (site.atom_index == atom) return site.orbital_labels.at(orbital).c_str();
    throw LIBRPA_INVALID_ARGUMENT("cRPA output does not contain the requested atom");
}
