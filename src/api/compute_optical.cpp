#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <vector>

#include "../core/epsilon.h"
#include "../core/optical_task.h"
#include "../mpi/collective_stage.h"
#include "../utils/constants.h"
#include "dataset_helper.h"
#include "instance_manager.h"
#include "librpa_input.h"
#include "librpa_optical.h"

namespace
{
using namespace librpa_int;

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
    const auto ds = api::get_dataset_instance(h);
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
    auto result = build_native_velocity(ds->mf, ds->pbc, comm, root_input);
    collective_all_stage(comm, "G0W0 native velocity setter",
                         [&]
                         {
                             librpa_set_velocity_matrix_packed(
                                 h, ds->mf.get_n_spins(), ds->mf.get_n_kpoints(),
                                 ds->mf.get_n_states(),
                                 reinterpret_cast<const double *>(result.velocity.data()));
                         });
    if (comm.is_root())
        global::lib_printf(
            "G0W0 native velocity: metric %.8e frame_Ha %.8e hermiticity_Ha %.8e "
            "diagonal_imag_HaBohr %.8e; SCF E/C/occupations retained\n",
            result.residuals[0], result.residuals[1], result.residuals[2], result.residuals[3]);
    auto &coulomb_normalization = result.coulomb_normalization;

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
