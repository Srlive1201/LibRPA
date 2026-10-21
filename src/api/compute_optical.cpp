#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <sstream>
#include <vector>

#include "../core/optical_task.h"
#include "../mpi/collective_stage.h"
#include "../utils/constants.h"
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

void librpa_prepare_native_velocity(LibrpaHandler *h, const LibrpaOptions *options,
                                    const LibrpaNativeVelocityInput *input, int prepare_optical)
{
    const auto dataset = librpa_int::api::get_dataset_instance(h);
    librpa_int::prepare_native_velocity(*dataset, options, input, prepare_optical);
}

void librpa_prepare_supplied_velocity(LibrpaHandler *h, const LibrpaOptions *options,
                                      const LibrpaSuppliedVelocityInput *input, const int need_wing)
{
    const auto ds = librpa_int::api::get_dataset_instance(h);
    librpa_int::prepare_supplied_velocity(*ds, options, input, need_wing);
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
    const auto &optical_mf = ds->p_headwing->get_meanfield_df();
    view.n_spins = optical_mf.get_n_spins();
    view.n_kpoints = optical_mf.get_n_kpoints();
    view.n_states = optical_mf.get_n_states();
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
