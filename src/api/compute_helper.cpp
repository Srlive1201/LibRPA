#include "compute_helper.h"

#include <librpa_enums.h>
#include "librpa_optical.h"

#include <algorithm>
#include <memory>
#include <set>
#include <string>

#include "../core/dielecmodel.h"
#include "../io/global_io.h"
#include "instance_manager.h"
#include "../io/stl_io_helper.h"
#include "../utils/error.h"

namespace librpa_int::api
{

// Build a deduplicated list of requested k-points across ranks. This lets
// extraction APIs exchange only the k-resolved quantities callers asked for.
std::vector<int> collect_requested_iks(const MpiCommHandler &comm_h, const int n_kpts_this,
                                       const int *iks_this, const int n_kpoints)
{
    std::vector<int> nks_all(comm_h.nprocs, 0);
    nks_all[comm_h.myid] = n_kpts_this;
    comm_h.allreduce(MPI_IN_PLACE, nks_all.data(), comm_h.nprocs, MPI_SUM);

    const int nk_max = *std::max_element(nks_all.cbegin(), nks_all.cend());
    if (nk_max == 0) return {};

    std::vector<int> iks_all(comm_h.nprocs * nk_max, -1);
    int bad_ik_local = 0;
    int bad_ik_value_local = -1;
    for (int i = 0; i != n_kpts_this; ++i)
    {
        const int ik = iks_this[i];
        if (ik < 0 || ik >= n_kpoints)
        {
            bad_ik_local = 1;
            bad_ik_value_local = ik;
            continue;
        }
        iks_all[comm_h.myid * nk_max + i] = ik;
    }

    int bad_ik = 0;
    int bad_ik_value = -1;
    comm_h.allreduce(&bad_ik_local, &bad_ik, 1, MPI_MAX);
    comm_h.allreduce(&bad_ik_value_local, &bad_ik_value, 1, MPI_MAX);
    if (bad_ik)
        throw LIBRPA_RUNTIME_ERROR("requested k-point index out of range: " +
                                  std::to_string(bad_ik_value));

    comm_h.allreduce(MPI_IN_PLACE, iks_all.data(), static_cast<int>(iks_all.size()), MPI_MAX);

    std::set<int> unique_iks;
    for (int pid = 0; pid != comm_h.nprocs; ++pid)
    {
        for (int i = 0; i != nks_all[pid]; ++i)
        {
            unique_iks.insert(iks_all[pid * nk_max + i]);
        }
    }
    return {unique_iks.cbegin(), unique_iks.cend()};
}

AtomPairBvKRemap<atom_t> build_band_bvk_remap(const Atoms &atoms,
                                              const PeriodicBoundaryData &pbc,
                                              const int remap_convention)
{
    using global::ofs_myid;

    if (remap_convention != 0 && remap_convention != 1)
        throw LIBRPA_RUNTIME_ERROR("Invalid BvK remap convention: " + std::to_string(remap_convention));

    const bool is_ready = pbc.is_latt_set() && atoms.is_set() && atoms.is_frac_set() && !pbc.Rlist.empty();
    AtomPairBvKRemap<atom_t> remap;
    if (is_ready)
    {
        remap.build(atoms.coords_frac, pbc.Rlist, pbc.period, pbc.latvec, remap_convention);
    }

    ofs_myid << "Final BvK remapping for band interpolation:\n";
    ofs_myid << "| remap convention  : " << remap_convention << "\n";
    ofs_myid << "| input ready       : " << (is_ready ? "true" : "false") << "\n";
    ofs_myid << "| lattice set       : " << (pbc.is_latt_set() ? "true" : "false") << "\n";
    ofs_myid << "| atoms set         : " << (atoms.is_set() ? "true" : "false") << "\n";
    ofs_myid << "| atom frac set     : " << (atoms.is_frac_set() ? "true" : "false") << "\n";
    ofs_myid << "| period            : " << pbc.period << "\n";
    ofs_myid << "| Rlist size        : " << pbc.Rlist.size() << "\n";
    ofs_myid << "| atom-pair entries : " << remap.size() << "\n";
    for (const auto &[atom_pair, R_map]: remap.data())
    {
        for (const auto &[R, R_bvks]: R_map)
        {
            // Only print 1-to-1 remapping with debug
            if (R_bvks.size() < 2 && !global::should_output(LIBRPA_VERBOSE_DEBUG)) continue;
            ofs_myid << "| " << atom_pair << " " << R << " -> " << R_bvks << "\n";
        }
    }
    ofs_myid << std::flush;
    return remap;
}

}

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
