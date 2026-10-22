#pragma once

#include <functional>
#include <memory>
#include <vector>

#include "../atom.h"
#include "../atomic_basis.h"
#include "../chi0.h"
#include "../dielecmodel.h"
#include "../exx.h"
#include "../geometry.h"
#include "../gw.h"
#include "../meanfield.h"
#include "../pbc.h"
#include "../symmetry_context.h"
#include "../timefreq.h"
#include "librpa_options.h"

namespace librpa_int::qsgw
{
using ResponseObserver = std::function<void(const Chi0QMap&)>;

//! Borrowed numerical state for the QSGW task. Its API owner retains the dataset.
//! Existing G0W0 API orchestration is supplied at this boundary and stays in place.
struct TaskContext
{
    Atoms& atoms;
    AtomicBasis& basis_wfc;
    BlacsCtxtHandler& blacs_h;
    MpiCommHandler& comm_h;
    bool& is_band_calc_done;
    std::vector<Vector3_Order<double>>& kfrac_band_list;
    MeanField& mf;
    MeanField& mf_band;
    std::unique_ptr<Exx>& p_exx;
    std::unique_ptr<G0W0>& p_g0w0;
    std::unique_ptr<diele_func>& p_headwing;
    PeriodicBoundaryData& pbc;
    SymmetryContext& symmetry_context;
    TFGrids& tfg;
    velocity_matrix_t& velocity_matrix;

    std::function<void()> invalidate_compute_objects;
    std::function<void(int, int, int, int, int)> initialize_scf_dimension;
    std::function<void(const LibrpaOptions&)> initialize_tfgrids;
    std::function<void(const LibrpaOptions&, bool)> initialize_headwing;
    std::function<void(bool)> initialize_symmetry;
    std::function<void(const LibrpaOptions&, const ResponseObserver&)> build_sigma;
};
}  // namespace librpa_int::qsgw
