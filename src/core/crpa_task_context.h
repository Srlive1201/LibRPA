#pragma once

#include <functional>
#include <memory>

#include "atom.h"
#include "atomic_basis.h"
#include "chi0.h"
#include "dielecmodel.h"
#include "geometry.h"
#include "librpa_options.h"
#include "meanfield.h"
#include "pbc.h"
#include "ri.h"

namespace librpa_int
{
//! Borrow the loaded numerical data without changing the API Dataset's ownership.
struct CrpaTaskContext
{
    MeanField& mf;
    PeriodicBoundaryData& pbc;
    Atoms& atoms;
    AtomicBasis &basis_wfc, &basis_aux, &basis_aux_shrink;
    BasisConvention& basis_convention;
    MpiCommHandler& comm_h;
    BlacsCtxtHandler& blacs_h;
    ArrayDesc &desc_abf, &desc_abf_shrink;
    std::vector<atpair_t>& atpairs_local;
    Cs_LRI& cs_data;
    std::map<Vector3_Order<double>, ComplexMatrix>& sinvS;
    atpair_k_cplx_mat_t& vq;
    std::unique_ptr<Chi0>& p_chi0;
    std::function<void(const LibrpaOptions&)> initialize_response;
    std::function<std::unique_ptr<diele_func>(const LibrpaOptions&, bool,
                                              const HeadwingBandSelection*)>
        make_headwing;
};
}  // namespace librpa_int
