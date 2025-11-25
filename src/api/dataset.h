#pragma once

#include <memory>

#include "../mpi/base_blacs.h"
#include "../mpi/base_mpi.h"
#include "../core/atomic_basis.h"
#include "../core/chi0.h"
#include "../core/exx.h"
#include "../core/geometry.h"
#include "../core/gw.h"
#include "../core/meanfield.h"
#include "../core/pbc.h"
#include "../core/ri.h"
#include "../core/timefreq.h"

namespace librpa_int
{

#define is_null_dataset_ptr(p) p == nullptr

/*!
 * Core object to hold input and output data
 */
class Dataset
{
public:
    /* Variables */
    // Environment control
    MpiCommHandler comm_h;
    BlacsCtxtHandler blacs_ctxt_h;

    // Physical system.
    //! Basic set functions for wave function expansion.
    AtomicBasis basis_wfc;
    //! Auxiliary basic set functions for RI
    AtomicBasis basis_abf;
    Atoms atoms;
    PeriodicBoundaryData pbc;

    // Input data.
    MeanField mf;
    TFGrids tfg;
    Cs_LRI cs_data;
    atpair_k_cplx_mat_t vq;
    atpair_k_cplx_mat_t vq_cut;

    // Computation data.
    // All computation data objects should be contained as pointers
    // and are created only when requested.
    std::unique_ptr<Exx> p_exx = nullptr;
    std::unique_ptr<Chi0> p_chi0 = nullptr;
    std::unique_ptr<G0W0> p_g0w0 = nullptr;

    /* Constructors and destructors */
    Dataset(MPI_Comm comm);
    ~Dataset() { free(); }
    void free();

    /* Disable copy */
    Dataset(const Dataset &) = delete;
    Dataset &operator=(const Dataset &) = delete;
};

}
