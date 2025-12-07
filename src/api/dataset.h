#pragma once

#include <memory>

#include "../core/atom.h"
#include "../core/atomic_basis.h"
#include "../core/chi0.h"
#include "../core/exx.h"
#include "../core/geometry.h"
#include "../core/gw.h"
#include "../core/meanfield.h"
#include "../core/pbc.h"
#include "../core/ri.h"
#include "../core/timefreq.h"
#include "../mpi/base_blacs.h"
#include "../mpi/base_mpi.h"

namespace librpa_int
{

/*!
 * Core object to hold runtime environemtn, input and output data
 */
class Dataset
{
public:
    /* Member variables */
    // Environment control
    MpiCommHandler comm_h;
    BlacsCtxtHandler blacs_ctxt_h;
    //! Array descriptor for matrices of auxiliary basis set size
    ArrayDesc desc_abf;
    //! Distribution of atom-pairs on current process for atomic-basis matrix data
    std::vector<atpair_t> atpairs_local;

    // Physical system.
    //! Handliing boject for basic set functions for wave function expansion.
    AtomicBasis basis_wfc;
    //! Handliing boject for auxiliary basic set functions for RI
    AtomicBasis basis_aux;
    Atoms atoms;
    PeriodicBoundaryData pbc;

    // Input data.
    MeanField mf;
    TFGrids tfg;
    Cs_LRI cs_data;
    // atom-pair distribution of Coulomb matrices
    atpair_k_cplx_mat_t vq;
    atpair_k_cplx_mat_t vq_cut;
    // 2D distribution of Coulomb matrices
    std::map<Vector3_Order<double>, ComplexMatrix> vq_block_loc;
    std::map<Vector3_Order<double>, ComplexMatrix> vq_cut_block_loc;

    // Output data, held by computation objects
    // All computation data objects should be contained here as pointers
    // and are created only when requested.
    std::unique_ptr<Exx> p_exx;
    std::unique_ptr<Chi0> p_chi0;
    std::unique_ptr<G0W0> p_g0w0;

    /* Constructors and destructors */
    Dataset(MPI_Comm comm);
    ~Dataset() { free(); }
    void free();

    /* Disable copy */
    Dataset(const Dataset &) = delete;
    Dataset &operator=(const Dataset &) = delete;
};

typedef std::shared_ptr<Dataset> dataset_ptr_t;

}
