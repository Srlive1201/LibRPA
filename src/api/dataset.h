#pragma once

#include <memory>
#include <unordered_map>

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
#include "../math/matrix_m.h"
#include "../mpi/base_blacs.h"
#include "../mpi/base_mpi.h"

namespace librpa_int
{

/*!
 * Core object to hold runtime environment, input and output data
 */
class Dataset
{
private:
    bool comm_blacs_coul_initialized_;
    bool coul_blacs2ap_redistributed_;
public:
    /* Member variables */
    // Environment control
    // Global MPI communicators and BLACS context handlers
    MpiCommHandler comm_h;
    BlacsCtxtHandler blacs_h;
    // Communicators for Coulomb matrices 2D input and re-distribution
    MpiCommHandler comm_coul_h;
    MpiCommHandler comm_coul_inter_q_h;
    MpiCommHandler comm_coul_intra_q_h;
    BlacsCtxtHandler blacs_coul_intra_q_h;
    ArrayDesc desc_coul_intra_q;
    // Communicators for (unordered) atom pairs.
    // All atom pairs are distributed among processes in the communicator.
    MpiCommHandler comm_ap_h;
    // Communicators among unit cell vectors.
    // All BvK cell vectors are distributed among processes in the communicator.
    MpiCommHandler comm_R_h;
    //! Array descriptor for matrices of wave-function basis (using blacs_h)
    ArrayDesc desc_wfc;
    //! Array descriptor for matrices of auxiliary basis set size (using blacs_h)
    ArrayDesc desc_abf;
    //! Atom pairs on current process for atomic-basis matrix data
    std::vector<atpair_t> atpairs_local;
    //! Distribution of unique atom-pairs on all processes for atomic-basis matrix data
    std::unordered_map<int, std::set<atpair_t>> atpairs_unique_all;
    //! Unit cell vectors on current process
    std::vector<Vector3_Order<int>> Rs_local;

    // Physical system.
    //! Handling object for basic set functions for wave function expansion.
    AtomicBasis basis_wfc;
    //! Handling object for auxiliary basic set functions for RI
    AtomicBasis basis_aux;
    //! Atomic structure
    Atoms atoms;
    //! Periodic boundary setting
    PeriodicBoundaryData pbc;

    // Input data.
    //! Mean-field starting point
    MeanField mf;
    //! Time-frequency grids
    TFGrids tfg;
    //! Real-space RI coefficient tensors (local RI)
    Cs_LRI cs_data;
    // atom-pair distribution of Coulomb matrices
    atpair_k_cplx_mat_t vq;
    atpair_k_cplx_mat_t vq_cut;
    // Local index boundaries and 2D distribution of Coulomb matrices
    // Indices are just for data parsing and blacs context initialization.
    // Further handling should rely on comm_coul_*/blacs_coul_* communicators.
    // Lower bounds are included, and upper bounds are excluded.
    // Bare and cut Coulombs are enforced to use the same layout.
    int vq_lbrow = -1, vq_ubrow = -1;
    int vq_lbcol = -1, vq_ubcol = -1;
    std::map<Vector3_Order<double>, Matz> vq_block_loc;
    std::map<Vector3_Order<double>, Matz> vq_cut_block_loc;
    // Macroscopic dielectric functions at imaginary frequencies
    // Only used for head correction of GW calculation
    std::vector<double> epsmacs_imagfreq;
    std::vector<double> omegas_imagfreq;

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

    void initialize_comm_blacs_coul();
    void redistribute_coulomb_blacs2ap();
    void finalize_comm_blacs_coul();

    // Splitting global communicators to atom pairs and unit cell vectors
    void initialize_comm_ap_r();
    void finalize_comm_ap_r();

    /* Disable copy */
    Dataset(const Dataset &) = delete;
    Dataset &operator=(const Dataset &) = delete;
};

typedef std::shared_ptr<Dataset> dataset_ptr_t;

}
