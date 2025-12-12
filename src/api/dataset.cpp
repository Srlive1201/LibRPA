#include "dataset.h"
#include <ios>

#include "../utils/profiler.h"
#include "../io/stl_io_helper.h"

namespace librpa_int
{

Dataset::Dataset(MPI_Comm comm)
    : comm_blacs_coul_initialized_(false),
      comm_h(comm, true),
      blacs_h(comm),
      comm_coul_h(),
      comm_coul_inter_q_h(),
      comm_coul_intra_q_h(),
      blacs_coul_intra_q_h(),
      desc_coul(),
      desc_wfc(),
      desc_abf(),
      atpairs_local(),
      basis_wfc(),
      basis_aux(),
      atoms(),
      pbc(),
      mf(),
      tfg(),
      cs_data(),
      vq(), vq_cut(),
      vq_lbrow(-1),  vq_ubrow(-1), vq_lbcol(-1), vq_ubcol(-1),
      vq_block_loc(), vq_cut_block_loc(),
      epsmacs_imagfreq(),
      omegas_imagfreq(),
      p_exx(nullptr),
      p_chi0(nullptr),
      p_g0w0(nullptr)
{
    blacs_h.init();
    // TODO more flexible process grid initialization
    blacs_h.set_square_grid();
}

void Dataset::free()
{
    finalize_comm_blacs_coul();
}

void Dataset::initialize_comm_blacs_coul()
{
    using std::cout;
    using std::endl;
    using global::ofs_myid;
    // Reset first
    if (comm_blacs_coul_initialized_) finalize_comm_blacs_coul();
    // Check if local Coulomb matrices are available on one of the processes
    const int has_coul = (vq_block_loc.size() > 0 || vq_cut_block_loc.size() > 0);
    int has_coul_somewhere;
    comm_h.allreduce(&has_coul, &has_coul_somewhere, 1, MPI_MAX);
    if (has_coul_somewhere == 0)
    {
        // No process has Coulomb local blocks, data not parsed
        // NOTE: this immediate return may be changed in the future,
        //       depending on how we want to use the block-cyclic local matrices.
        return;
    }

    global::profiler.start(__FUNCTION__);
    if (!this->basis_aux.initialized())
    {
        throw LIBRPA_RUNTIME_ERROR("Auxiliary basis must be set first");
    }

    // Below we try to restore the two-level parallel distribution of matrices from input
    MPI_Comm comm_coul;
    MPI_Comm_split(comm_h.comm, has_coul, comm_h.myid, &comm_coul);
    if (has_coul)
        comm_coul_h.reset_comm(comm_coul, true);
    else
        MPI_Comm_free(&comm_coul);

    ofs_myid << "Initializing Coulomb BLACS communicators, context and array descriptor" << endl;
    if (comm_coul_h.is_initialized())
    {
        MPI_Comm comm_intra_q, comm_inter_q;
        const auto &block_loc = vq_block_loc.size() > 0? vq_block_loc : vq_cut_block_loc;
        std::vector<int> iqs;
        for (const auto &[q, _]: block_loc)
            iqs.emplace_back(this->pbc.get_k_index_ibz(q));
        // Assuming q-points are correctly grouped among the processes,
        // we use the first q-point index to identify different groups
        std::sort(iqs.begin(), iqs.end());
        MPI_Comm_split(comm_coul_h.comm, iqs[0], comm_coul_h.myid, &comm_intra_q);
        comm_coul_intra_q_h.reset_comm(comm_intra_q, true);
        ofs_myid << "Coulomb intra-q comm/myid/nprocs: "
                 << comm_coul_intra_q_h.comm << " "
                 << comm_coul_intra_q_h.myid << " "
                 << comm_coul_intra_q_h.nprocs << endl;
        MPI_Comm_split(comm_coul_h.comm, comm_coul_intra_q_h.myid, comm_coul_h.myid, &comm_inter_q);
        comm_coul_inter_q_h.reset_comm(comm_inter_q, true);
        ofs_myid << "Coulomb inter-q comm/myid/nprocs: "
                 << comm_coul_inter_q_h.comm << " "
                 << comm_coul_inter_q_h.myid << " "
                 << comm_coul_inter_q_h.nprocs << endl;
    }
    // Till now, the communicators are initialized.
    // Now try to restore the BLACS context and the array descriptor
    if (comm_coul_intra_q_h.is_initialized())
    {
        const int myid = comm_coul_intra_q_h.myid;
        const int nprocs = comm_coul_intra_q_h.nprocs;
        std::vector<int> lbrow_all(nprocs, 0), lbcol_all(nprocs, 0);
        lbrow_all[myid] = vq_lbrow;
        lbcol_all[myid] = vq_lbcol;
        comm_coul_intra_q_h.allreduce(MPI_IN_PLACE, lbrow_all.data(), nprocs, MPI_SUM);
        comm_coul_intra_q_h.allreduce(MPI_IN_PLACE, lbcol_all.data(), nprocs, MPI_SUM);
        ofs_myid << "Collected lbrow: " << lbrow_all << endl;
        ofs_myid << "Collected lbcol: " << lbcol_all << endl;
        // Get number of rows and columns, as well as the process which stores the head of the matrix
        int nprows = 0, npcols = 0, myid_src = 0;
        for (int pid = 0; pid < nprocs; pid++)
        {
            if (lbrow_all[pid] == 0)
            {
                nprows++;
                if (lbcol_all[pid] == 0) myid_src = pid;
            }
            if (lbcol_all[pid] == 0) npcols++;
        }
        // Checking the boundary of the first two processes for the 2D layout
        CTXT_LAYOUT layout = CTXT_LAYOUT::R;
        if (nprocs > 1)
        {
            if (lbrow_all[0] != lbrow_all[1]) layout = CTXT_LAYOUT::C;
        }
        ofs_myid << "Coulomb BLACS context nprows/npcols/layout: " << nprows << " " << npcols;
        if (layout == CTXT_LAYOUT::C)
            ofs_myid << " column-major" << endl;
        else
            ofs_myid << " row-major" << endl;
        blacs_coul_intra_q_h.reset_comm(comm_coul_intra_q_h.comm, true);
        blacs_coul_intra_q_h.set_grid(nprows, npcols, layout);
        // Now we have restored the BLACS communicator

        // Set up the array descriptor.
        // Decide the row and column blocks size
        int nbrows = vq_ubrow - vq_lbrow;
        int nbcols = vq_ubcol - vq_lbcol;
        int mb, nb;
        comm_coul_intra_q_h.allreduce(&nbrows, &mb, 1, MPI_MAX);
        comm_coul_intra_q_h.allreduce(&nbcols, &nb, 1, MPI_MAX);
        ofs_myid << "Coulomb array descriptor block size (mb/nb): " << mb << " " << nb << endl;
        int irsrc = 0, icsrc = 0;
        blacs_coul_intra_q_h.get_pcoord(myid_src, irsrc, icsrc);
        ofs_myid << "Coulomb array descriptor source pid (ir/ic): " << irsrc << " " << icsrc << endl;
        desc_coul.reset_handler(blacs_coul_intra_q_h);
        const int n_aux = this->basis_aux.nb_total;
        desc_coul.init(n_aux, n_aux, mb, nb, irsrc, icsrc);
    }

    ofs_myid << "comm_coul_h.is_initialized() = " << std::boolalpha << comm_coul_h.is_initialized() << std::endl;

    if (comm_coul_h.is_initialized())
    {
        if (comm_coul_h.is_root())
        {
            cout << "BLACS environment of 2D Coulomb matrices is restored" << endl;
        }
        comm_coul_h.barrier();
        for (int pid = 0; pid < comm_coul_h.nprocs; pid++)
        {
            // if (pid == comm_coul_h.myid)
            // {
            //     std::cout << blacs_coul_intra_q_h.info() << std::endl;
            // }
            comm_coul_h.barrier();
        }
    }

    comm_blacs_coul_initialized_ = true;
    global::profiler.stop(__FUNCTION__);
}

void Dataset::finalize_comm_blacs_coul()
{
    if (!comm_blacs_coul_initialized_) return;
    global::profiler.start(__FUNCTION__);
    desc_coul.reset_handler();
    blacs_coul_intra_q_h.reset_comm();
    comm_coul_inter_q_h.free_comm();
    comm_coul_intra_q_h.free_comm();
    comm_coul_h.free_comm();
    comm_blacs_coul_initialized_ = false;
    global::profiler.stop(__FUNCTION__);
}

}
