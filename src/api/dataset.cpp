#include "dataset.h"
#include <ios>
#include <memory>
#include <vector>

#include "../core/utils_atomic_basis_blacs.h"
#include "../math/utils_matrix_m_mpi.h"
#include "../utils/profiler.h"
#include "../io/stl_io_helper.h"

namespace librpa_int
{

Dataset::Dataset(MPI_Comm comm)
    : comm_blacs_coul_initialized_(false),
      coul_blacs2ap_redistributed_(false),
      comm_h(comm, true),
      blacs_h(comm),
      comm_coul_h(),
      comm_coul_inter_q_h(),
      comm_coul_intra_q_h(),
      blacs_coul_intra_q_h(),
      desc_coul_intra_q(),
      comm_ap_h(),
      comm_R_h(),
      desc_wfc(),
      desc_abf(),
      atpairs_local(), atpairs_unique_all(),
      Rs_local(),
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
    // TODO: two-level-parallelization initialization
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
    // Already initialized, do not do it again since 2D blocks may have been released.
    // To run it again, manually finalize it first
    if (comm_blacs_coul_initialized_) return;
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
                npcols++;
                if (lbcol_all[pid] == 0) myid_src = pid;
            }
            if (lbcol_all[pid] == 0) nprows++;
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
        desc_coul_intra_q.reset_handler(blacs_coul_intra_q_h);
        const int n_aux = this->basis_aux.nb_total;
        desc_coul_intra_q.init(n_aux, n_aux, mb, nb, irsrc, icsrc);
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

void Dataset::redistribute_coulomb_blacs2ap()
{
    using std::endl;
    using global::ofs_myid;
    using global::profiler;

    // Already redistributed, or blacs data has been parsed
    if (coul_blacs2ap_redistributed_) return;
    initialize_comm_blacs_coul();
    profiler.start(__FUNCTION__);
    // Step 1: redistribute to global blacs_h if
    // - some process does not have block data
    // - there are more than one comm_coul_inter_q_h, meaning Nabs x Nabs matrix is distributed only on a subset of all processes
    // - the shape and layout of blacs_coul_intra_q_h is different from blacs_h
    // - the descriptor is different from the global one
    // Otherwise, blacs_coul_intra_q_h is equivalent to blacs_h, no need to redistribute
    const int same_shape_blacs = blacs_h.nprows == blacs_coul_intra_q_h.nprows &&
                                 blacs_h.npcols == blacs_coul_intra_q_h.npcols &&
                                 blacs_h.layout == blacs_coul_intra_q_h.layout;
    int same_shape_blacs_all;
    comm_h.allreduce(&same_shape_blacs, &same_shape_blacs_all, 1, MPI_MIN);
    const int same_desc = desc_abf.m() == desc_coul_intra_q.m() &&
                          desc_abf.n() == desc_coul_intra_q.n() &&
                          desc_abf.mb() == desc_coul_intra_q.mb() &&
                          desc_abf.nb() == desc_coul_intra_q.nb() &&
                          desc_abf.irsrc() == desc_coul_intra_q.irsrc() &&
                          desc_abf.icsrc() == desc_coul_intra_q.icsrc();
    int same_desc_all;
    comm_h.allreduce(&same_desc, &same_desc_all, 1, MPI_MIN);

    const int n_comms_q = comm_coul_inter_q_h.nprocs;
    int n_comms_q_min, n_comms_q_max;
    comm_h.allreduce(&n_comms_q, &n_comms_q_min, 1, MPI_MIN);
    comm_h.allreduce(&n_comms_q, &n_comms_q_max, 1, MPI_MAX);
    // No block matrix data parsed, no need to redistribute
    if (n_comms_q_max == 0)
    {
        profiler.stop(__FUNCTION__);
        return;
    }

    const int n_aux = basis_aux.nb_total;
    // Major for communication with p?gemr2d
    const MAJOR major_comm = MAJOR::COL;
    std::map<Vector3_Order<double>, Matz> vq_redist;
    std::map<Vector3_Order<double>, Matz> vq_cut_redist;
    if (n_comms_q_min == 0 || n_comms_q_max > 1 || !same_shape_blacs_all || !same_desc_all)
    {
        auto mat_local = init_local_mat<cplxdb>(desc_abf, major_comm);
        for (int i_comm_q = 0; i_comm_q < n_comms_q_max; i_comm_q++)
        {
            std::vector<int> desc(9, n_aux);
            desc[6] = 0;
            desc[7] = 0;
            if (comm_coul_inter_q_h.nprocs > 0 && comm_coul_inter_q_h.myid == i_comm_q)
            {
                // source, copy the descriptor
                desc = std::vector<int>(desc_coul_intra_q.desc, desc_coul_intra_q.desc + 9);
            }
            else
            {
                // not source, use a dummy descriptor: dtype = 1, ctxt = -1;
                desc[0] = 1;
                desc[1] = -1;
            }
            std::vector<std::pair<std::map<Vector3_Order<double>, Matz>&, std::map<Vector3_Order<double>, Matz>&>>
                vq_src_dist_pair{ { vq_block_loc, vq_redist}, {vq_cut_block_loc, vq_cut_redist} };
            for (auto [vq_src, vq_dst]: vq_src_dist_pair)
            {
                // Get the number and coordinates of q-points for this communicator
                int nq_this = vq_src.size();
                comm_coul_inter_q_h.bcast(&nq_this, 1, i_comm_q);
                // Let processes out of comm_coul know what are going to be processed
                int pid_comm_coul_root_in_global = 0;
                if (comm_coul_h.is_initialized() && comm_coul_h.myid == 0)
                    pid_comm_coul_root_in_global = comm_h.myid;
                comm_h.allreduce(MPI_IN_PLACE, &pid_comm_coul_root_in_global, 1, MPI_MAX);
                if (n_comms_q_min == 0) comm_h.bcast(&nq_this, 1, pid_comm_coul_root_in_global);
                if (nq_this < 1) continue;
                ofs_myid << "vq_src nq_this for i_inter_q " << i_comm_q << " = " << nq_this << endl; 
                ofs_myid << "Array desc for source         : " << desc << endl;
                ofs_myid << "Array desc for dest (desc_abf): " << desc_abf.info_desc() << endl;
                std::vector<double> qs(nq_this * 3);
                map<int, Matz> mat_q;
                Matz mat_comm(0, 0, major_comm);
                if (i_comm_q == comm_coul_inter_q_h.myid)
                {
                    int iq = 0;
                    for (const auto &[q, mat]: vq_src)
                    {
                        qs[iq * 3] = q.x;
                        qs[iq * 3 + 1] = q.y;
                        qs[iq * 3 + 2] = q.z;
                        mat_q[iq] = mat;
                        iq++;
                    }
                }
                else
                {
                    // Dummy matrix
                    for (int iq = 0; iq < nq_this; iq++)
                        mat_q[iq] = Matz(0, 0, major_comm);
                }
                comm_coul_inter_q_h.bcast(qs.data(), nq_this * 3, i_comm_q);
                if (n_comms_q_min == 0) comm_h.bcast(qs.data(), nq_this * 3, pid_comm_coul_root_in_global);
                ofs_myid << "broadcasting coulomb matrices on q-points: " << qs << endl; 
                for (int iq = 0; iq < nq_this; iq++)
                {
                    const auto &mat_lo_src = mat_q.at(iq);
                    ScalapackConnector::pgemr2d(n_aux, n_aux, mat_lo_src.ptr(), 1, 1, desc.data(),
                                                mat_local.ptr(), 1, 1, desc_abf.desc,
                                                blacs_h.ictxt);
                    Vector3_Order<double> q{qs[iq*3], qs[iq*3+1], qs[iq*3+2]};
                    vq_dst[q] = mat_local.copy();
                    mat_q.erase(iq);
                }
            }
            // Do the same for cut Coulomb
            // int nq_cut = vq_cut_block_loc.size();
        }
    }
    else
    {
        vq_redist = std::move(vq_block_loc);
        vq_cut_redist = std::move(vq_cut_block_loc);
    }
    // Till now, every process has Coulomb matrices at all q-points, but only one BLACS block
    // Consistency check
    ofs_myid << "vq_redist size after pzgemr2d    : " << vq_redist.size() << endl;
    ofs_myid << "vq_cut_redist size after pzgemr2d: " << vq_cut_redist.size() << endl;
    assert(vq_redist.size() == 0 || vq_redist.size() == this->pbc.klist_ibz.size());
    assert(vq_cut_redist.size() == 0 || vq_cut_redist.size() == this->pbc.klist_ibz.size());

    // Step 2: redistribute BLACS 2D layout to atom-pair layout, and copy from column-major Matz to ComplexMatrix
    IndexScheduler sched;
    ofs_myid << "atpairs_unique_all:" << this->atpairs_unique_all << endl;
    sched.init(this->atpairs_unique_all, this->basis_aux, this->basis_aux, desc_abf, false);
    std::vector<std::pair<std::map<Vector3_Order<double>, Matz> &, atpair_k_cplx_mat_t &>>
        vq_src_dst_pair{{vq_redist, vq}, {vq_cut_redist, vq_cut}};
    // both bare and cut coulomb
    for (auto &[vq_src, vq_dst]: vq_src_dst_pair)
    {
        for (auto it_q = vq_src.begin(); it_q != vq_src.end();)
        {
            const auto &q = it_q->first;
            const auto &mat_loc = it_q->second;
            auto IJmap = get_ap_map_from_blacs_dist_scheduler(mat_loc, sched, basis_aux, basis_aux, desc_abf);
            for (auto it_ap = IJmap.begin(); it_ap != IJmap.end();)
            {
                const auto IJ = it_ap->first;
                const auto I = IJ.first;
                const auto n_I = basis_aux[I];
                const auto J = IJ.second;
                const auto n_J = basis_aux[J];
                auto &cmat_new = vq_dst[I][J][q];
                cmat_new = std::make_shared<ComplexMatrix>(n_I, n_J);
                const size_t n = n_I * n_J;
                auto &matz = it_ap->second;
                assert(as_size(matz.nr()) == n_I);
                assert(as_size(matz.nc()) == n_J);
                matz.swap_to_row_major();
                memcpy(cmat_new->c, matz.ptr(), n * sizeof(cplxdb));
                assert(as_size(matz.nr()) == n_I);
                assert(as_size(matz.nc()) == n_J);
                assert(as_size(matz.size()) == n);
                it_ap = IJmap.erase(it_ap);
            }
            it_q = vq_src.erase(it_q);
        }
    }

    coul_blacs2ap_redistributed_ = true;
    profiler.stop(__FUNCTION__);
}

void Dataset::finalize_comm_blacs_coul()
{
    if (!comm_blacs_coul_initialized_) return;
    global::profiler.start(__FUNCTION__);
    desc_coul_intra_q.reset_handler();
    blacs_coul_intra_q_h.reset_comm();
    comm_coul_inter_q_h.free_comm();
    comm_coul_intra_q_h.free_comm();
    comm_coul_h.free_comm();
    comm_blacs_coul_initialized_ = false;
    global::profiler.stop(__FUNCTION__);
}

}
