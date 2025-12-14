#include "exx.h"

#include <omp.h>

#include "../io/global_io.h"
#include "../io/stl_io_helper.h"
#include "../math/lapack_connector.h"
#include "../math/utils_matrix_m_mpi.h"
#include "../math/vector3_order.h"
#include "../mpi/global_mpi.h"
#include "../utils/constants.h"
#include "../utils/libri_utils.h"
#include "../utils/profiler.h"
#include "atomic_basis.h"
#include "meanfield_mpi.h"
#include "geometry.h"
#include "librpa_enums.h"
// #include "params.h"
#include "pbc.h"
#include "utils_atomic_basis_blacs.h"
#ifdef LIBRPA_USE_LIBRI
#include <RI/physics/Exx.h>
#include <RI/ri/Cell_Nearest.h>
#else
#include "../utils/libri_stub.h"
#endif

namespace librpa_int
{

Exx::Exx(const MeanField &mf_in,
         const AtomicBasis &atbasis_wfc_in,
         const PeriodicBoundaryData &pbc_in,
         const MpiCommHandler &comm_h_in,
         bool is_mf_eigvec_k_distributed)
    : mf(mf_in), atbasis_wfc(atbasis_wfc_in), pbc(pbc_in), comm_h(comm_h_in)
{
    is_mf_eigvec_k_distributed_ = is_mf_eigvec_k_distributed;
    is_rspace_built_ = false;
    is_kspace_built_ = false;
    is_rspace_redist_for_KS_ = false;
    is_rspace_redist_blacs_ = false;

    // Runtime options
    libri_threshold_C = 0.0;
    libri_threshold_V = 0.0;
    libri_threshold_D = 0.0;
};

static ComplexMatrix extract_dmat_cplx_R_IJblock(const ComplexMatrix& dmat_cplx, const AtomicBasis &ab_wfc, const atom_t& I, const atom_t& J)
{
    const auto I_num = ab_wfc.get_atom_nb(I);
    const auto J_num = ab_wfc.get_atom_nb(J);
    ComplexMatrix dmat_cplx_IJR(I_num, J_num);
    for (size_t i = 0; i != I_num; i++)
    {
        size_t i_glo = ab_wfc.get_global_index(I, i);
        for (size_t j = 0; j != J_num; j++)
        {
            size_t j_glo = ab_wfc.get_global_index(J, j);
            dmat_cplx_IJR(i, j) = dmat_cplx(i_glo, j_glo);
        }
    }
    return dmat_cplx_IJR;
}


static void warn_dmat_IJR_nonzero_imag(const ComplexMatrix& dmat_cplx, const int& ispin, const atom_t& I, const atom_t& J, const Vector3_Order<int> R)
{
    if (dmat_cplx.get_max_abs_imag() > 1e-2)
        global::lib_printf("Warning: complex-valued density matrix, spin %d IJR %zu %zu (%d, %d, %d)\n", ispin, I, J, R.x, R.y, R.z);
}


#ifdef LIBRPA_USE_LIBRI
static void build_dmat_libri_serial(
    const MeanField &mf,
    const AtomicBasis &atbasis_wfc,
    int ispin,
    const vector<Vector3_Order<double>> &kfrac_list,
    const std::vector<std::pair<atpair_t, Vector3_Order<int>>> IJRs,
    std::map<int, std::map<std::pair<int,std::array<int,3>>,RI::Tensor<double>>> &dmat_libri)
{
    // for (const auto &[ik, mat]: mf.get_eigenvectors().at(ispin))
    // {
    //     global::ofs_myid << ik << std::endl;
    //     print_complex_matrix("eigenvector", mat, global::ofs_myid, true);
    // }
    global::profiler.start("exx_build_dmat_libri_serial");
    std::map<Vector3_Order<int>, std::vector<atpair_t>> map_R_IJs;
    for (const auto &IJR: IJRs)
    {
        const auto &R = IJR.second;
        map_R_IJs[R].push_back(IJR.first);
    }
    for (const auto &R_IJs: map_R_IJs)
    {
        const auto &R = R_IJs.first;
        const auto &IJs = R_IJs.second;
        std::array<int,3> Ra{R.x,R.y,R.z};
        const auto dmat_cplx = mf.get_dmat_cplx_R(ispin, kfrac_list, R);
        // global::ofs_myid << R << std::endl;
        // print_complex_matrix("dmat_cplx[R]", dmat_cplx, global::ofs_myid, true);
        omp_lock_t dmat_lock;
        omp_init_lock(&dmat_lock);
#pragma omp parallel for schedule(dynamic)
        for (const auto &IJ: IJs)
        {
            const auto &I = IJ.first;
            const auto &J = IJ.second;
            const auto dmat_IJR = extract_dmat_cplx_R_IJblock(dmat_cplx, atbasis_wfc, I, J);
            warn_dmat_IJR_nonzero_imag(dmat_IJR, ispin, I, J, R);
            std::valarray<double> dmat_va(dmat_IJR.real().c, dmat_IJR.size);
            auto pdmat = std::make_shared<std::valarray<double>>();
            *pdmat = dmat_va;
            omp_set_lock(&dmat_lock);
            dmat_libri[I][{J, Ra}] = RI::Tensor<double>({size_t(dmat_IJR.nr), size_t(dmat_IJR.nc)}, pdmat);
            omp_unset_lock(&dmat_lock);
        }
#pragma omp barrier
        omp_destroy_lock(&dmat_lock);
    }
    global::profiler.stop("exx_build_dmat_libri_serial");
}

static void build_dmat_libri_kpara(
    const MeanField &mf,
    const MpiCommHandler &comm_h,
    const AtomicBasis &atbasis_wfc,
    int ispin,
    const vector<Vector3_Order<double>> &kfrac_list,
    const std::vector<std::pair<atpair_t, Vector3_Order<int>>> IJRs,
    std::map<int, std::map<std::pair<int,std::array<int,3>>,RI::Tensor<double>>> &dmat_libri)
{
    global::profiler.start("exx_build_dmat_libri_kpara");

    std::map<Vector3_Order<int>, std::vector<atpair_t>> map_R_IJs;
    for (const auto &[IJ, R]: IJRs)
    {
        map_R_IJs[R].emplace_back(IJ);
    }
    std::vector<Vector3_Order<int>> Rs_this;
    for (const auto &[R, _]: map_R_IJs)
        Rs_this.emplace_back(R);
    const int n_Rs_this = map_R_IJs.size();
    int n_Rs_max = n_Rs_this;
    MPI_Allreduce(MPI_IN_PLACE, &n_Rs_max, 1, MPI_INT, MPI_MAX, comm_h.comm);
    const auto dmat_Rs_cplx = get_dmat_cplx_Rs_kpara(ispin, mf, kfrac_list, Rs_this, comm_h);
    // global::ofs_myid << kfrac_list << std::endl;
    for (const auto &R_dmat_cplx: dmat_Rs_cplx)
    {
        const auto &R = R_dmat_cplx.first;
        const auto &dmat_cplx = R_dmat_cplx.second;
        // global::ofs_myid << R << std::endl;
        // print_complex_matrix("test", dmat_cplx, global::ofs_myid, true);
        std::array<int,3> Ra{R.x,R.y,R.z};
        const auto &map_IJs = map_R_IJs.at(R);
        omp_lock_t dmat_lock;
        omp_init_lock(&dmat_lock);
#pragma omp parallel for schedule(dynamic)
        for (const auto &IJ: map_IJs)
        {
            const auto &I = IJ.first;
            const auto &J = IJ.second;
            const auto dm_block = extract_dmat_cplx_R_IJblock(dmat_cplx, atbasis_wfc, I, J);
            warn_dmat_IJR_nonzero_imag(dm_block, ispin, I, J, R);
            std::valarray<double> dmat_va(dm_block.real().c, dm_block.size);
            auto pdmat = std::make_shared<std::valarray<double>>();
            *pdmat = dmat_va;
            omp_set_lock(&dmat_lock);
            dmat_libri[I][{J, Ra}] = RI::Tensor<double>({size_t(dm_block.nr), size_t(dm_block.nc)}, pdmat);
            omp_unset_lock(&dmat_lock);
        }
#pragma omp barrier
        omp_destroy_lock(&dmat_lock);
    }
    global::profiler.stop("exx_build_dmat_libri_kpara");
}
#endif

void Exx::build(const LibrpaParallelRouting routing,
                const AtomicBasis &atbasis_abf, const Cs_LRI &Cs,
                const atpair_R_mat_t &coul_mat)
{
    using std::endl;
    using global::profiler;
    using global::ofs_myid;

    assert(routing == LibrpaParallelRouting::LIBRI);

    if (this->is_rspace_built_)
    {
        return;
    }

    const auto &n_spins = this->mf.get_n_spins();
    const auto &ab_wfc = this->atbasis_wfc;
    const int n_atoms = ab_wfc.n_atoms;
    const auto &Rlist = this->pbc.Rlist;

#ifdef LIBRPA_USE_LIBRI
    if (comm_h.is_root())
    {
        global::lib_printf("Computing EXX orbital energy using LibRI\n");
    }
    comm_h.barrier();

    RI::Exx<int, int, 3, double> exx_libri;
    map<int,std::array<double,3>> atoms_pos;
    for (int i = 0; i < n_atoms; i++)
        atoms_pos.insert(pair<int, std::array<double, 3>>{i, {0, 0, 0}});

    exx_libri.set_parallel(comm_h.comm, atoms_pos, this->pbc.latvec_array, this->pbc.period_array);

    // Initialize Cs libRI container on each process
    // Note: we use different treatment in different routings
    //     R-tau routing:
    //         Each process has a full Cs copy.
    //         Thus in each process we only pass a few to LibRI container.
    //     atom-pair routing:
    //         Cs is already distributed across all processes.
    //         Pass the all Cs to libRI container.

    profiler.start("build_real_space_exx_1", "Prepare C libRI object");
    ofs_myid << "Number of Cs keys: " << get_num_keys(Cs.data_libri) << endl;
    // print_keys(global::ofs_myid, Cs.data_libri);
    global::ofs_myid << "Before set_Cs" << std::endl;
    // global::ofs_myid << Cs.data_libri << std::endl;
    exx_libri.set_Cs(Cs.data_libri, libri_threshold_C);
    global::ofs_myid << "Finish  set_Cs" << std::endl;
    // exx_libri.set_Cs({}, libri_threshold_C);
    profiler.stop("build_real_space_exx_1");
    ofs_myid << "Finished setup Cs for EXX\n";
    std::flush(global::ofs_myid);

    // initialize Coulomb matrix
    global::profiler.start("build_real_space_exx_2", "Prepare V libRI object");
    std::map<int, std::map<std::pair<int,std::array<int,3>>, RI::Tensor<double>>> V_libri;
    global::profiler.start("build_real_space_exx_2_1");

    if (routing == LibrpaParallelRouting::RTAU)
    {
        // Full Coulomb case, have to re-distribute
        for (auto IJR: dispatch_vector_prod(get_atom_pair(coul_mat), Rlist, comm_h.myid, comm_h.nprocs, true, true))
        {
            const auto I = IJR.first.first;
            const auto J = IJR.first.second;
            const auto R = IJR.second;
            const auto& VIJR = coul_mat.at(I).at(J).at(R);
            // debug
            // printf("I J R %zu %zu %d %d %d, max(V) %f\n", I, J, R.x, R.y, R.z, VIJR->max());
            std::array<int,3> Ra{R.x,R.y,R.z};
            std::valarray<double> VIJR_va(VIJR->c, VIJR->size);
            auto pv = std::make_shared<std::valarray<double>>();
            *pv = VIJR_va;
            V_libri[I][{J, Ra}] = RI::Tensor<double>({size_t(VIJR->nr), size_t(VIJR->nc)}, pv);
        }
    }
    else
    {
        for (const auto &I_JRV: coul_mat)
        {
            const auto I = I_JRV.first;
            for (const auto &J_RV: I_JRV.second)
            {
                const auto J = J_RV.first;
                for (const auto &R_V: J_RV.second)
                {
                    const auto &R = R_V.first;
                    const auto &V = R_V.second;
                    std::array<int,3> Ra{R.x,R.y,R.z};
                    std::valarray<double> VIJR_va(V->c, V->size);
                    auto pv = std::make_shared<std::valarray<double>>();
                    *pv = VIJR_va;
                    V_libri[I][{J, Ra}] = RI::Tensor<double>({size_t(V->nr), size_t(V->nc)}, pv);
                }
            }
        }
    }
    global::profiler.stop("build_real_space_exx_2_1");
    global::ofs_myid << "Number of V keys: " << get_num_keys(V_libri) << endl;
    global::profiler.start("build_real_space_exx_2_2");
    exx_libri.set_Vs(V_libri, libri_threshold_V);
    V_libri.clear();
    global::profiler.stop("build_real_space_exx_2_2");
    global::profiler.stop("build_real_space_exx_2");
    global::lib_printf("Task %4d: V setup for EXX\n", comm_h.myid);
    // cout << V_libri << endl;

    // initialize density matrix
    vector<atpair_t> atpair_dmat;
    for (atom_t I = 0; I < as_atom(n_atoms); I++)
        for (atom_t J = 0; J < as_atom(n_atoms); J++)
            atpair_dmat.push_back({I, J});
    const auto dmat_IJRs_local = dispatch_vector_prod(atpair_dmat, Rlist, comm_h.myid, comm_h.nprocs, true, true);

    for (auto isp = 0; isp != n_spins; isp++)
    {
        global::profiler.start("build_real_space_exx_3", "Prepare DM libRI object");
        std::map<int, std::map<std::pair<int,std::array<int,3>>,RI::Tensor<double>>> dmat_libri;
        if (is_mf_eigvec_k_distributed_)
        {
            build_dmat_libri_kpara(mf, comm_h, atbasis_wfc, isp, this->pbc.kfrac_list, dmat_IJRs_local, dmat_libri);
        }
        else
        {
            build_dmat_libri_serial(mf, atbasis_wfc, isp, this->pbc.kfrac_list, dmat_IJRs_local, dmat_libri);
        }
        global::ofs_myid << "Number of Dmat keys: " << get_num_keys(dmat_libri) << "\n";
        // global::ofs_myid << dmat_libri << std::endl;
        // print_keys(global::ofs_myid, dmat_libri);
        exx_libri.set_Ds(dmat_libri, libri_threshold_D);
        global::profiler.stop("build_real_space_exx_3");
        global::lib_printf("Task %4d: DM setup for EXX\n", comm_h.myid);

        global::profiler.start("build_real_space_exx_4", "Call libRI Hexx calculation");
        exx_libri.cal_Hs();
        global::profiler.stop("build_real_space_exx_4");

        global::lib_printf("Task %4d: cal_Hs elapsed time: %f\n", comm_h.myid, global::profiler.get_wall_time_last("build_real_space_exx_4"));
        global::ofs_myid << "Number of exx_libri.Hs keys: " << get_num_keys(exx_libri.Hs) << "\n";
        // print_keys(global::ofs_myid, exx_libri.Hs);
        // ofs_myid << "exx_libri.Hs:\n" << exx_libri.Hs << endl;

        for (const auto &I_JR_exx: exx_libri.Hs)
        {
            const auto &I = I_JR_exx.first;
            const auto &n_I = ab_wfc.get_atom_nb(I);
            for (const auto &JR_exx: I_JR_exx.second)
            {
                const auto &J = JR_exx.first.first;
                const auto &n_J = ab_wfc.get_atom_nb(J);
                const auto &Ra = JR_exx.first.second;
                const auto R = Vector3_Order<int>{Ra[0], Ra[1], Ra[2]};
                Matd exx_temp(n_I, n_J, JR_exx.second.ptr(), MAJOR::ROW);
                this->exx_IJR[isp][I][J][R] = exx_temp;
            }
        }
    }
    // debug, print the Hexx matrices
    // for (const auto& isp_IJkH: this->Hexx)
    // {
    //     const auto& isp = isp_IJkH.first;
    //     for (const auto& I_JkH: isp_IJkH.second)
    //     {
    //         const auto& I = I_JkH.first;
    //         for (const auto& J_kH: I_JkH.second)
    //         {
    //             const auto& J = J_kH.first;
    //             for (const auto& k_H: J_kH.second)
    //             {
    //                 cout << isp << " " << I << " " << J << " {" << k_H.first << "} whole size: " << k_H.second->size << endl;
    //                 print_complex_matrix("", *k_H.second);
    //             }
    //         }
    //     }
    // }

#else
    if (comm_h.is_root())
    {
        global::lib_printf("Error: trying build EXX orbital energy with LibRI, but the program is not compiled against LibRI\n");
    }
    throw std::logic_error("compilation");
    comm_h.barrier();
#endif

    is_rspace_built_= true;
}

void Exx::build_KS(const std::map<int, std::map<int, ComplexMatrix>> &wfc_target,
                   const std::vector<Vector3_Order<double>> &kfrac_target,
                   const Atoms &geometry)
{
    throw LIBRPA_RUNTIME_ERROR("not implemented");
}

void Exx::build_KS_blacs(const std::map<int, std::map<int, ComplexMatrix>> &wfc_target,
                         const std::vector<Vector3_Order<double>> &kfrac_target,
                         const Atoms &geometry,
                         const BlacsCtxtHandler &blacs_ctxt_h)
{
    using RI::Communicate_Tensors_Map_Judge::comm_map2_first;

    // Ensure the communicator of BLACS context is the same as that parsed when constructed
    assert(blacs_ctxt_h.comm() == this->comm_h.comm);
    assert(this->is_rspace_built_);

    // Reset k-space matrices built from last call
    if (this->is_kspace_built_)
    {
        global::lib_printf("Warning: reset EXX k-space matrices\n");
        this->reset_kspace();
    }

    const auto& n_aos = this->mf.get_n_aos();
    const auto& n_spins = this->mf.get_n_spins();
    const auto& n_bands = this->mf.get_n_bands();

    // prepare scalapack array descriptors
    ArrayDesc desc_nao_nao(blacs_ctxt_h);
    ArrayDesc desc_nband_nao(blacs_ctxt_h);
    ArrayDesc desc_nband_nband(blacs_ctxt_h);

    // For communication of eigenvectors and final KS matrices
    ArrayDesc desc_nao_nband_fb(blacs_ctxt_h);  // emulate ComplexMatrix storage
    ArrayDesc desc_nband_nband_fb(blacs_ctxt_h);

    desc_nao_nao.init_1b1p(n_aos, n_aos, 0, 0);
    desc_nband_nao.init_1b1p(n_bands, n_aos, 0, 0);
    desc_nband_nband.init_1b1p(n_bands, n_bands, 0, 0);

    // local 2D-block submatrices
    auto Hexx_nao_nao = init_local_mat<complex<double>>(desc_nao_nao, MAJOR::COL);
    auto Hexx_nband_nband = init_local_mat<complex<double>>(desc_nband_nband, MAJOR::COL);
    auto temp_nband_nao = init_local_mat<complex<double>>(desc_nband_nao, MAJOR::COL);

    if (!is_rspace_redist_for_KS_)
    {
        const auto set_IJ_naonao = get_necessary_IJ_from_block_2D(
            this->atbasis_wfc, this->atbasis_wfc, desc_nao_nao);
        const auto Iset_Jset = convert_IJset_to_Iset_Jset(set_IJ_naonao);

        for (int isp = 0; isp < n_spins; isp++)
        {
            // collect necessary data
            global::profiler.start("build_real_space_exx_5", "Collect Hexx IJ from world");
            std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<double>>> exx_is_tensor;

            auto it = this->exx_IJR.find(isp);
            if (it != this->exx_IJR.cend())
            {
                const auto &exx_is = it->second;
                for (const auto &[I, J_Rexx]: exx_is)
                {
                    const auto &n_I = this->atbasis_wfc.get_atom_nb(I);
                    for (const auto &[J, R_exx]: J_Rexx)
                    {
                        for (const auto &[R, mat]: R_exx)
                        {
                            const auto &n_J = this->atbasis_wfc.get_atom_nb(J);
                            const std::array<int, 3> Ra{R.x, R.y, R.z};
                            exx_is_tensor[I][{J, Ra}] = RI::Tensor<double>({n_I, n_J}, mat.sptr());
                        }
                    }
                }
            }
            // Collect the IJ pair of Hs with all R for Fourier transform
            auto exx_is_IJR_for_blacs = comm_map2_first(comm_h.comm, exx_is_tensor, Iset_Jset.first, Iset_Jset.second);
            exx_is_tensor.clear();
            if (exx_is_IJR_for_blacs.size() == 0) continue;
            // Now each process should have all Rs corresponding to atom-pairs that required for BLACS matrix operation
            // Swap with the original
            map<atom_t, map<atom_t, map<Vector3_Order<int>, Matd>>> exx_is_new;
            for (const auto &[I, JRmat]: exx_is_IJR_for_blacs)
            {
                const int n_I = this->atbasis_wfc.get_atom_nb(I);
                for (const auto &[JR, mat]: JRmat)
                {
                    const atom_t J = JR.first;
                    const int n_J = this->atbasis_wfc.get_atom_nb(J);
                    const auto &Ra = JR.second;
                    const Vector3_Order<int> R{Ra[0], Ra[1], Ra[2]};
                    exx_is_new[I][J][R] = Matd{n_I, n_J, mat.data, MAJOR::ROW};
                }
            }
            this->exx_IJR.at(isp).swap(exx_is_new);
            exx_is_IJR_for_blacs.clear();
        }
        is_rspace_redist_for_KS_ = true;
        is_rspace_redist_blacs_ = true;
    }
    else
    {
        if (!is_rspace_redist_blacs_)
            throw LIBRPA_RUNTIME_ERROR("");
    }

    for (int isp = 0; isp < n_spins; isp++)
    {
        this->exx_KS[isp] = {};
        this->Eexx[isp] = {};
        map<atom_t, map<atom_t, map<Vector3_Order<int>, Matd>>> exx_is_local;
        // Convert each <I,<J, R>> pair to the nearest neighbour to speed up later Fourier transform
        // while keep the accuracy in further band interpolation.
        // Reuse the cleared-up exx_I_JR_local object
        auto it = exx_IJR.find(isp);
        if (it != exx_IJR.cend())
        {
            if (geometry.is_frac_set())
            {
                const auto &coords_frac = geometry.coords_frac;
                const auto &period = this->pbc.period;

                for (const auto &[I, J_Rexx]: it->second)
                {
                    for (const auto &[J, R_exx]: J_Rexx)
                    {
                        for (const auto &[R, mat]: R_exx)
                        {
                            auto distsq = std::numeric_limits<double>::max();
                            Vector3<int> R_IJ;
                            Vector3<int> R_bvk;
                            for (int i = -1; i < 2; i++)
                            {
                                R_IJ.x = i * period.x + R.x;
                                for (int j = -1; j < 2; j++)
                                {
                                    R_IJ.y = j * period.y + R.x;
                                    for (int k = -1; k < 2; k++)
                                    {
                                        R_IJ.z = k * period.z + R.z;
                                        const auto diff =
                                            (Vector3<double>(coords_frac.at(I)[0], coords_frac.at(I)[1],
                                                            coords_frac.at(I)[2]) -
                                            Vector3<double>(coords_frac.at(J)[0], coords_frac.at(J)[1],
                                                            coords_frac.at(J)[2]) -
                                            Vector3<double>(R_IJ.x, R_IJ.y, R_IJ.z)) * this->pbc.latvec;
                                        const auto norm2 = diff.norm2();
                                        if (norm2 < distsq)
                                        {
                                            distsq = norm2;
                                            R_bvk = R_IJ;
                                        }
                                    }
                                }
                            }
                            exx_is_local[I][J][R_bvk] = mat;
                        }
                    }
                }
            }
            else
            {
                exx_is_local = it->second;
            }
        }
        global::profiler.stop("build_real_space_exx_5");
        global::lib_printf("Task %4d: tensor communicate elapsed time: %f\n", comm_h.myid, global::profiler.get_wall_time_last("build_real_space_exx_5"));

        if (is_mf_eigvec_k_distributed_)
        {
            // collect parsed eigenvectors all all processes
            // global::ofs_myid << "isp " << isp << std::endl;
            std::vector<int> nks_all(comm_h.nprocs, 0);
            std::vector<int> iks_local;
            auto it = wfc_target.find(isp);
            if (it != wfc_target.cend())
            {
                const auto &wfc_sp = it->second;
                nks_all[comm_h.myid] = wfc_sp.size();
                for (const auto &[ik, _]: wfc_sp)
                {
                    iks_local.emplace_back(ik);
                }
            }
            MPI_Allreduce(MPI_IN_PLACE, nks_all.data(), comm_h.nprocs, mpi_datatype<int>::value, MPI_SUM, comm_h.comm);
            const int nk_max = *std::max_element(nks_all.cbegin(), nks_all.cend());
            std::vector<int> iks_all(nk_max * comm_h.nprocs, 0);
            for (int i = 0; i < nks_all[comm_h.myid]; i++)
            {
                iks_all[comm_h.myid * nk_max + i] = iks_local[i];
            }
            MPI_Allreduce(MPI_IN_PLACE, iks_all.data(), comm_h.nprocs * nk_max, mpi_datatype<int>::value, MPI_SUM, comm_h.comm);

            for (int pid = 0; pid < comm_h.nprocs; pid++)
            {
                const int nk_this = nks_all[pid];
                if (nk_this == 0) continue;  // no eigenvector on this process
                auto [irsrc, icsrc] = blacs_ctxt_h.get_pcoord(pid);
                desc_nao_nband_fb.init(n_aos, n_bands, n_aos, n_bands, irsrc, icsrc);
                desc_nband_nband_fb.init(n_bands, n_bands, n_bands, n_bands, irsrc, icsrc);
                auto Hexx_nband_nband_fb = init_local_mat<complex<double>>(desc_nband_nband_fb, MAJOR::COL);
                for (int ik_this = 0; ik_this < nk_this; ik_this++)
                {
                    global::profiler.start("build_real_space_exx_6", "Hexx IJ -> 2D block");
                    Hexx_nao_nao.zero_out();
                    const int ik = iks_all[pid * nk_max + ik_this];
                    const auto& kfrac = kfrac_target[ik];
                    const std::function<complex<double>(const atom_t, const atom_t, const Vector3_Order<int> &)>
                        fourier = [kfrac](const atom_t I, const atom_t J, const Vector3_Order<int> &R)
                        {
                            const auto ang = (kfrac * R) * TWO_PI;
                            return complex<double>{std::cos(ang), std::sin(ang)};
                        };
                    collect_block_from_IJ_storage_matrix_transform(Hexx_nao_nao, desc_nao_nao,
                            this->atbasis_wfc, this->atbasis_wfc, fourier, exx_is_local);
                    global::profiler.stop("build_real_space_exx_6");

                    std::vector<complex<double>> dummy(1);
                    global::profiler.start("build_real_space_exx_7", "Rotate Hexx ij -> KS");
                    if (pid == comm_h.myid)
                    {
                        // global::ofs_myid << "isp " << isp << " ik " << ik << std::endl;
                        const auto &wfc_sk = wfc_target.at(isp).at(ik);
                        // processing the k-point eigenvector on this process
                        ScalapackConnector::pgemm_f('C', 'N', n_bands, n_aos, n_aos, 1.0,
                                                    wfc_sk.c, 1, 1, desc_nao_nband_fb.desc,
                                                    Hexx_nao_nao.ptr(), 1, 1, desc_nao_nao.desc,
                                                    0.0,
                                                    temp_nband_nao.ptr(), 1, 1, desc_nband_nao.desc);
                        ScalapackConnector::pgemm_f('N', 'N', n_bands, n_bands, n_aos, -1.0,
                                                    temp_nband_nao.ptr(), 1, 1, desc_nband_nao.desc,
                                                    wfc_sk.c, 1, 1, desc_nao_nband_fb.desc,
                                                    0.0,
                                                    Hexx_nband_nband_fb.ptr(), 1, 1, desc_nband_nband_fb.desc);
                        this->exx_KS[isp][ik] = Hexx_nband_nband_fb.copy();
                        // cout << "Hexx_nband_nband_fb isp " << isp  << " ik " << ik << endl << Hexx_nband_nband_fb;
                        for (int ib = 0; ib != n_bands; ib++)
                            this->Eexx[isp][ik][ib] = Hexx_nband_nband_fb(ib, ib).real();
                    }
                    else
                    {
                        // processing the k-point eigenvector at other processes
                        ScalapackConnector::pgemm_f('C', 'N', n_bands, n_aos, n_aos, 1.0,
                                                    dummy.data(), 1, 1, desc_nao_nband_fb.desc,
                                                    Hexx_nao_nao.ptr(), 1, 1, desc_nao_nao.desc,
                                                    0.0,
                                                    temp_nband_nao.ptr(), 1, 1, desc_nband_nao.desc);
                        ScalapackConnector::pgemm_f('N', 'N', n_bands, n_bands, n_aos, -1.0,
                                                    temp_nband_nao.ptr(), 1, 1, desc_nband_nao.desc,
                                                    dummy.data(), 1, 1, desc_nao_nband_fb.desc,
                                                    0.0,
                                                    Hexx_nband_nband_fb.ptr(), 1, 1, desc_nband_nband_fb.desc);
                    }
                    global::profiler.stop("build_real_space_exx_7");
                }
            }
        }
        else
        {
            // Everything will be collected to rank 0
            desc_nband_nband_fb.init(n_bands, n_bands, n_bands, n_bands, 0, 0);
            auto Hexx_nband_nband_fb = init_local_mat<complex<double>>(desc_nband_nband_fb, MAJOR::COL);

            // Each process have all KS eigenvectors, extract local blocks
            for (size_t ik = 0; ik < kfrac_target.size(); ik++)
            {
                global::profiler.start("build_real_space_exx_6", "Hexx IJ -> 2D block");
                Hexx_nao_nao.zero_out();
                const auto& kfrac = kfrac_target[ik];
                const std::function<complex<double>(const atom_t, const atom_t, const Vector3_Order<int> &)>
                    fourier = [kfrac](const atom_t I, const atom_t J, const Vector3_Order<int> &R)
                    {
                        const auto ang = (kfrac * R) * TWO_PI;
                        return complex<double>{std::cos(ang), std::sin(ang)};
                    };
                collect_block_from_IJ_storage_matrix_transform(Hexx_nao_nao, desc_nao_nao,
                        this->atbasis_wfc, this->atbasis_wfc, fourier, exx_is_local);
                global::profiler.stop("build_real_space_exx_6");
                // global::lib_printf("%s\n", str(Hexx_nao_nao).c_str());
                const auto &wfc_isp_k = wfc_target.at(isp).at(ik);
                blacs_ctxt_h.barrier();
                const auto wfc_block = get_local_mat(wfc_isp_k.c, MAJOR::ROW, desc_nband_nao, MAJOR::COL).conj();
                // global::lib_printf("%s\n", str(wfc_block).c_str());
                // global::lib_printf("%s\n", desc_nao_nao.info_desc().c_str());
                // global::lib_printf("%s\n", desc_nband_nao.info_desc().c_str());
                global::profiler.start("build_real_space_exx_7", "Rotate Hexx ij -> KS");
                ScalapackConnector::pgemm_f('N', 'N', n_bands, n_aos, n_aos, 1.0,
                                            wfc_block.ptr(), 1, 1, desc_nband_nao.desc,
                                            Hexx_nao_nao.ptr(), 1, 1, desc_nao_nao.desc,
                                            0.0,
                                            temp_nband_nao.ptr(), 1, 1, desc_nband_nao.desc);
                ScalapackConnector::pgemm_f('N', 'C', n_bands, n_bands, n_aos, -1.0,
                                            temp_nband_nao.ptr(), 1, 1, desc_nband_nao.desc,
                                            wfc_block.ptr(), 1, 1, desc_nband_nao.desc,
                                            0.0,
                                            Hexx_nband_nband.ptr(), 1, 1, desc_nband_nband.desc);
                global::profiler.stop("build_real_space_exx_7");

                // collect to master
                global::profiler.start("build_real_space_exx_8", "Collect Eexx to root process");
                ScalapackConnector::pgemr2d_f(n_bands, n_bands,
                                            Hexx_nband_nband.ptr(), 1, 1, desc_nband_nband.desc,
                                            Hexx_nband_nband_fb.ptr(), 1, 1, desc_nband_nband_fb.desc,
                                            desc_nband_nband_fb.ictxt());
                this->exx_KS[isp][ik] = Hexx_nband_nband_fb;
                // cout << "Hexx_nband_nband_fb isp " << isp  << " ik " << ik << endl << Hexx_nband_nband_fb;
                if (blacs_ctxt_h.myid == 0)
                {
                    for (int ib = 0; ib != n_bands; ib++)
                        this->Eexx[isp][ik][ib] = Hexx_nband_nband_fb(ib, ib).real();
                }
                global::profiler.stop("build_real_space_exx_8");
            }
        }
    }
}

void Exx::build_KS_kgrid()
{
    this->build_KS(this->mf.get_eigenvectors(), this->pbc.kfrac_list, {});
}

void Exx::build_KS_band(const std::map<int, std::map<int, ComplexMatrix>> &wfc_band,
                        const std::vector<Vector3_Order<double>> &kfrac_band, const Atoms &geometry)
{
    this->build_KS(wfc_band, kfrac_band, geometry);
}

void Exx::build_KS_kgrid_blacs(const BlacsCtxtHandler &blacs_ctxt_h)
{
    this->build_KS_blacs(this->mf.get_eigenvectors(), this->pbc.kfrac_list, {}, blacs_ctxt_h);
}

void Exx::build_KS_band_blacs(const std::map<int, std::map<int, ComplexMatrix>> &wfc_band,
                              const std::vector<Vector3_Order<double>> &kfrac_band, const Atoms &geometry,
                              const BlacsCtxtHandler &blacs_ctxt_h)
{
    this->build_KS_blacs(wfc_band, kfrac_band, geometry, blacs_ctxt_h);
}

void Exx::reset_rspace()
{
    this->exx_IJR.clear();
    this->is_rspace_built_ = false;
}

void Exx::reset_kspace()
{
    this->exx_KS.clear();
    this->Eexx.clear();
    this->is_kspace_built_ = false;
}

} /* end of namespace librpa_int */
