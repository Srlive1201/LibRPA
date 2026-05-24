#include "gw.h"

// Public API headers
#include "librpa_enums.h"

#include <cstddef>
#include <fstream>
#include <functional>
#include <map>
#include <sstream>
#include <vector>

#include "../io/fs.h"
#include "../io/global_io.h"
#include "../io/output_aims.h"
#include "../io/stl_io_helper.h"
// #include "../math/utils_matrix_m.h"
#include "../math/utils_matrix_m_mpi.h"
#include "../mpi/global_mpi.h"
#include "../utils/constants.h"
#include "../utils/error.h"
#include "../utils/libri_utils.h"
#include "../utils/profiler.h"
#include "../utils/utils_mem.h"
#include "atom.h"
#include "atomic_basis.h"
#include "epsilon.h"
#include "geometry.h"
#include "meanfield_mpi.h"
#include "pbc.h"
#include "ri.h"
#include "utils_atomic_basis_blacs.h"
#include "utils_complexmatrix.h"

#ifdef LIBRPA_USE_LIBRI
#include <RI/global/Tensor.h>
#include <RI/physics/GW.h>
using RI::Tensor;
using RI::Communicate_Tensors_Map_Judge::comm_map2_first;
#endif

namespace librpa_int
{

G0W0::G0W0(const MeanField &mf_in, const AtomicBasis &atbasis_wfc_in,
           const PeriodicBoundaryData &pbc_in, const TFGrids &tfg_in,
           const MpiCommHandler &comm_h_in, bool is_mf_eigvec_k_distributed)
    : mf(mf_in),
      atbasis_wfc(atbasis_wfc_in),
      pbc(pbc_in),
      tfg(tfg_in), comm_h(comm_h_in)
{
    comm_h.check_initialized();

    is_mf_eigvec_k_distributed_ = is_mf_eigvec_k_distributed;
    is_rspace_built_ = false;
    is_kspace_built_ = false;
    is_rspace_redist_for_KS_ = false;
    is_rspace_redist_blacs_ = false;

    // Public runtime options
    libri_threshold_C = 0.0;
    libri_threshold_Wc = 0.0;
    libri_threshold_G = 0.0;
    output_dir = "./";  // POSIX
    output_sigc_mat = false;
    output_sigc_mat_rt = false;
    output_sigc_mat_rf = false;
}

void G0W0::reset_rspace()
{
    sigc_is_f_IJ_R.clear();
    is_rspace_built_ = false;
    is_rspace_redist_for_KS_ = false;
    is_rspace_redist_blacs_ = false;
}

void G0W0::reset_kspace()
{
    sigc_is_ik_f_KS.clear(); is_kspace_built_ = false;
}

#ifdef LIBRPA_USE_LIBRI
template <typename Tdata>
static void build_gf_libri_kpara(
    const MeanField &mf,
    const MpiCommHandler &comm_h,
    const AtomicBasis &atbasis_wfc,
    int ispin, int ispinor_bra, int ispinor_ket,
    const vector<Vector3_Order<double>> &kfrac_list,
    const std::vector<double> &taus,
    const std::vector<std::pair<atpair_t, Vector3_Order<int>>> IJRs,
    std::map<double, std::map<int, std::map<std::pair<int,std::array<int,3>>,RI::Tensor<Tdata>>>> &tau_gf_libri)
{
    global::profiler.start("g0w0_build_gf_libri_kpara");
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
    comm_h.allreduce(MPI_IN_PLACE, &n_Rs_max, 1, MPI_MAX);
    auto gf_taus_Rs_cplx = get_gf_cplx_imagtimes_Rs_kpara(ispin, ispinor_bra, ispinor_ket, mf, kfrac_list, taus, Rs_this, comm_h);
    // global::ofs_myid << "gf_taus_Rs_cplx " << gf_taus_Rs_cplx << std::endl;
    for (const auto &tau_gf_R_cplx: gf_taus_Rs_cplx)
    {
        const auto &tau = tau_gf_R_cplx.first;
        const auto &gf_R_cplx = tau_gf_R_cplx.second;
        for (const auto &R_gf_cplx: gf_R_cplx)
        {
            const auto &R = R_gf_cplx.first;
            const auto &gf_cplx = R_gf_cplx.second;
            std::array<int,3> Ra{R.x,R.y,R.z};
            const auto &map_IJs = map_R_IJs.at(R);
            omp_lock_t gf_lock;
            omp_init_lock(&gf_lock);
#pragma omp parallel for schedule(dynamic)
            for (const auto &IJ: map_IJs)
            {
                const auto &I = IJ.first;
                const auto &J = IJ.second;
                const auto gf_block = get_ap_block_from_global(gf_cplx, IJ, atbasis_wfc, atbasis_wfc);
                if constexpr (is_complex<Tdata>())
                {
                    auto p_gfmat = std::make_shared<std::valarray<Tdata>>(gf_block.c, gf_block.size);
                    omp_set_lock(&gf_lock);
                    tau_gf_libri[tau][I][{J, Ra}] = RI::Tensor<Tdata>({as_size(gf_block.nr), as_size(gf_block.nc)}, p_gfmat);
                    omp_unset_lock(&gf_lock);
                }
                else
                {
                    auto p_gfmat = std::make_shared<std::valarray<Tdata>>(gf_block.real().c, gf_block.size);
                    omp_set_lock(&gf_lock);
                    tau_gf_libri[tau][I][{J, Ra}] = RI::Tensor<Tdata>({as_size(gf_block.nr), as_size(gf_block.nc)}, p_gfmat);
                    omp_unset_lock(&gf_lock);
                }
            }
#pragma omp barrier
            omp_destroy_lock(&gf_lock);
        }
        // global::ofs_myid << R << std::endl;
        // print_complex_matrix("test", dmat_cplx, global::ofs_myid, true);
    }

    global::profiler.stop("g0w0_build_gf_libri_kpara");
}

template <typename Tdata>
static void build_gf_libri_kserial(
    const MeanField &mf,
    const AtomicBasis &atbasis_wfc,
    int ispin, int ispinor_bra, int ispinor_ket,
    const vector<Vector3_Order<double>> &kfrac_list,
    const std::vector<double> &taus,
    const std::vector<std::pair<atpair_t, Vector3_Order<int>>> IJRs,
    std::map<double, std::map<int, std::map<std::pair<int,std::array<int,3>>,RI::Tensor<Tdata>>>> &tau_gf_libri)
{
    global::profiler.start("g0w0_build_gf_libri_kserial");
    std::set<Vector3_Order<int>> Rs_local;
    for (const auto &IJR: IJRs)
    {
        Rs_local.insert(IJR.second);
    }
    auto gf = mf.get_gf_cplx_imagtimes_Rs(ispin, ispinor_bra, ispinor_ket, kfrac_list, taus, {Rs_local.cbegin(), Rs_local.cend()});
    // global::ofs_myid << "gf " << gf << std::endl;
    tau_gf_libri.clear();
    // TODO: enable threading below
    for (auto t: taus)
    {
        tau_gf_libri[t] = {};
    }
    for (const auto &IJR: IJRs)
    {
        const auto &IJ = IJR.first;
        const auto &I = IJ.first;
        const auto &n_I = atbasis_wfc.get_atom_nb(I);
        const auto &J = IJ.second;
        const auto &n_J = atbasis_wfc.get_atom_nb(J);
        const auto &R = IJR.second;
        for (auto t: taus)
        {
            // skip when this process does not have any local GF data, i.e. Rs_local is empty
            if (gf.count(t) == 0 || gf.at(t).count(R) == 0) continue;
            const auto &gf_global = gf.at(t).at(R);
            std::shared_ptr<std::valarray<Tdata>> mat_ptr = std::make_shared<std::valarray<Tdata>>(n_I * n_J);
            if constexpr (is_complex<Tdata>())
            {
                for (size_t i = 0; i != n_I; i++)
                    for (size_t j = 0; j != n_J; j++)
                    {
                        (*mat_ptr)[i * n_J + j] = gf_global(atbasis_wfc.get_global_index(I, i),
                                                            atbasis_wfc.get_global_index(J, j));
                    }
            }
            else
            {
                for (size_t i = 0; i != n_I; i++)
                    for (size_t j = 0; j != n_J; j++)
                    {
                        (*mat_ptr)[i * n_J + j] = gf_global(atbasis_wfc.get_global_index(I, i),
                                                            atbasis_wfc.get_global_index(J, j)).real();
                    }
            }
            tau_gf_libri[t][as_int(I)][{as_int(J), {R.x, R.y, R.z}}] = RI::Tensor<Tdata>({n_I, n_J}, mat_ptr);
        }
    }
    gf.clear();
    global::profiler.stop("g0w0_build_gf_libri_kserial");
}
#endif

void G0W0::build_spacetime(
    const LibrpaParallelRouting parallel_routing,
    const AtomicBasis& atbasis_abf,
    const Cs_LRI &LRI_Cs,
    std::map<double, std::map<Vector3_Order<double>, Matz>> &Wc_freq_q,
    const ArrayDesc &ad_Wc)
{
    using global::profiler;
    using global::lib_printf_root;
    using global::ofs_myid;

    if (parallel_routing != LibrpaParallelRouting::LIBRI)
    {
        comm_h.barrier();
        throw LIBRPA_RUNTIME_ERROR("not implemented");
    }

    lib_printf_root("Calculating correlation self-energy by space-time method\n");
    if (!tfg.has_time_grids())
    {
        lib_printf_root("Parsed time-frequency object do not have time grids, exiting\n");
        comm_h.barrier();
        throw LIBRPA_RUNTIME_ERROR("input TFGrids object has no time grids");
    }
    comm_h.barrier();

    const bool use_complex_tensor = mf.get_n_spinor() > 1;
    if (use_complex_tensor)
        lib_printf_root("Using complex tensor\n");
    else
        lib_printf_root("Using real tensor\n");

    assert(ad_Wc.initialized());
    // assert(blacs_sigc_h.initialized());
    // blacs_sigc_h_ = blacs_sigc_h;

    const int natom = this->atbasis_wfc.n_atoms;

    std::ofstream ofs_sigmac_r;

#ifndef LIBRPA_USE_LIBRI
    if (comm_h.myid == 0)
    {
        std::cout << "LIBRA::G0W0::build_spacetime is only implemented on top of LibRI" << std::endl;
        std::cout << "Please recompiler LibRPA with -DUSE_LIBRI and configure include path" << std::endl;
    }
    comm_h.barrier();
    throw LIBRPA_RUNTIME_ERROR("compilation");
#else
    // Transform from frequency/reciprocal to time/real-space
    profiler.start("g0w0_build_spacetime_ct_ft_wc", "Tranform Wc (q,w) -> (R,t)");
    profiler.start("g0w0_build_spacetime_ct_ft_real_work", "Perform transformation");
    // global::ofs_myid << Wc_freq_q.at(tfg.get_freq_nodes()[0]) << endl;
    // for (const auto &[freq, qWc]: Wc_freq_q)
    // {
    //     int iomega = tfg.get_freq_index(freq);
    //     for (const auto &[q, Wc]: qWc)
    //     {
    //         std::stringstream ss;
    //         ss << "Wc_omega_q"
    //             << "_iw_" << std::setfill('0') << std::setw(5) << iomega
    //             << "_iq_" << std::setfill('0') << std::setw(5) << pbc.get_k_index_full(q) << ".csc";
    //         write_matrix_elsi_csc_parallel(ss.str(), Wc, ad_Wc);
    //     }
    // }
    auto Wc_tau_R_blacs = CT_FT_Wc_freq_q(comm_h, Wc_freq_q, pbc, tfg, true);
    release_free_mem();
    profiler.stop("g0w0_build_spacetime_ct_ft_real_work");

    // Build the Wc LibRI object by redistributing the BLACS submatrices
    typedef std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<double>>> dtensor_map;
    typedef std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<cplxdb>>> ztensor_map;
    std::map<double, dtensor_map> tau_Wc_libri;
    std::map<double, ztensor_map> tau_Wc_libri_cplx;

    // NOTE: for case of a few atoms, some process may have much more memory load than others
    profiler.start("g0w0_build_spacetime_get_ap", "Compute balance atom-pairs distribution");
    const auto map_atpairs_balanced = get_balanced_ap_distribution_for_consec_descriptor(
        atbasis_abf, atbasis_abf, ad_Wc);
    const auto it_ap_myid = map_atpairs_balanced.find(ad_Wc.myid());
    if (it_ap_myid != map_atpairs_balanced.cend())
        ofs_myid << it_ap_myid->second << std::endl;
    else
        ofs_myid << "No atom pairs for Wc on this process" << std::endl;
    profiler.stop("g0w0_build_spacetime_get_ap");

    int wc_major_mask = 0;
    for (const auto &[tau, map_R_mat]: Wc_tau_R_blacs)
    {
        for (const auto &[R, mat_blacs]: map_R_mat)
        {
            if (mat_blacs.major() == MAJOR::ROW)
                wc_major_mask |= 1;
            else if (mat_blacs.major() == MAJOR::COL)
                wc_major_mask |= 2;
            else
                wc_major_mask |= 4;
        }
    }
    MPI_Allreduce(MPI_IN_PLACE, &wc_major_mask, 1, mpi_datatype<int>::value, MPI_BOR, ad_Wc.comm());

    MAJOR wc_major = MAJOR::AUTO;
    if (wc_major_mask == 1)
        wc_major = MAJOR::ROW;
    else if (wc_major_mask == 2)
        wc_major = MAJOR::COL;
    else if (wc_major_mask == 0)
    {
        throw LIBRPA_RUNTIME_ERROR("Wc(R,t) is empty");
    }
    else
    {
        throw LIBRPA_RUNTIME_ERROR("Inconsistent storage order among Wc(R,t) blocks");
    }

    ofs_myid << "wc_major == MAJOR::ROW ? " << std::boolalpha << (wc_major == MAJOR::ROW) << std::endl;

    IndexScheduler sched;
    // The scheduler indexes the existing BLACS buffers, so this must match Wc_tau_R_blacs.
    // LibRI's row-major requirement is handled below when each atom-pair block is wrapped.
    sched.init(map_atpairs_balanced, atbasis_abf, atbasis_abf, ad_Wc, wc_major == MAJOR::ROW);
    for (auto &[tau, map_R_mat]: Wc_tau_R_blacs)
    {
        for (auto &[R, mat_blacs]: map_R_mat)
        {
            auto pair_mat = get_ap_map_from_blacs_dist_scheduler(mat_blacs, sched, atbasis_abf, atbasis_abf, ad_Wc);
            auto itau = tfg.get_time_index(tau);
            // if (itau == 0)
            // {
            //     librpa_int::global::ofs_myid << "R: " << R << " tau: " << tau << std::endl;
            //     librpa_int::global::ofs_myid << "BLACS: " << std::endl << mat_blacs << std::endl;
            //     librpa_int::global::ofs_myid << "pair: " << std::endl << pair_mat << std::endl;
            //     std::stringstream ss;
            //     ss << "Wc_tau_R"
            //         << "_itau_" << std::setfill('0') << std::setw(5) << itau
            //         << "_iR_" << std::setfill('0') << std::setw(5) << get_R_index(pbc.Rlist, R) << ".csc";
            //     write_matrix_elsi_csc_parallel(ss.str(), mat_blacs, ad_Wc);
            // }
            // if (Params::debug)
            // {
            //     std::stringstream ss;
            //     ss << Params::output_dir << "Wc_tau_R"
            //         << "_itau_" << std::setfill('0') << std::setw(5) << tfg.get_time_index(tau)
            //         << "_iR_" << std::setfill('0') << std::setw(5) << get_R_index(Rlist, R) << ".csc";
            //     utils::write_matrix_elsi_csc_parallel(ss.str(), mat_blacs, envs::array_desc_abf_global);
            // }
            profiler.start("g0w0_build_spacetime_prep_Wc_all", "Prepare LibRI Wc object");
            for (auto &[pair, mat_ap]: pair_mat)
            {
                const auto &I = as_int(pair.first);
                const auto &J = as_int(pair.second);
                const auto &nabf_I = atbasis_abf.get_atom_nb(I);
                const auto &nabf_J = atbasis_abf.get_atom_nb(J);
                // LibRI Tensor is fixed to row major
                if (mat_ap.is_col_major()) mat_ap.swap_to_row_major();
                if (use_complex_tensor)
                    tau_Wc_libri_cplx[tau][I][{J, {R.x, R.y, R.z}}] = RI::Tensor<cplxdb>({nabf_I, nabf_J}, mat_ap.sptr());
                else
                    tau_Wc_libri[tau][I][{J, {R.x, R.y, R.z}}] = RI::Tensor<double>({nabf_I, nabf_J}, mat_ap.get_real().sptr());
                // std::stringstream ss;
                // ss << "tau_Wc_libri_"
                //     << "_itau_" << std::setfill('0') << std::setw(5) << itau
                //     << "_I_" << std::setfill('0') << std::setw(5) << I
                //     << "_J_" << std::setfill('0') << std::setw(5) << J
                //     << "_iR_" << std::setfill('0') << std::setw(5) << get_R_index(pbc.Rlist, R) << ".dat";
                // librpa_int::global::ofs_myid << "Writing to " << ss.str() << " " << tau << " " << I << " " << J << " " << R << std::endl;
                // std::ofstream ofs(ss.str());
                // // ofs << tau_Wc_libri[tau][I][{J, {R.x, R.y, R.z}}] << std::endl;
                // global::ofs_myid << "mat_ap row major? " << std::boolalpha << mat_ap.is_row_major() << std::endl;
                // ofs << mat_ap.get_real() << std::endl;
                // ofs.close();
            }
            // global::ofs_myid << "tau " << tau << endl << tau_Wc_libri[tau] << endl;
            mat_blacs.clear();
            profiler.stop("g0w0_build_spacetime_prep_Wc_all");
        }
    }
    profiler.stop("g0w0_build_spacetime_ct_ft_wc");

    lib_printf_root("Time for Fourier transform of Wc in GW (seconds, Wall/CPU): %f %f\n",
            profiler.get_wall_time_last("g0w0_build_spacetime_ct_ft_wc"),
            profiler.get_cpu_time_last("g0w0_build_spacetime_ct_ft_wc"));

    RI::GW<int, int, 3, double> gw_libri;
    RI::GW<int, int, 3, cplxdb> gw_libri_cplx;

    map<int,std::array<double,3>> atoms_pos;
    // Dummy atoms position
    for (int i = 0; i != natom; i++)
        atoms_pos.insert(pair<int, std::array<double, 3>>{i, {0, 0, 0}});

    global::profiler.start("g0w0_build_spacetime_2", "Setup LibRI G0W0 object and C data");
    if (use_complex_tensor)
    {
        gw_libri_cplx.set_parallel(comm_h.comm, atoms_pos, pbc.latvec_array, pbc.period_array);
        ztensor_map data_libri;
        for (const auto &I_JR_C : LRI_Cs.data_libri)
        {
            const auto I = I_JR_C.first;
            for (const auto &JR_C : I_JR_C.second)
            {
                const auto J = JR_C.first.first;
                const auto R = JR_C.first.second;
                const auto &C = JR_C.second;
                auto JR = std::pair<int, std::array<int, 3>>(J, R);
                data_libri[I][JR] = RI::Global_Func::convert<cplxdb>(C);
            }
        }
        gw_libri_cplx.set_Cs(data_libri, this->libri_threshold_C);
    }
    else
    {
        gw_libri.set_parallel(comm_h.comm, atoms_pos, pbc.latvec_array, pbc.period_array);
        gw_libri.set_Cs(LRI_Cs.data_libri, this->libri_threshold_C);
    }
    profiler.stop("g0w0_build_spacetime_2");
    lib_printf_root("Time for LibRI G0W0 setup (seconds, Wall/CPU): %f %f\n",
            profiler.get_wall_time_last("g0w0_build_spacetime_2"),
            profiler.get_cpu_time_last("g0w0_build_spacetime_2"));

    const auto tot_atpair_ordered = generate_atom_pair_from_nat(natom, true);
    const auto &Rlist = this->pbc.Rlist;
    auto IJR_local_gf = dispatch_vector_prod(tot_atpair_ordered, Rlist, ad_Wc.myid(), ad_Wc.nprocs(), true, false);
    ofs_myid << "IJR_local_gf: " << IJR_local_gf << std::endl;

    const int nfreq = tfg.get_n_grids();

    for (auto itau = 0; itau != nfreq; itau++)
    {
        // librpa_int::global::lib_printf("task %d itau %d start\n", mpi_comm_global_h.myid, itau);
        const auto tau = tfg.get_time_nodes()[itau];
        // librpa_int::global::lib_printf("task %d Wc_tau_R.count(tau) %zu\n", mpi_comm_global_h.myid, Wc_tau_R.count(tau));
        profiler.start("g0w0_build_spacetime_3", "Setup LibRI Wc");
        size_t n_obj_wc_libri = 0;
        if (use_complex_tensor)
        {
            auto it = tau_Wc_libri_cplx.find(tau);
            const auto &Wc_libri = (it == tau_Wc_libri_cplx.end())? ztensor_map{} : it->second;
            for (const auto &w: Wc_libri)
            {
                n_obj_wc_libri += w.second.size();
            }
            gw_libri_cplx.set_Ws(Wc_libri, this->libri_threshold_Wc);
            if (it != tau_Wc_libri_cplx.end()) tau_Wc_libri_cplx.erase(it);
        }
        else
        {
            auto it = tau_Wc_libri.find(tau);
            const auto &Wc_libri = (it == tau_Wc_libri.end())? dtensor_map{} : it->second;
            for (const auto &w: Wc_libri)
            {
                n_obj_wc_libri += w.second.size();
            }
            gw_libri.set_Ws(Wc_libri, this->libri_threshold_Wc);
            if (it != tau_Wc_libri.end()) tau_Wc_libri.erase(it);
        }
        profiler.stop("g0w0_build_spacetime_3");

        const int n_spinor = mf.get_n_spinor();

        for (int ispin = 0; ispin != mf.get_n_spins(); ispin++)
        {
            for (int ispinor_bra = 0; ispinor_bra < n_spinor; ispinor_bra++)
            {
                for (int ispinor_ket = 0; ispinor_ket < n_spinor; ispinor_ket++)
                {
                    dtensor_map sigc_posi_tau, sigc_nega_tau;
                    ztensor_map sigc_posi_tau_cplx, sigc_nega_tau_cplx;

                    profiler.start("g0w0_build_spacetime_4", "Compute G(R,t) and G(R,-t)");
                    // global::ofs_myid << "gf size " << gf.size() << endl;
                    // global::ofs_myid << "t " << gf[tau].size() << " ; -t " << gf[-tau].size() << endl;
                    std::map<double, dtensor_map> tau_gf_libri;
                    std::map<double, ztensor_map> tau_gf_libri_cplx;
                    const std::vector<double> taus{tau, -tau};
                    if (use_complex_tensor)
                    {
                        if (this->is_mf_eigvec_k_distributed_)
                        {
                            build_gf_libri_kpara(mf, comm_h, atbasis_wfc, ispin, ispinor_bra, ispinor_ket, this->pbc.kfrac_list,
                                                taus, IJR_local_gf, tau_gf_libri_cplx);
                        }
                        else
                        {
                            build_gf_libri_kserial(mf, atbasis_wfc, ispin, ispinor_bra, ispinor_ket, this->pbc.kfrac_list,
                                                taus, IJR_local_gf, tau_gf_libri_cplx);
                        }
                    }
                    else
                    {
                        if (this->is_mf_eigvec_k_distributed_)
                            build_gf_libri_kpara(mf, comm_h, atbasis_wfc, ispin, ispinor_bra, ispinor_ket, this->pbc.kfrac_list,
                                                taus, IJR_local_gf, tau_gf_libri);
                        else
                            build_gf_libri_kserial(mf, atbasis_wfc, ispin, ispinor_bra, ispinor_ket, this->pbc.kfrac_list,
                                                taus, IJR_local_gf, tau_gf_libri);
                    }
                    // if (itau == 0 && ispin == 0)  // debug
                    // {
                    //     global::ofs_myid << "tau_gf_libri itau=0 ispin=0" << std::endl;
                    //     global::ofs_myid << tau_gf_libri << std::endl;
                    // }
                    profiler.stop("g0w0_build_spacetime_4");

                    if (n_spinor > 1)
                    {
                        lib_printf_root("Time for Green's function, i_spin %d i_spinor_bra %d i_spinor_ket %d i_tau %d (seconds, Wall/CPU): %f %f\n",
                                ispin + 1, ispinor_bra + 1, ispinor_ket + 1, itau + 1,
                                profiler.get_wall_time_last("g0w0_build_spacetime_4"),
                                profiler.get_cpu_time_last("g0w0_build_spacetime_4"));
                    }
                    else
                    {
                        lib_printf_root("Time for Green's function, i_spin %d i_tau %d (seconds, Wall/CPU): %f %f\n",
                                ispin + 1, itau + 1,
                                profiler.get_wall_time_last("g0w0_build_spacetime_4"),
                                profiler.get_cpu_time_last("g0w0_build_spacetime_4"));
                    }

                    for (auto t: taus)
                    {
                        size_t n_obj_gf_libri = 0;
                        double wtime_g0w0_cal_sigc = omp_get_wtime();
                        if (use_complex_tensor)
                        {
                            const auto &gf_libri = tau_gf_libri_cplx.at(t);
                            for (const auto &gf: gf_libri)
                                n_obj_gf_libri += gf.second.size();

                            gw_libri_cplx.set_Gs(gf_libri, this->libri_threshold_G);
                            global::profiler.start("g0w0_build_spacetime_5", "Call libRI cal_Sigc");
                            gw_libri_cplx.cal_Sigmas();
                            global::profiler.stop("g0w0_build_spacetime_5");
                            global::profiler.start("g0w0_build_spacetime_5_clean");
                            gw_libri_cplx.free_Gs();
                            global::profiler.stop("g0w0_build_spacetime_5_clean");

                            // Check size of data
                            double mem_mb = get_tensor_map_bytes(gw_libri_cplx.Sigmas) * 1e-6;
                            global::ofs_myid << "Temporary Sigc_tau size for time " << t << " [MB]: " << mem_mb << std::endl;

                            if (t > 0)
                                sigc_posi_tau_cplx = std::move(gw_libri_cplx.Sigmas);
                            else
                                sigc_nega_tau_cplx = std::move(gw_libri_cplx.Sigmas);
                            gw_libri_cplx.Sigmas.clear();
                        }
                        else
                        {
                            // ofs_myid << tau_gf_libri << std::endl;
                            const auto &gf_libri = tau_gf_libri.at(t);
                            for (const auto &gf: gf_libri)
                                n_obj_gf_libri += gf.second.size();
                            // if (t > 0)  // debug
                            // {
                            //     global::ofs_myid << "gf_libri posi ispin=0 itau " << itau << " " << get_num_keys(gf_libri)  << std::endl;
                            //     print_keys(global::ofs_myid, gf_libri);
                            //     // global::ofs_myid << gf_libri << std::endl;
                            //     for (const auto &[I, JR_gf]: gf_libri)
                            //     {
                            //         for (const auto &[JR, gf]: JR_gf)
                            //         {
                            //             std::stringstream ss;
                            //             Vector3_Order<int> R(JR.second[0], JR.second[1], JR.second[2]);
                            //             ss << "gf_libri_posi"
                            //                 << "_itau_" << std::setfill('0') << std::setw(5) << itau
                            //                 << "_I_" << std::setfill('0') << std::setw(5) << I
                            //                 << "_J_" << std::setfill('0') << std::setw(5) << JR.first
                            //                 << "_iR_" << std::setfill('0') << std::setw(5) << get_R_index(pbc.Rlist, R) << ".dat";
                            //             librpa_int::global::ofs_myid << "Writing GF to " << ss.str() << " " << tau << " " << I << " " << JR.second << " " << R << std::endl;
                            //             std::ofstream ofs(ss.str());
                            //             ofs << gf << std::endl;
                            //             ofs.close();
                            //         }
                            //     }
                            // }
                            // else
                            // {
                            //     global::ofs_myid << "gf_libri nega ispin=0 itau " << itau << " " << get_num_keys(gf_libri) << std::endl;
                            //     print_keys(global::ofs_myid, gf_libri);
                            //     // global::ofs_myid << gf_libri << std::endl;
                            //     for (const auto &[I, JR_gf]: gf_libri)
                            //     {
                            //         for (const auto &[JR, gf]: JR_gf)
                            //         {
                            //             std::stringstream ss;
                            //             Vector3_Order<int> R(JR.second[0], JR.second[1], JR.second[2]);
                            //             ss << "gf_libri_nega"
                            //                 << "_itau_" << std::setfill('0') << std::setw(5) << itau
                            //                 << "_I_" << std::setfill('0') << std::setw(5) << I
                            //                 << "_J_" << std::setfill('0') << std::setw(5) << JR.first
                            //                 << "_iR_" << std::setfill('0') << std::setw(5) << get_R_index(pbc.Rlist, R) << ".dat";
                            //             librpa_int::global::ofs_myid << "Writing GF to " << ss.str() << " " << tau << " " << I << " " << JR.second << " " << R << std::endl;
                            //             std::ofstream ofs(ss.str());
                            //             ofs << gf << std::endl;
                            //             ofs.close();
                            //         }
                            //     }
                            // }
                            gw_libri.set_Gs(gf_libri, this->libri_threshold_G);
                            global::profiler.start("g0w0_build_spacetime_5", "Call libRI cal_Sigc");
                            gw_libri.cal_Sigmas();
                            global::profiler.stop("g0w0_build_spacetime_5");
                            global::profiler.start("g0w0_build_spacetime_5_clean");
                            gw_libri.free_Gs();
                            global::profiler.stop("g0w0_build_spacetime_5_clean");

                            // Check size of data
                            double mem_mb = get_tensor_map_bytes(gw_libri.Sigmas) * 1e-6;
                            global::ofs_myid << "Temporary Sigc_tau size for time " << t << " [MB]: " << mem_mb << std::endl;

                            if (t > 0)
                                sigc_posi_tau = std::move(gw_libri.Sigmas);
                            else
                                sigc_nega_tau = std::move(gw_libri.Sigmas);
                            gw_libri.Sigmas.clear();
                            // if (itau == 0 && ispin == 0)
                            // {
                            //     if (t > 0)
                            //     {
                            //         global::ofs_myid << "sigc_posi_tau itau=0 ispin=0 " << get_num_keys(sigc_posi_tau) << std::endl;
                            //         print_keys(global::ofs_myid, sigc_posi_tau);
                            //         global::ofs_myid << sigc_posi_tau << std::endl;
                            //     }
                            //     if (t < 0)
                            //     {
                            //         global::ofs_myid << "sigc_nega_tau itau=0 ispin=0 " << get_num_keys(sigc_nega_tau) << std::endl;
                            //         print_keys(global::ofs_myid, sigc_nega_tau);
                            //         global::ofs_myid << sigc_nega_tau << std::endl;
                            //     }
                            // }
                        }
                        wtime_g0w0_cal_sigc = omp_get_wtime() - wtime_g0w0_cal_sigc;
                        if (n_spinor > 1)
                            librpa_int::global::lib_printf(
                                "Task %4d. libRI G0W0, spin %1d, bra %1d, ket %1d, time grid %12.6f. Wc size %zu, GF "
                                "size %zu. Wall time %f\n",
                                comm_h.myid, ispin, ispinor_bra, ispinor_ket, t, n_obj_wc_libri, n_obj_gf_libri,
                                wtime_g0w0_cal_sigc);
                        else
                            librpa_int::global::lib_printf(
                                "Task %4d. libRI G0W0, spin %1d, time grid %12.6f. Wc size %zu, GF "
                                "size %zu. Wall time %f\n",
                                comm_h.myid, ispin, t, n_obj_wc_libri, n_obj_gf_libri,
                                wtime_g0w0_cal_sigc);
                    }

                    size_t n_IJR_myid = 0; // for sigcmat output

                    if (this->output_sigc_mat_rt)
                    {
                        std::stringstream ss;
                        ss << path_as_directory(this->output_dir) << "SigcRT"
                            << "_ispin_" << std::setfill('0') << std::setw(5) << ispin;
                        if (n_spinor > 1)
                        {
                           ss << "_spinor_" << ispinor_bra << "_" << ispinor_ket;
                        }
                        ss << "_itau_" << std::setfill('0') << std::setw(5) << itau
                           << "_myid_" << std::setfill('0') << std::setw(5) << global::myid_global << ".dat";
                        ofs_sigmac_r.open(ss.str(), std::ios::out | std::ios::binary);
                        ofs_sigmac_r.write((char *) &n_IJR_myid, sizeof(size_t)); // placeholder
                    }

                    auto accumulate_sigc_tau_to_freq = [&](const auto &sigc_posi_tau_in,
                                                           const auto &sigc_nega_tau_in,
                                                           auto block_scale_factor,
                                                           auto make_freq_value, size_t elem_size)
                    {
                        for (const auto &[I, J_R_sigc_posi] : sigc_posi_tau_in)
                        {
                            const auto n_I = this->atbasis_wfc.get_atom_nb(I);

                            auto it_nega_I = sigc_nega_tau_in.find(I);
                            if (it_nega_I == sigc_nega_tau_in.cend()) continue;

                            for (const auto &[JR, sigc_posi_block] : J_R_sigc_posi)
                            {
                                const auto J = JR.first;
                                const auto n_J = this->atbasis_wfc.get_atom_nb(J);
                                const auto &Ra = JR.second;

                                auto it_nega_JR = it_nega_I->second.find(JR);
                                if (it_nega_JR == it_nega_I->second.cend()) continue;

                                const auto &sigc_nega_block = it_nega_JR->second;

                                const Vector3_Order<int> R{Ra[0], Ra[1], Ra[2]};

                                const auto it_R = std::find(Rlist.cbegin(), Rlist.cend(), R);
                                const auto iR = std::distance(Rlist.cbegin(), it_R);

                                const auto sigc_cos = block_scale_factor * (sigc_posi_block + sigc_nega_block);
                                const auto sigc_sin = block_scale_factor * (sigc_posi_block - sigc_nega_block);

                                ++n_IJR_myid;

                                if (this->output_sigc_mat_rt)
                                {
                                    size_t dims[5];
                                    dims[0] = as_size(iR);
                                    dims[1] = as_size(I);
                                    dims[2] = as_size(J);
                                    dims[3] = as_size(n_I);
                                    dims[4] = as_size(n_J);

                                    ofs_sigmac_r.write(reinterpret_cast<const char *>(dims),
                                                       5 * sizeof(size_t));

                                    ofs_sigmac_r.write(
                                        reinterpret_cast<const char *>(sigc_cos.ptr()),
                                        n_I * n_J * elem_size);

                                    ofs_sigmac_r.write(
                                        reinterpret_cast<const char *>(sigc_sin.ptr()),
                                        n_I * n_J * elem_size);
                                }

                                for (size_t iomega = 0; iomega != tfg.get_n_grids(); ++iomega)
                                {
                                    const auto omega = tfg.get_freq_nodes()[iomega];
                                    const auto t2f_sin = tfg.get_sintrans_t2f()(iomega, itau);
                                    const auto t2f_cos = tfg.get_costrans_t2f()(iomega, itau);

                                    Matz sigc_temp(n_I, n_J, MAJOR::ROW);

                                    for (size_t i = 0; i != static_cast<size_t>(n_I); ++i)
                                    {
                                        for (size_t j = 0; j != static_cast<size_t>(n_J); ++j)
                                        {
                                            sigc_temp(i, j) = make_freq_value(
                                                sigc_cos(i, j), sigc_sin(i, j), t2f_cos, t2f_sin);
                                        }
                                    }

                                    const atpair_t IJ{I, J};
                                    auto &m_IJ =
                                        sigc_is_f_IJ_R[ispin][ispinor_bra][ispinor_ket][omega][IJ];

                                    auto it = m_IJ.find(R);
                                    if (it == m_IJ.cend())
                                    {
                                        m_IJ.emplace(R, std::move(sigc_temp));
                                    }
                                    else
                                    {
                                        it->second += sigc_temp;
                                    }
                                }
                            }
                        }
                    };

                    // symmetrize and perform transformation
                    global::profiler.start("g0w0_build_spacetime_6", "Transform Sigc (R,t) -> (R,w)");
                    if (use_complex_tensor)
                    {
                        accumulate_sigc_tau_to_freq(
                            sigc_posi_tau_cplx, sigc_nega_tau_cplx,
                            cplxdb(0.5, 0.0),
                            [](const cplxdb &sigc_cos, const cplxdb &sigc_sin, double t2f_cos,
                               double t2f_sin) -> cplxdb
                            {
                                return sigc_cos * t2f_cos + sigc_sin * cplxdb(0.0, t2f_sin);
                            },
                            sizeof(cplxdb));
                        sigc_posi_tau_cplx.clear();
                        sigc_nega_tau_cplx.clear();
                    }
                    else
                    {
                        accumulate_sigc_tau_to_freq(
                            sigc_posi_tau, sigc_nega_tau,
                            0.5,
                            [](double sigc_cos, double sigc_sin, double t2f_cos,
                               double t2f_sin) -> cplxdb
                            { return cplxdb{sigc_cos * t2f_cos, sigc_sin * t2f_sin}; },
                            sizeof(double));
                        sigc_posi_tau.clear();
                        sigc_nega_tau.clear();
                    }

                    global::profiler.stop("g0w0_build_spacetime_6");

                    if (this->output_sigc_mat_rt)
                    {
                        ofs_sigmac_r.seekp(0);
                        ofs_sigmac_r.write((char *) &n_IJR_myid, sizeof(size_t)); // overwrite
                        ofs_sigmac_r.close();
                    }
                }
            }
        }
        global::profiler.start("g0w0_build_spacetime_free_Ws");
        gw_libri.free_Ws();
        global::profiler.stop("g0w0_build_spacetime_free_Ws");
        // Release freed memory to OS, to resolve memory fragments in LibRI
        release_free_mem();
    }
    is_rspace_built_ = true;
#endif

    // Export real-space imaginary-frequency NAO sigma_c matrices
    if (is_rspace_built_ && this->output_sigc_mat_rf)
    {
        const int n_spinor = mf.get_n_spinor();
        for (int ispin = 0; ispin != mf.get_n_spins(); ispin++)
        {
            for (int ispinor_bra = 0; ispinor_bra < n_spinor; ispinor_bra++)
            {
                for (int ispinor_ket = 0; ispinor_ket < n_spinor; ispinor_ket++)
                {
                    for (size_t iomega = 0; iomega != tfg.get_n_grids(); iomega++)
                    {
                        size_t n_IJR_myid = 0;
                        std::stringstream ss;
                        ss << path_as_directory(this->output_dir) << "SigcRF"
                           << "_ispin_" << std::setfill('0') << std::setw(5) << ispin;
                        if (n_spinor > 1)
                        {
                           ss << "_spinor_" << ispinor_bra << "_" << ispinor_ket;
                        }
                        ss << "_iomega_" << std::setfill('0') << std::setw(5) << iomega
                           << "_myid_" << std::setfill('0') << std::setw(5) << global::myid_global << ".dat";
                        ofs_sigmac_r.open(ss.str(), std::ios::out | std::ios::binary);
                        ofs_sigmac_r.write((char *) &n_IJR_myid, sizeof(size_t)); // placeholder

                        const auto omega = tfg.get_freq_nodes()[iomega];
                        const auto &sigc_IJ_R = sigc_is_f_IJ_R[ispin][ispinor_bra][ispinor_ket][omega];
                        for (const auto &[IJ, R_sigc]: sigc_IJ_R)
                        {
                            const int I = IJ.first;
                            const int J = IJ.second;
                            for (const auto &[R, sigc]: R_sigc)
                            {
                                const auto iR = this->pbc.get_R_index(R);
                                const auto &n_I = this->atbasis_wfc.get_atom_nb(I);
                                const auto &n_J = this->atbasis_wfc.get_atom_nb(J);
                                n_IJR_myid++;
                                size_t dims[5];
                                dims[0] = iR;
                                dims[1] = I;
                                dims[2] = J;
                                dims[3] = n_I;
                                dims[4] = n_J;
                                assert (sigc.size() == n_I * n_J);
                                ofs_sigmac_r.write((char *) dims, 5 * sizeof(size_t));
                                ofs_sigmac_r.write((char *) sigc.ptr(), sigc.size() * sizeof(cplxdb));
                            }
                        }
                        ofs_sigmac_r.seekp(0);
                        ofs_sigmac_r.write((char *) &n_IJR_myid, sizeof(size_t)); // placeholder
                        ofs_sigmac_r.close();
                    }
                }
            }
        }
    }
}

void G0W0::build_sigc_matrix_KS(const std::map<int, std::map<int, std::map<int, ComplexMatrix>>> &wfc_target,
                                const std::vector<Vector3_Order<double>> &kfrac_target,
                                const Atoms &geometry)
{
    throw LIBRPA_RUNTIME_ERROR("Not implemented yet");
}

void G0W0::build_sigc_matrix_KS_blacs(const std::map<int, std::map<int, std::map<int, ComplexMatrix>>> &wfc_target,
                                      const std::vector<Vector3_Order<double>> &kfrac_target,
                                      const Atoms &geometry, const BlacsCtxtHandler &blacs_ctxt_h)
{
    assert(blacs_ctxt_h.comm() == this->comm_h.comm);
    assert(this->is_rspace_built_);

    if (this->is_kspace_built_)
    {
        global::lib_printf("Warning: reset Sigmac_c k-space matrices\n");
        this->reset_kspace();
    }

    const int n_aos = mf.get_n_aos();
    const int n_bands = mf.get_n_bands();
    const int n_spins = mf.get_n_spins();
    const int n_spinor = mf.get_n_spinor();
    const auto &period = pbc.period;
    const auto &coord_frac = geometry.coords_frac;

    global::profiler.start("g0w0_build_sigc_KS");

#ifndef LIBRPA_USE_LIBRI
    if (global::mpi_comm_global_h.myid == 0)
    {
        std::cout << "LIBRA::G0W0::build_sigc_matrix_KS is only implemented on top of LibRI" << std::endl;
        std::cout << "Please recompile LibRPA with -DLIBRPA_USE_LIBRI and optionally configure include path" << std::endl;
    }
    global::mpi_comm_global_h.barrier();
    throw LIBRPA_RUNTIME_ERROR("G0W0 needs compilation with LibRI");
#else
    ArrayDesc desc_nband_nao(blacs_ctxt_h);
    desc_nband_nao.init_1b1p(n_bands, n_aos, 0, 0);
    ArrayDesc desc_nao_nao(blacs_ctxt_h);
    desc_nao_nao.init_1b1p(n_aos, n_aos, 0, 0);
    ArrayDesc desc_nband_nband(blacs_ctxt_h);
    desc_nband_nband.init_1b1p(n_bands, n_bands, 0, 0);

    ArrayDesc desc_nao_nband_fb(blacs_ctxt_h);  // For k-parallel
    ArrayDesc desc_nband_nband_fb(blacs_ctxt_h);

    const auto set_IJ_nao_nao = get_necessary_IJ_from_block_2D(
        this->atbasis_wfc, this->atbasis_wfc, desc_nao_nao);
    const auto s0_s1 = get_s0_s1_for_comm_map2_first(set_IJ_nao_nao);

    // Check Sigmac matrix distribution
    if (!is_rspace_redist_for_KS_)
    {
        global::profiler.start("g0w0_build_sigc_KS_rspace_redist");
        for (int isp = 0; isp < n_spins; isp++)
        {
            for (int ispn_bra = 0; ispn_bra < n_spinor; ispn_bra++)
            {
                for (int ispn_ket = 0; ispn_ket < n_spinor; ispn_ket++)
                {
                    auto sigc_orig = find_nested_int_map_3(sigc_is_f_IJ_R, isp, ispn_bra, ispn_ket);
                    for (const auto& freq: this->tfg.get_freq_nodes())
                    {
                        std::map<int, std::map<std::pair<int, std::array<int, 3>>, Tensor<cplxdb>>> sigc_I_JR_local;
                        if (sigc_orig != nullptr)
                        {
                            auto it_sp_f = sigc_orig->find(freq);
                            if (it_sp_f != sigc_orig->cend())
                            {
                                for (const auto &IJ_R_sigc: it_sp_f->second)
                                {
                                    const auto IJ = IJ_R_sigc.first;
                                    const int I = IJ.first;
                                    const int J = IJ.second;
                                    const auto &n_I = this->atbasis_wfc.get_atom_nb(I);
                                    const auto &n_J = this->atbasis_wfc.get_atom_nb(J);
                                    for (const auto &[R, sigc]: IJ_R_sigc.second)
                                    {
                                        const std::array<int, 3> Ra{R.x, R.y, R.z};
                                        sigc_I_JR_local[I][{J, Ra}] = Tensor<complex<double>>({n_I, n_J}, sigc.sptr());
                                    }
                                }
                            }
                        }
                        // for certain spin and frequency
                        auto sigc_I_JR = comm_map2_first(blacs_ctxt_h.comm(), sigc_I_JR_local, s0_s1.first, s0_s1.second);
                        sigc_I_JR_local.clear();
                        if (sigc_I_JR.size() == 0) continue;
                        ap_p_map<map<Vector3_Order<int>, Matz>> sigc_new;
                        for (const auto &[I, JRmat]: sigc_I_JR)
                        {
                            const int n_I = this->atbasis_wfc.get_atom_nb(I);
                            for (const auto &[JR, mat]: JRmat)
                            {
                                const atom_t J = JR.first;
                                atpair_t IJ{I, J};
                                const int n_J = this->atbasis_wfc.get_atom_nb(J);
                                const auto &Ra = JR.second;
                                const Vector3_Order<int> R{Ra[0], Ra[1], Ra[2]};
                                sigc_new[IJ][R] = Matz{n_I, n_J, mat.data, MAJOR::ROW};
                            }
                        }
                        if (sigc_orig == nullptr)
                            sigc_is_f_IJ_R[isp][ispn_bra][ispn_ket][freq] = sigc_new;
                        else
                        {
                            auto it_sp_f = sigc_orig->find(freq);
                            if (it_sp_f == sigc_orig->cend())
                                sigc_orig->emplace(freq, sigc_new);
                            else
                                it_sp_f->second.swap(sigc_new);
                        }
                        sigc_I_JR.clear();
                    }
                }
            }
        }

        is_rspace_redist_for_KS_ = true;
        is_rspace_redist_blacs_ = true;
        global::profiler.stop("g0w0_build_sigc_KS_rspace_redist");
    }
    else
    {
        if (!is_rspace_redist_blacs_)
            throw LIBRPA_RUNTIME_ERROR("has redistributed for non-BLACS");
    }

    sigc_is_ik_f_KS.clear();

    // local 2D-block submatrices
    auto sigc_nao_nao = init_local_mat<complex<double>>(desc_nao_nao, MAJOR::COL);
    auto sigc_nband_nband = init_local_mat<complex<double>>(desc_nband_nband, MAJOR::COL);

    for (int isp = 0; isp < n_spins; isp++)
    {
	// Initialize, make sure the map on every access has the isp key
        this->sigc_is_ik_f_KS[isp] = {};

        for (int ispn_bra = 0; ispn_bra < n_spinor; ispn_bra++)
        {
            for (int ispn_ket = 0; ispn_ket < n_spinor; ispn_ket++)
            {
                map<double, map<atom_t, map<atom_t, map<Vector3_Order<int>, Matz>>>> sigc_isp_local;
                global::profiler.start("g0w0_build_sigc_KS_find_bvk");
                for (const auto& freq: this->tfg.get_freq_nodes())
                {
                    const auto &sigc_IJ_R = sigc_is_f_IJ_R.at(isp).at(ispn_bra).at(ispn_ket).at(freq);
                    sigc_isp_local[freq] = {};
                    // Convert each <I,<J, R>> pair to the nearest neighbour to speed up later Fourier transform
                    // while keep the accuracy in further band interpolation.
                    // Reuse the cleared-up sigc_I_JR_local object
                    if (geometry.is_frac_set())
                    {
                        // NOTE: from redistribution, sigc_is_f_IJ_R is ensured to have all [spin][freq] keys,
                        // but each value (atom-pair map) can be empty.
                        for (const auto &[IJ, R_sigc]: sigc_IJ_R)
                        {
                            const auto &I = IJ.first;
                            const auto &J = IJ.second;
                            for (auto &[R, sigc]: R_sigc)
                            {
                                auto distsq = std::numeric_limits<double>::max();
                                Vector3<int> R_IJ;
                                Vector3_Order<int> R_bvk;
                                for (int i = -1; i < 2; i++)
                                {
                                    R_IJ.x = i * period.x + R.x;
                                    for (int j = -1; j < 2; j++)
                                    {
                                        R_IJ.y = j * period.y + R.y;
                                        for (int k = -1; k < 2; k++)
                                        {
                                            R_IJ.z = k * period.z + R.z;
                                            const auto diff =
                                                (Vector3<double>(coord_frac.at(I)[0], coord_frac.at(I)[1],
                                                                coord_frac.at(I)[2]) -
                                                Vector3<double>(coord_frac.at(J)[0], coord_frac.at(J)[1],
                                                                coord_frac.at(J)[2]) -
                                                Vector3<double>(R_IJ.x, R_IJ.y, R_IJ.z)) * pbc.latvec;
                                            const auto norm2 = diff.norm2();
                                            if (norm2 < distsq)
                                            {
                                                distsq = norm2;
                                                R_bvk = R_IJ;
                                            }
                                        }
                                    }
                                }
                                sigc_isp_local[freq][I][J][R_bvk] = sigc;
                            }
                        }
                    }
                    else
                    {
                        for (const auto &[IJ, R_sigc]: sigc_IJ_R)
                        {
                            const auto &I = IJ.first;
                            const auto &J = IJ.second;
                            sigc_isp_local[freq][I][J] = R_sigc;
                        }
                    }
                }
                global::profiler.stop("g0w0_build_sigc_KS_find_bvk");

                if (is_mf_eigvec_k_distributed_)
                {
                    global::profiler.start("g0w0_build_sigc_KS_rotate_kpara");
                    auto temp_nband_nao = init_local_mat<complex<double>>(desc_nband_nao, MAJOR::COL);
                    // collect parsed eigenvectors all all processes
                    std::vector<int> nks_all(comm_h.nprocs, 0);
                    std::vector<int> iks_local;
                    auto wfc_sp = find_nested_int_map_2(wfc_target, isp, ispn_bra);
                    if (wfc_sp != nullptr)
                    {
                        nks_all[comm_h.myid] = wfc_sp->size();
                        for (const auto &[ik, _]: *wfc_sp)
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
                        auto sigc_nband_nband_fb = init_local_mat<complex<double>>(desc_nband_nband_fb, MAJOR::COL);
                        for (int ik_this = 0; ik_this < nk_this; ik_this++)
                        {
                            const int ik = iks_all[pid * nk_max + ik_this];
                            const auto& kfrac = kfrac_target[ik];
                            const std::function<complex<double>(const atom_t, const atom_t, const Vector3_Order<int> &)>
                                fourier = [kfrac](const atom_t I, const atom_t J, const Vector3_Order<int> &R)
                                {
                                    const auto ang = (kfrac * R) * TWO_PI;
                                    return complex<double>{std::cos(ang), std::sin(ang)};
                                };
                            std::vector<complex<double>> dummy(1);
                            for (const auto& freq: this->tfg.get_freq_nodes())
                            {
                                sigc_nao_nao.zero_out();
                                sigc_nband_nband_fb.zero_out();
                                collect_block_from_IJ_storage_matrix_transform(sigc_nao_nao, desc_nao_nao,
                                        this->atbasis_wfc, this->atbasis_wfc, fourier, sigc_isp_local.at(freq));
                                if (pid == comm_h.myid)
                                {
                                    const auto &wfc_bra = wfc_target.at(isp).at(ispn_bra).at(ik);
                                    const auto &wfc_ket = wfc_target.at(isp).at(ispn_ket).at(ik);
                                    // processing the k-point eigenvector on this process
                                    ScalapackConnector::pgemm_f('C', 'N', n_bands, n_aos, n_aos, 1.0,
                                                                wfc_bra.c, 1, 1, desc_nao_nband_fb.desc,
                                                                sigc_nao_nao.ptr(), 1, 1, desc_nao_nao.desc,
                                                                0.0,
                                                                temp_nband_nao.ptr(), 1, 1, desc_nband_nao.desc);
                                    ScalapackConnector::pgemm_f('N', 'N', n_bands, n_bands, n_aos, 1.0,
                                                                temp_nband_nao.ptr(), 1, 1, desc_nband_nao.desc,
                                                                wfc_ket.c, 1, 1, desc_nao_nband_fb.desc,
                                                                0.0,
                                                                sigc_nband_nband_fb.ptr(), 1, 1, desc_nband_nband_fb.desc);
                                    auto sigc_freq = find_nested_int_map_2(this->sigc_is_ik_f_KS, isp, ik);
                                    if (sigc_freq == nullptr || sigc_freq->count(freq) == 0)
                                        this->sigc_is_ik_f_KS[isp][ik][freq] = Matz(n_bands, n_bands, MAJOR::COL);
                                    this->sigc_is_ik_f_KS[isp][ik][freq] += sigc_nband_nband_fb;
                                }
                                else
                                {
                                    // processing the k-point eigenvector at other processes
                                    ScalapackConnector::pgemm_f('C', 'N', n_bands, n_aos, n_aos, 1.0,
                                                                dummy.data(), 1, 1, desc_nao_nband_fb.desc,
                                                                sigc_nao_nao.ptr(), 1, 1, desc_nao_nao.desc,
                                                                0.0,
                                                                temp_nband_nao.ptr(), 1, 1, desc_nband_nao.desc);
                                    ScalapackConnector::pgemm_f('N', 'N', n_bands, n_bands, n_aos, 1.0,
                                                                temp_nband_nao.ptr(), 1, 1, desc_nband_nao.desc,
                                                                dummy.data(), 1, 1, desc_nao_nband_fb.desc,
                                                                0.0,
                                                                sigc_nband_nband_fb.ptr(), 1, 1, desc_nband_nband_fb.desc);
                                }
                            }
                        }
                    }
                    global::profiler.stop("g0w0_build_sigc_KS_rotate_kpara");
                }
                else
                {
                    global::profiler.start("g0w0_build_sigc_KS_rotate_kserial");
                    desc_nband_nband_fb.init(n_bands, n_bands, n_bands, n_bands, 0, 0);
                    auto sigc_nband_nband_fb = init_local_mat<complex<double>>(desc_nband_nband_fb, MAJOR::COL);
                    for (const auto& freq: this->tfg.get_freq_nodes())
                    {
                        // Perform Fourier transform and rotate
                        for (size_t ik = 0; ik < kfrac_target.size(); ik++)
                        {
                            const auto kfrac = kfrac_target[ik];
                            // const std::function<complex<double>(const atom_t, const atom_t, const Vector3_Order<int> &)>
                            const std::function<complex<double>(const atom_t, const atom_t, const Vector3_Order<int> &)>
                                fourier = [kfrac](const atom_t I, const atom_t J, const Vector3_Order<int> &R)
                                {
                                    const auto ang = (kfrac * R) * TWO_PI;
                                    return complex<double>{std::cos(ang), std::sin(ang)};
                                };

                            sigc_nao_nao.zero_out();
                            sigc_nband_nband_fb.zero_out();
                            collect_block_from_IJ_storage_matrix_transform(sigc_nao_nao, desc_nao_nao,
                                    this->atbasis_wfc, this->atbasis_wfc, fourier, sigc_isp_local.at(freq));
                            // prepare wave function BLACS
                            // FIXME: check if bra and ket are correct
                            const auto &wfc_bra = wfc_target.at(isp).at(ispn_bra).at(ik);
                            const auto &wfc_ket = wfc_target.at(isp).at(ispn_ket).at(ik);
                            blacs_ctxt_h.barrier();
                            const auto wfc_block = get_local_mat(wfc_bra.c, MAJOR::ROW, desc_nband_nao, MAJOR::COL).conj();
                            auto temp_nband_nao = multiply_scalapack(wfc_block, desc_nband_nao, sigc_nao_nao, desc_nao_nao, desc_nband_nao);
                            if (ispn_bra == ispn_ket)
                                ScalapackConnector::pgemm_f('N', 'C', n_bands, n_bands, n_aos, 1.0,
                                                            temp_nband_nao.ptr(), 1, 1, desc_nband_nao.desc,
                                                            wfc_block.ptr(), 1, 1, desc_nband_nao.desc, 0.0,
                                                            sigc_nband_nband.ptr(), 1, 1, desc_nband_nband.desc);
                            else
                            {
                                const auto wfc_ket_block = get_local_mat(wfc_ket.c, MAJOR::ROW, desc_nband_nao, MAJOR::COL);
                                ScalapackConnector::pgemm_f('N', 'T', n_bands, n_bands, n_aos, 1.0,
                                                            temp_nband_nao.ptr(), 1, 1, desc_nband_nao.desc,
                                                            wfc_ket_block.ptr(), 1, 1, desc_nband_nao.desc, 0.0,
                                                            sigc_nband_nband.ptr(), 1, 1, desc_nband_nband.desc);
                            }
                            // collect the full matrix to master
                            // TODO: would need a different strategy for large system
                            ScalapackConnector::pgemr2d_f(n_bands, n_bands,
                                                          sigc_nband_nband.ptr(), 1, 1, desc_nband_nband.desc,
                                                          sigc_nband_nband_fb.ptr(), 1, 1, desc_nband_nband_fb.desc,
                                                          desc_nband_nband_fb.ictxt());
                            // NOTE: only the matrices at master process is meaningful
                            auto sigc_freq = find_nested_int_map_2(this->sigc_is_ik_f_KS, isp, ik);
                            if (sigc_freq == nullptr || sigc_freq->count(freq) == 0)
                            {
                                this->sigc_is_ik_f_KS[isp][ik][freq] = sigc_nband_nband_fb.copy();
                            }
                            else
                            {
                                sigc_freq->at(freq) += sigc_nband_nband_fb;
                            }
                        }
                    }
                    global::profiler.stop("g0w0_build_sigc_KS_rotate_kserial");
                }
            }
        }
    }
#endif
    this->is_kspace_built_ = true;
    global::profiler.stop("g0w0_build_sigc_KS");
}

void G0W0::build_sigc_matrix_KS_kgrid()
{
    comm_h.barrier();
    librpa_int::global::ofs_myid << "build_sigc_matrix_KS_kgrid: constructing self-energy matrix for SCF k-grid" << std::endl;
    this->build_sigc_matrix_KS(this->mf.get_eigenvectors(), this->pbc.kfrac_list, {});
    // if (comm_h.is_root())
    // {
    //     write_self_energy_omega("self_energy_omega.dat", *this);
    // }
}

void G0W0::build_sigc_matrix_KS_band(const std::map<int, std::map<int, std::map<int, ComplexMatrix>>> &wfc,
                                     const std::vector<Vector3_Order<double>> &kfrac_band,
                                     const Atoms &geometry)
{
    comm_h.barrier();
    if (comm_h.myid == 0)
    {
        librpa_int::global::lib_printf("build_sigc_matrix_KS_kgrid: constructing self-energy matrix for band k-path\n");
    }
    this->build_sigc_matrix_KS(wfc, kfrac_band, geometry);
}

void G0W0::build_sigc_matrix_KS_kgrid_blacs(const BlacsCtxtHandler &blacs_ctxt_h)
{
    comm_h.barrier();
    librpa_int::global::ofs_myid << "build_sigc_matrix_KS_kgrid: constructing self-energy matrix for SCF k-grid with BLACS" << std::endl;
    this->build_sigc_matrix_KS_blacs(this->mf.get_eigenvectors(), this->pbc.kfrac_list, {}, blacs_ctxt_h);
    // if (comm_h.is_root())
    // {
    //     write_self_energy_omega("self_energy_omega.dat", *this);
    // }
}

void G0W0::build_sigc_matrix_KS_band_blacs(
    const std::map<int, std::map<int, std::map<int, ComplexMatrix>>> &wfc,
    const std::vector<Vector3_Order<double>> &kfrac_band, const Atoms &geometry,
    const BlacsCtxtHandler &blacs_ctxt_h)
{
    comm_h.barrier();
    if (comm_h.myid == 0)
    {
        librpa_int::global::lib_printf("build_sigc_matrix_KS_band: constructing self-energy matrix for band k-path with BLACS\n");
    }
    this->build_sigc_matrix_KS_blacs(wfc, kfrac_band, geometry, blacs_ctxt_h);
}

} // namespace librpa_int
