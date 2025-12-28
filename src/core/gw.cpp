#include "gw.h"

// Public API headers
#include "librpa_enums.h"

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
static void build_gf_libri_kpara(
    const MeanField &mf,
    const MpiCommHandler &comm_h,
    const AtomicBasis &atbasis_wfc,
    int ispin,
    const vector<Vector3_Order<double>> &kfrac_list,
    const std::vector<double> &taus,
    const std::vector<std::pair<atpair_t, Vector3_Order<int>>> IJRs,
    std::map<double, std::map<int, std::map<std::pair<int,std::array<int,3>>,RI::Tensor<double>>>> &tau_gf_libri)
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
    auto gf_taus_Rs_cplx = get_gf_cplx_imagtimes_Rs_kpara(ispin, mf, kfrac_list, taus, Rs_this, comm_h);
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
                auto p_gfmat = std::make_shared<std::valarray<double>>(gf_block.real().c, gf_block.size);
                omp_set_lock(&gf_lock);
                tau_gf_libri[tau][I][{J, Ra}] = RI::Tensor<double>({as_size(gf_block.nr), as_size(gf_block.nc)}, p_gfmat);
                omp_unset_lock(&gf_lock);
            }
#pragma omp barrier
            omp_destroy_lock(&gf_lock);
        }
        // global::ofs_myid << R << std::endl;
        // print_complex_matrix("test", dmat_cplx, global::ofs_myid, true);
    }

    global::profiler.stop("g0w0_build_gf_libri_kpara");
}

static void build_gf_libri_kserial(
    const MeanField &mf,
    const AtomicBasis &atbasis_wfc,
    int ispin,
    const vector<Vector3_Order<double>> &kfrac_list,
    const std::vector<double> &taus,
    const std::vector<std::pair<atpair_t, Vector3_Order<int>>> IJRs,
    std::map<double, std::map<int, std::map<std::pair<int,std::array<int,3>>,RI::Tensor<double>>>> &tau_gf_libri)
{
    global::profiler.start("g0w0_build_gf_libri_kserial");
    std::set<Vector3_Order<int>> Rs_local;
    for (const auto &IJR: IJRs)
    {
        Rs_local.insert(IJR.second);
    }
    auto gf = mf.get_gf_real_imagtimes_Rs(ispin, kfrac_list, taus, {Rs_local.cbegin(), Rs_local.cend()});
    // global::ofs_myid << "gf " << gf << std::endl;
    tau_gf_libri.clear();
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
            matrix gf_IJ_block(n_I, n_J);
            {
                for (size_t i = 0; i != n_I; i++)
                    for (size_t j = 0; j != n_J; j++)
                    {
                        gf_IJ_block(i, j) = gf_global(atbasis_wfc.get_global_index(I, i),
                                                        atbasis_wfc.get_global_index(J, j));
                    }
                std::shared_ptr<std::valarray<double>> mat_ptr = std::make_shared<std::valarray<double>>(gf_IJ_block.c, gf_IJ_block.size);
                tau_gf_libri[t][as_int(I)][{as_int(J), {R.x, R.y, R.z}}] = RI::Tensor<double>({n_I, n_J}, mat_ptr);
            }
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
    using librpa_int::global::mpi_comm_global_h;

    if(parallel_routing != LibrpaParallelRouting::LIBRI)
    {
        mpi_comm_global_h.barrier();
        throw LIBRPA_RUNTIME_ERROR("not implemented");
    }

    librpa_int::global::lib_printf_root("Calculating correlation self-energy by space-time method\n");
    if (!tfg.has_time_grids())
    {
        librpa_int::global::lib_printf_root("Parsed time-frequency object do not have time grids, exiting\n");
        mpi_comm_global_h.barrier();
        throw LIBRPA_RUNTIME_ERROR("input TFGrids object has no time grids");
    }
    mpi_comm_global_h.barrier();

    assert(ad_Wc.initialized());
    // assert(blacs_sigc_h.initialized());
    // blacs_sigc_h_ = blacs_sigc_h;

    const int natom = this->atbasis_wfc.n_atoms;

    std::ofstream ofs_sigmac_r;

#ifndef LIBRPA_USE_LIBRI
    if (mpi_comm_global_h.myid == 0)
    {
        std::cout << "LIBRA::G0W0::build_spacetime is only implemented on top of LibRI" << std::endl;
        std::cout << "Please recompiler LibRPA with -DUSE_LIBRI and configure include path" << std::endl;
    }
    global::mpi_comm_global_h.barrier();
    throw LIBRPA_RUNTIME_ERROR("compilation");
#else
    // Transform from frequency/reciprocal to time/real-space
    global::profiler.start("g0w0_build_spacetime_ct_ft_wc", "Tranform Wc (q,w) -> (R,t)");
    global::profiler.start("g0w0_build_spacetime_ct_ft_real_work", "Perform transformation");
    // global::ofs_myid << Wc_freq_q.at(tfg.get_freq_nodes()[0]) << endl;
    auto Wc_tau_R_blacs = CT_FT_Wc_freq_q(comm_h, Wc_freq_q, pbc, tfg, true, MAJOR::ROW);
    release_free_mem();
    global::profiler.stop("g0w0_build_spacetime_ct_ft_real_work");

    // Build the Wc LibRI object by redistributing the BLACS submatrices
    typedef std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<double>>> tensor_map;
    std::map<double, tensor_map> tau_Wc_libri;
    // NOTE: for case of a few atoms, some process may have much more memory load than others
    global::profiler.start("g0w0_build_spacetime_get_ap", "Compute balance atom-pairs distribution");
    const auto &map_atpairs_balanced = get_balanced_ap_distribution_for_consec_descriptor(
        atbasis_abf, atbasis_abf, ad_Wc);
    global::profiler.stop("g0w0_build_spacetime_get_ap");
    IndexScheduler sched;
    sched.init(map_atpairs_balanced, atbasis_abf, atbasis_abf, ad_Wc, MAJOR::ROW);
    for (auto &[tau, map_R_mat]: Wc_tau_R_blacs)
    {
        for (auto &[R, mat_blacs]: map_R_mat)
        {
            auto pair_mat = get_ap_map_from_blacs_dist_scheduler(mat_blacs, sched, atbasis_abf, atbasis_abf, ad_Wc);
            // if (tau == tfg.get_time_nodes()[0])
            // {
            //     librpa_int::global::ofs_myid << pair_mat << endl;
            // }
            // if (Params::debug)
            // {
            //     std::stringstream ss;
            //     ss << Params::output_dir << "Wc_tau_R"
            //         << "_itau_" << std::setfill('0') << std::setw(5) << tfg.get_time_index(tau)
            //         << "_iR_" << std::setfill('0') << std::setw(5) << get_R_index(Rlist, R) << ".csc";
            //     utils::write_matrix_elsi_csc_parallel(ss.str(), mat_blacs, envs::array_desc_abf_global);
            // }
            global::profiler.start("g0w0_build_spacetime_prep_Wc_all", "Prepare LibRI Wc object");
            for (auto &[pair, mat_ap]: pair_mat)
            {
                const auto &I = as_int(pair.first);
                const auto &J = as_int(pair.second);
                const auto &nabf_I = atbasis_abf.get_atom_nb(I);
                const auto &nabf_J = atbasis_abf.get_atom_nb(J);
                tau_Wc_libri[tau][I][{J, {R.x, R.y, R.z}}] = RI::Tensor<double>({nabf_I, nabf_J}, mat_ap.get_real().sptr());
            }
            // global::ofs_myid << "tau " << tau << endl << tau_Wc_libri[tau] << endl;
            mat_blacs.clear();
            global::profiler.stop("g0w0_build_spacetime_prep_Wc_all");
        }
    }
    global::profiler.stop("g0w0_build_spacetime_ct_ft_wc");

    librpa_int::global::lib_printf_root("Time for Fourier transform of Wc in GW (seconds, Wall/CPU): %f %f\n",
            global::profiler.get_wall_time_last("g0w0_build_spacetime_ct_ft_wc"),
            global::profiler.get_cpu_time_last("g0w0_build_spacetime_ct_ft_wc"));

    RI::GW<int, int, 3, double> gw_libri;
    map<int,std::array<double,3>> atoms_pos;
    // Dummy atoms position
    for (int i = 0; i != natom; i++)
        atoms_pos.insert(pair<int, std::array<double, 3>>{i, {0, 0, 0}});

    global::profiler.start("g0w0_build_spacetime_2", "Setup LibRI G0W0 object and C data");
    gw_libri.set_parallel(ad_Wc.comm(), atoms_pos, pbc.latvec_array, pbc.period_array);
    gw_libri.set_Cs(LRI_Cs.data_libri, this->libri_threshold_C);
    global::profiler.stop("g0w0_build_spacetime_2");
    librpa_int::global::lib_printf_root("Time for LibRI G0W0 setup (seconds, Wall/CPU): %f %f\n",
            global::profiler.get_wall_time_last("g0w0_build_spacetime_2"),
            global::profiler.get_cpu_time_last("g0w0_build_spacetime_2"));

    const auto tot_atpair_ordered = generate_atom_pair_from_nat(natom, true);
    const auto &Rlist = this->pbc.Rlist;
    auto IJR_local_gf = dispatch_vector_prod(tot_atpair_ordered, Rlist, ad_Wc.myid(), ad_Wc.nprocs(), true, false);
    global::ofs_myid << "IJR_local_gf: " << IJR_local_gf << std::endl;

    const int nfreq = tfg.get_n_grids();

    for (auto itau = 0; itau != nfreq; itau++)
    {
        // librpa_int::global::lib_printf("task %d itau %d start\n", mpi_comm_global_h.myid, itau);
        const auto tau = tfg.get_time_nodes()[itau];
        // librpa_int::global::lib_printf("task %d Wc_tau_R.count(tau) %zu\n", mpi_comm_global_h.myid, Wc_tau_R.count(tau));
        auto it = tau_Wc_libri.find(tau);
        const auto &Wc_libri = (it == tau_Wc_libri.end())? tensor_map{} : it->second;
        size_t n_obj_wc_libri = 0;
        for (const auto &w: Wc_libri)
        {
            n_obj_wc_libri += w.second.size();
        }
        global::profiler.start("g0w0_build_spacetime_3", "Setup LibRI Wc");
        gw_libri.set_Ws(Wc_libri, this->libri_threshold_Wc);
        if (it != tau_Wc_libri.end()) tau_Wc_libri.erase(it);
        global::profiler.stop("g0w0_build_spacetime_3");

        for (int ispin = 0; ispin != mf.get_n_spins(); ispin++)
        {
            std::map<int, std::map<std::pair<int, std::array<int, 3>>, Tensor<double>>> sigc_posi_tau;
            std::map<int, std::map<std::pair<int, std::array<int, 3>>, Tensor<double>>> sigc_nega_tau;

            global::profiler.start("g0w0_build_spacetime_4", "Compute G(R,t) and G(R,-t)");
            // global::ofs_myid << "gf size " << gf.size() << endl;
            // global::ofs_myid << "t " << gf[tau].size() << " ; -t " << gf[-tau].size() << endl;
            std::map<double, std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<double>>>> tau_gf_libri;
            const std::vector<double> taus{tau, -tau};
            if (this->is_mf_eigvec_k_distributed_)
            {
                build_gf_libri_kpara(mf, comm_h, atbasis_wfc, ispin, this->pbc.kfrac_list,
                                     taus, IJR_local_gf, tau_gf_libri);
            }
            else
            {
                build_gf_libri_kserial(mf, atbasis_wfc, ispin, this->pbc.kfrac_list,
                                       taus, IJR_local_gf, tau_gf_libri);
            }
            // if (itau == 0 && ispin == 0)  // debug
            // {
            //     global::ofs_myid << tau_gf_libri << std::endl;
            // }
            global::profiler.stop("g0w0_build_spacetime_4");
            librpa_int::global::lib_printf_root("Time for Green's function, i_spin %d i_tau %d (seconds, Wall/CPU): %f %f\n",
                    ispin + 1, itau + 1,
                    global::profiler.get_wall_time_last("g0w0_build_spacetime_4"),
                    global::profiler.get_cpu_time_last("g0w0_build_spacetime_4"));

            for (auto t: taus)
            {
                const auto &gf_libri = tau_gf_libri.at(t);
                size_t n_obj_gf_libri = 0;
                for (const auto &gf: gf_libri)
                    n_obj_gf_libri += gf.second.size();

                double wtime_g0w0_cal_sigc = omp_get_wtime();
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

                wtime_g0w0_cal_sigc = omp_get_wtime() - wtime_g0w0_cal_sigc;
                librpa_int::global::lib_printf("Task %4d. libRI G0W0, spin %1d, time grid %12.6f. Wc size %zu, GF size %zu. Wall time %f\n",
                       ad_Wc.myid(), ispin, t, n_obj_wc_libri, n_obj_gf_libri, wtime_g0w0_cal_sigc);
            }

            size_t n_IJR_myid = 0; // for sigcmat output

            if (this->output_sigc_mat_rt)
            {
                std::stringstream ss;
                ss << path_as_directory(this->output_dir) << "SigcRT"
                    << "_ispin_" << std::setfill('0') << std::setw(5) << ispin
                    << "_itau_" << std::setfill('0') << std::setw(5) << itau
                    << "_myid_" << std::setfill('0') << std::setw(5) << global::myid_global << ".dat";
                ofs_sigmac_r.open(ss.str(), std::ios::out | std::ios::binary);
                ofs_sigmac_r.write((char *) &n_IJR_myid, sizeof(size_t)); // placeholder
            }

            // symmetrize and perform transformation
            global::profiler.start("g0w0_build_spacetime_6", "Transform Sigc (R,t) -> (R,w)");
            for (const auto &I_JR_sigc_posi: sigc_posi_tau)
            {
                const auto &I = I_JR_sigc_posi.first;
                const auto &n_I = this->atbasis_wfc.get_atom_nb(I);
                for (const auto &JR_sigc_posi: I_JR_sigc_posi.second)
                {
                    const auto &J = JR_sigc_posi.first.first;
                    const auto &n_J = this->atbasis_wfc.get_atom_nb(J);
                    const auto &Ra = JR_sigc_posi.first.second;
                    const auto &sigc_posi_block = JR_sigc_posi.second;
                    if (sigc_nega_tau.count(I) == 0 || sigc_nega_tau.at(I).count(JR_sigc_posi.first) == 0)
                    {
                        continue;
                    }
                    const auto &sigc_nega_block = sigc_nega_tau.at(I).at(JR_sigc_posi.first);
                    const auto R = Vector3_Order<int>(Ra[0], Ra[1], Ra[2]);
                    const auto iR = std::distance(Rlist.cbegin(), std::find(Rlist.cbegin(), Rlist.cend(), R));

                    const auto sigc_cos = 0.5 * (sigc_posi_block + sigc_nega_block);
                    const auto sigc_sin = 0.5 * (sigc_posi_block - sigc_nega_block);

                    n_IJR_myid++;
                    if (this->output_sigc_mat_rt)
                    {
                        // write indices and dimensions
                        size_t dims[5];
                        dims[0] = iR;
                        dims[1] = I;
                        dims[2] = J;
                        dims[3] = n_I;
                        dims[4] = n_J;
                        ofs_sigmac_r.write((char *) dims, 5 * sizeof(size_t));
                        // cos: contribute to the real part
                        ofs_sigmac_r.write((char *) sigc_cos.ptr(), n_I * n_J * sizeof(double));
                        // sin: contribute to the imaginary part
                        ofs_sigmac_r.write((char *) sigc_sin.ptr(), n_I * n_J * sizeof(double));
                    }

                    for (size_t iomega = 0; iomega != tfg.get_n_grids(); iomega++)
                    {
                        const auto omega = tfg.get_freq_nodes()[iomega];
                        const auto t2f_sin = tfg.get_sintrans_t2f()(iomega, itau);
                        const auto t2f_cos = tfg.get_costrans_t2f()(iomega, itau);
                        // row-major used here for libRI communication when building sigc_KS
                        Matz sigc_temp(n_I, n_J, MAJOR::ROW);
                        for (size_t i = 0; i != n_I; i++)
                            for (size_t j = 0; j != n_J; j++)
                            {
                                sigc_temp(i, j) = std::complex<double>{sigc_cos(i, j) * t2f_cos, sigc_sin(i, j) * t2f_sin};
                            }
                        atpair_t IJ{I, J};
                        auto &m_IJ = sigc_is_f_IJ_R[ispin][omega][IJ];
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

            sigc_posi_tau.clear();
            sigc_nega_tau.clear();

            global::profiler.stop("g0w0_build_spacetime_6");

            if (this->output_sigc_mat_rt)
            {
                ofs_sigmac_r.seekp(0);
                ofs_sigmac_r.write((char *) &n_IJR_myid, sizeof(size_t)); // overwrite
                ofs_sigmac_r.close();
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
        for (int ispin = 0; ispin != mf.get_n_spins(); ispin++)
        {
            for (size_t iomega = 0; iomega != tfg.get_n_grids(); iomega++)
            {
                size_t n_IJR_myid = 0;
                std::stringstream ss;
                ss << path_as_directory(this->output_dir) << "SigcRF"
                    << "_ispin_" << std::setfill('0') << std::setw(5) << ispin
                    << "_iomega_" << std::setfill('0') << std::setw(5) << iomega
                    << "_myid_" << std::setfill('0') << std::setw(5) << global::myid_global << ".dat";
                ofs_sigmac_r.open(ss.str(), std::ios::out | std::ios::binary);
                ofs_sigmac_r.write((char *) &n_IJR_myid, sizeof(size_t)); // placeholder

                const auto omega = tfg.get_freq_nodes()[iomega];
                const auto &sigc_IJ_R = sigc_is_f_IJ_R[ispin][omega];
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
                        ofs_sigmac_r.write((char *) sigc.ptr(), sigc.size() * 2 * sizeof(double));
                    }
                }
                ofs_sigmac_r.seekp(0);
                ofs_sigmac_r.write((char *) &n_IJR_myid, sizeof(size_t)); // placeholder
                ofs_sigmac_r.close();
            }
        }
    }
}

void G0W0::build_sigc_matrix_KS(const std::map<int, std::map<int, ComplexMatrix>> &wfc_target,
                                const std::vector<Vector3_Order<double>> &kfrac_target,
                                const Atoms &geometry)
{
    throw LIBRPA_RUNTIME_ERROR("Not implemented yet");
}

void G0W0::build_sigc_matrix_KS_blacs(const std::map<int, std::map<int, ComplexMatrix>> &wfc_target,
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
        for (int isp = 0; isp < n_spins; isp++)
        {
            for (const auto& freq: this->tfg.get_freq_nodes())
            {
                std::map<int, std::map<std::pair<int, std::array<int, 3>>, Tensor<complex<double>>>> sigc_I_JR_local;
                auto it_sp = sigc_is_f_IJ_R.find(isp);
                if (it_sp != sigc_is_f_IJ_R.cend())
                {
                    auto it_sp_f = it_sp->second.find(freq);
                    if (it_sp_f != it_sp->second.cend())
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
                ap_p_map<map<Vector3_Order<int>, Matz>> sigc_sp_f_new;
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
                        sigc_sp_f_new[IJ][R] = Matz{n_I, n_J, mat.data, MAJOR::ROW};
                    }
                }
                if (it_sp == this->sigc_is_f_IJ_R.cend())
                    sigc_is_f_IJ_R[isp][freq] = sigc_sp_f_new;
                else
                {
                    auto it_sp_f = it_sp->second.find(freq);
                    if (it_sp_f == it_sp->second.cend())
                        it_sp->second.emplace(freq, sigc_sp_f_new);
                    else
                        it_sp_f->second.swap(sigc_sp_f_new);
                }
                sigc_I_JR.clear();
            }
        }

        is_rspace_redist_for_KS_ = true;
        is_rspace_redist_blacs_ = true;
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
        map<double, map<atom_t, map<atom_t, map<Vector3_Order<int>, Matz>>>> sigc_isp_local;
        for (const auto& freq: this->tfg.get_freq_nodes())
        {
            const auto &sigc_IJ_R = sigc_is_f_IJ_R.at(isp).at(freq);
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

        this->sigc_is_ik_f_KS[isp] = {};

        if (is_mf_eigvec_k_distributed_)
        {
            auto temp_nband_nao = init_local_mat<complex<double>>(desc_nband_nao, MAJOR::COL);
            // collect parsed eigenvectors all all processes
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
                            const auto &wfc_sk = wfc_target.at(isp).at(ik);
                            // processing the k-point eigenvector on this process
                            ScalapackConnector::pgemm_f('C', 'N', n_bands, n_aos, n_aos, 1.0,
                                                        wfc_sk.c, 1, 1, desc_nao_nband_fb.desc,
                                                        sigc_nao_nao.ptr(), 1, 1, desc_nao_nao.desc,
                                                        0.0,
                                                        temp_nband_nao.ptr(), 1, 1, desc_nband_nao.desc);
                            ScalapackConnector::pgemm_f('N', 'N', n_bands, n_bands, n_aos, 1.0,
                                                        temp_nband_nao.ptr(), 1, 1, desc_nband_nao.desc,
                                                        wfc_sk.c, 1, 1, desc_nao_nband_fb.desc,
                                                        0.0,
                                                        sigc_nband_nband_fb.ptr(), 1, 1, desc_nband_nband_fb.desc);
                            this->sigc_is_ik_f_KS[isp][ik][freq] = sigc_nband_nband_fb.copy();
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
        }
        else
        {
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
                    const auto &wfc_isp_k = wfc_target.at(isp).at(ik);
                    blacs_ctxt_h.barrier();
                    const auto wfc_block = get_local_mat(wfc_isp_k.c, MAJOR::ROW, desc_nband_nao, MAJOR::COL).conj();
                    auto temp_nband_nao = multiply_scalapack(wfc_block, desc_nband_nao, sigc_nao_nao, desc_nao_nao, desc_nband_nao);
                    ScalapackConnector::pgemm_f('N', 'C', n_bands, n_bands, n_aos, 1.0,
                                                temp_nband_nao.ptr(), 1, 1, desc_nband_nao.desc,
                                                wfc_block.ptr(), 1, 1, desc_nband_nao.desc, 0.0,
                                                sigc_nband_nband.ptr(), 1, 1, desc_nband_nband.desc);
                    // collect the full matrix to master
                    // TODO: would need a different strategy for large system
                    ScalapackConnector::pgemr2d_f(n_bands, n_bands,
                                                sigc_nband_nband.ptr(), 1, 1, desc_nband_nband.desc,
                                                sigc_nband_nband_fb.ptr(), 1, 1, desc_nband_nband_fb.desc,
                                                desc_nband_nband_fb.ictxt());
                    // NOTE: only the matrices at master process is meaningful
                    sigc_is_ik_f_KS[isp][ik][freq] = sigc_nband_nband_fb.copy();
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

void G0W0::build_sigc_matrix_KS_band(const std::map<int, std::map<int, ComplexMatrix>> &wfc,
                                     const std::vector<Vector3_Order<double>> &kfrac_band,
                                     const Atoms &geometry)
{
    librpa_int::global::mpi_comm_global_h.barrier();
    if (librpa_int::global::mpi_comm_global_h.myid == 0)
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

void G0W0::build_sigc_matrix_KS_band_blacs(const std::map<int, std::map<int, ComplexMatrix>> &wfc,
                                     const std::vector<Vector3_Order<double>> &kfrac_band,
                                     const Atoms &geometry, const BlacsCtxtHandler &blacs_ctxt_h)
{
    librpa_int::global::mpi_comm_global_h.barrier();
    if (librpa_int::global::mpi_comm_global_h.myid == 0)
    {
        librpa_int::global::lib_printf("build_sigc_matrix_KS_kgrid: constructing self-energy matrix for band k-path\n");
    }
    this->build_sigc_matrix_KS_blacs(wfc, kfrac_band, geometry, blacs_ctxt_h);
}

} // namespace librpa_int
