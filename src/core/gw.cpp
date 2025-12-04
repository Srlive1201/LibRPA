#include "gw.h"

#include <fstream>
#include <functional>
#include <map>
#include <sstream>

#include "../io/fs.h"
#include "../io/global_io.h"
#include "../io/stl_io_helper.h"
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
#include "librpa_enums.h"
#include "pbc.h"
#include "ri.h"
#include "utils_atomic_basis_blacs.h"

#ifdef LIBRPA_USE_LIBRI
#include <RI/global/Tensor.h>
#include <RI/physics/GW.h>
using RI::Tensor;
using RI::Communicate_Tensors_Map_Judge::comm_map2_first;
#endif

namespace librpa_int
{

G0W0::G0W0(const MeanField &mf_in, const AtomicBasis &atbasis_wfc_in,
           const PeriodicBoundaryData &pbc_in, const TFGrids &tfg_in, const MpiCommHandler &comm_h_in)
    : mf(mf_in),
      atbasis_wfc(atbasis_wfc_in),
      pbc(pbc_in),
      tfg(tfg_in), comm_h(comm_h_in)
{
    is_rspace_built_ = false;
    is_kspace_built_ = false;
    libri_threshold_C = 0.0;
    libri_threshold_Wc = 0.0;
    libri_threshold_G = 0.0;

    output_sigc_mat = false;
    output_sigc_mat_rt = false;
    output_sigc_mat_rf = false;
}

void G0W0::reset_rspace()
{
    sigc_is_f_R_IJ.clear();
    is_rspace_built_ = false;
}

void G0W0::reset_kspace()
{
    sigc_is_ik_f_KS.clear();
    is_kspace_built_ = false;
}

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
    std::set<Vector3_Order<int>> Rs_local;
    for (const auto &IJR: IJR_local_gf)
    {
        Rs_local.insert(IJR.second);
    }

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
            auto gf = mf.get_gf_real_imagtimes_Rs(ispin, this->pbc.kfrac_list, {tau, -tau}, {Rs_local.cbegin(), Rs_local.cend()});
            // global::ofs_myid << "gf size " << gf.size() << endl;
            // global::ofs_myid << "t " << gf[tau].size() << " ; -t " << gf[-tau].size() << endl;
            std::map<double, std::map<int, std::map<std::pair<int, std::array<int, 3>>, RI::Tensor<double>>>> tau_gf_libri;
            for (auto t: {tau, -tau})
            {
                tau_gf_libri[t] = {};
            }
            for (const auto &IJR: IJR_local_gf)
            {
                const auto &I = IJR.first.first;
                const auto &n_I = this->atbasis_wfc.get_atom_nb(I);
                const auto &J = IJR.first.second;
                const auto &n_J = this->atbasis_wfc.get_atom_nb(J);
                const auto &R = IJR.second;
                for (auto t: {tau, -tau})
                {
                    // skip when this process does not have any local GF data, i.e. Rs_local is empty
                    if (gf.count(t) == 0 || gf.at(t).count(R) == 0) continue;
                    const auto &gf_global = gf.at(t).at(R);
                    matrix gf_IJ_block(n_I, n_J);
                    {
                        for (size_t i = 0; i != n_I; i++)
                            for (size_t j = 0; j != n_J; j++)
                            {
                                gf_IJ_block(i, j) = gf_global(this->atbasis_wfc.get_global_index(I, i),
                                                              this->atbasis_wfc.get_global_index(J, j));
                            }
                        std::shared_ptr<std::valarray<double>> mat_ptr = std::make_shared<std::valarray<double>>(gf_IJ_block.c, gf_IJ_block.size);
                        tau_gf_libri[t][as_int(I)][{as_int(J), {R.x, R.y, R.z}}] = RI::Tensor<double>({n_I, n_J}, mat_ptr);
                    }
                }
            }
            gf.clear();
            global::profiler.stop("g0w0_build_spacetime_4");
            librpa_int::global::lib_printf_root("Time for Green's function, i_spin %d i_tau %d (seconds, Wall/CPU): %f %f\n",
                    ispin + 1, itau + 1,
                    global::profiler.get_wall_time_last("g0w0_build_spacetime_4"),
                    global::profiler.get_cpu_time_last("g0w0_build_spacetime_4"));

            for (auto t: {tau, -tau})
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
                        matrix_m<complex<double>> sigc_temp(n_I, n_J, MAJOR::ROW);
                        for (size_t i = 0; i != n_I; i++)
                            for (size_t j = 0; j != n_J; j++)
                            {
                                sigc_temp(i, j) = std::complex<double>{sigc_cos(i, j) * t2f_cos, sigc_sin(i, j) * t2f_sin};
                            }
                        if (sigc_is_f_R_IJ.count(ispin) == 0 ||
                            sigc_is_f_R_IJ.at(ispin).count(omega) == 0 ||
                            sigc_is_f_R_IJ.at(ispin).at(omega).count(R) == 0 ||
                            sigc_is_f_R_IJ.at(ispin).at(omega).at(R).count(I) == 0 ||
                            sigc_is_f_R_IJ.at(ispin).at(omega).at(R).at(I).count(J) == 0)
                        {
                            sigc_is_f_R_IJ[ispin][omega][R][I][J] = std::move(sigc_temp);
                        }
                        else
                        {
                            sigc_is_f_R_IJ[ispin][omega][R][I][J] += sigc_temp;
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
                ofs_sigmac_r.write((char *) &n_IJR_myid, sizeof(size_t)); // placeholder
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
                const auto &sigc_RIJ = sigc_is_f_R_IJ[ispin][omega];
                for (const auto &[R, IJsigc]: sigc_RIJ)
                {
                    const auto iR = this->pbc.get_R_index(R);
                    for (const auto &[I, Jsigc]: IJsigc)
                    {
                        const auto &n_I = this->atbasis_wfc.get_atom_nb(I);
                        for (const auto &[J, sigc]: Jsigc)
                        {
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
                                const Atoms &geometry,
                                const BlacsCtxtHandler &blacs_ctxt_h)
{
    assert(this->is_rspace_built_);
    if (this->is_kspace_built_)
    {
        global::lib_printf("Warning: reset Sigmac_c k-space matrices\n");
        this->reset_kspace();
    }

    const int n_aos = mf.get_n_aos();
    const int n_bands = mf.get_n_bands();
    const auto period = pbc.period;
    const auto &coord_frac = geometry.coords_frac;

    global::profiler.start("g0w0_build_sigc_KS");

#ifndef LIBRPA_USE_LIBRI
    if (global::mpi_comm_global_h.myid == 0)
    {
        std::cout << "LIBRA::G0W0::build_sigc_matrix_KS is only implemented on top of LibRI" << std::endl;
        std::cout << "Please recompile LibRPA with -DUSE_LIBRI and optionally configure include path" << std::endl;
    }
    global::mpi_comm_global_h.barrier();
    throw LIBRPA_RUNTIME_ERROR("compilation");
#else
    // char fn[80];
    ArrayDesc desc_nband_nao(blacs_ctxt_h);
    desc_nband_nao.init_1b1p(n_bands, n_aos, 0, 0);
    ArrayDesc desc_nao_nao(blacs_ctxt_h);
    desc_nao_nao.init_1b1p(n_aos, n_aos, 0, 0);
    ArrayDesc desc_nband_nband(blacs_ctxt_h);
    desc_nband_nband.init_1b1p(n_bands, n_bands, 0, 0);
    ArrayDesc desc_nband_nband_fb(blacs_ctxt_h);
    desc_nband_nband_fb.init(n_bands, n_bands, n_bands, n_bands, 0, 0);

    // local 2D-block submatrices
    auto sigc_nao_nao = init_local_mat<complex<double>>(desc_nao_nao, MAJOR::COL);
    auto sigc_nband_nband = init_local_mat<complex<double>>(desc_nband_nband, MAJOR::COL);

    const auto set_IJ_nao_nao = get_necessary_IJ_from_block_2D(
        this->atbasis_wfc, this->atbasis_wfc, desc_nao_nao);
    const auto s0_s1 = get_s0_s1_for_comm_map2_first(set_IJ_nao_nao);

    // NOTE: With many MPI tasks, matrix at a particular k point
    // can have contributions from many tasks. The same happens for EXX.
    for (int ispin = 0; ispin < this->mf.get_n_spins(); ispin++)
    {
        for (const auto& freq: this->tfg.get_freq_nodes())
        {
            // Communicate to obtain sub-matrices for necessary I-J pairs at all Rs
            std::map<int, std::map<std::pair<int, std::array<int, 3>>, Tensor<complex<double>>>>
                sigc_I_JR_local;
            if (this->sigc_is_f_R_IJ.count(ispin))
            {
                const auto& sigc_is_freq = this->sigc_is_f_R_IJ.at(ispin).at(freq);
                for (const auto &R_IJ_sigc: sigc_is_freq)
                {
                    const auto R = R_IJ_sigc.first;
                    for (const auto &I_J_sigc: R_IJ_sigc.second)
                    {
                        const auto I = I_J_sigc.first;
                        const auto &n_I = this->atbasis_wfc.get_atom_nb(I);
                        for (const auto &J_sigc: I_J_sigc.second)
                        {
                            const auto J = J_sigc.first;
                            const auto &n_J = this->atbasis_wfc.get_atom_nb(J);
                            const std::array<int, 3> Ra{R.x, R.y, R.z};
                            sigc_I_JR_local[I][{J, Ra}] = Tensor<complex<double>>({n_I, n_J}, J_sigc.second.sptr());
                        }
                    }
                }
            }
            auto sigc_I_JR = comm_map2_first(blacs_ctxt_h.comm(), sigc_I_JR_local, s0_s1.first, s0_s1.second);
            sigc_I_JR_local.clear();
            
            // Convert each <I,<J, R>> pair to the nearest neighbour to speed up later Fourier transform
            // while keep the accuracy in further band interpolation.
            // Reuse the cleared-up sigc_I_JR_local object
            if (coord_frac.size() > 0)
            {
                for (auto &I_sigcJR: sigc_I_JR)
                {
                    const auto &I = I_sigcJR.first;
                    for (auto &JR_sigc: I_sigcJR.second)
                    {
                        const auto &J = JR_sigc.first.first;
                        const auto &R = JR_sigc.first.second;

                        auto distsq = std::numeric_limits<double>::max();
                        Vector3<int> R_IJ;
                        std::array<int, 3> R_bvk;
                        for (int i = -1; i < 2; i++)
                        {
                            R_IJ.x = i * period.x + R[0];
                            for (int j = -1; j < 2; j++)
                            {
                                R_IJ.y = j * period.y + R[1];
                                for (int k = -1; k < 2; k++)
                                {
                                    R_IJ.z = k * period.z + R[2];
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
                                        R_bvk[0] = R_IJ.x;
                                        R_bvk[1] = R_IJ.y;
                                        R_bvk[2] = R_IJ.z;
                                    }
                                }
                            }
                        }
                        sigc_I_JR_local[I][{J, R_bvk}] = std::move(JR_sigc.second);
                    }
                }
            }
            else
            {
                sigc_I_JR_local = std::move(sigc_I_JR);
            }

            // Perform Fourier transform
            for (size_t ik = 0; ik < kfrac_target.size(); ik++)
            {
                const auto kfrac = kfrac_target[ik];

                const std::function<complex<double>(const int &, const std::pair<int, std::array<int, 3>> &)>
                    fourier = [kfrac](const int &I, const std::pair<int, std::array<int, 3>> &J_Ra)
                    {
                        const auto &Ra = J_Ra.second;
                        Vector3<double> R_IJ(Ra[0], Ra[1], Ra[2]);
                        const auto ang = (kfrac * R_IJ) * TWO_PI;
                        return complex<double>{std::cos(ang), std::sin(ang)};
                    };

                sigc_nao_nao.zero_out();
                collect_block_from_IJ_storage_tensor_transform(sigc_nao_nao, desc_nao_nao, 
                        this->atbasis_wfc, this->atbasis_wfc,
                        fourier, sigc_I_JR_local);
                // prepare wave function BLACS
                const auto &wfc_isp_k = wfc_target.at(ispin).at(ik);
                blacs_ctxt_h.barrier();
                const auto wfc_block = get_local_mat(wfc_isp_k.c, MAJOR::ROW, desc_nband_nao, MAJOR::COL).conj();
                auto temp_nband_nao = multiply_scalapack(wfc_block, desc_nband_nao, sigc_nao_nao, desc_nao_nao, desc_nband_nao);
                ScalapackConnector::pgemm_f('N', 'C', n_bands, n_bands, n_aos, 1.0,
                                            temp_nband_nao.ptr(), 1, 1, desc_nband_nao.desc,
                                            wfc_block.ptr(), 1, 1, desc_nband_nao.desc, 0.0,
                                            sigc_nband_nband.ptr(), 1, 1, desc_nband_nband.desc);
                // collect the full matrix to master
                // TODO: would need a different strategy for large system
                auto sigc_nband_nband_fb = init_local_mat<complex<double>>(desc_nband_nband_fb, MAJOR::COL);
                ScalapackConnector::pgemr2d_f(n_bands, n_bands,
                                              sigc_nband_nband.ptr(), 1, 1, desc_nband_nband.desc,
                                              sigc_nband_nband_fb.ptr(), 1, 1, desc_nband_nband_fb.desc,
                                              desc_nband_nband_fb.ictxt());
                // NOTE: only the matrices at master process is meaningful
                sigc_is_ik_f_KS[ispin][ik][freq] = sigc_nband_nband_fb;
            }
        }
    }
#endif
    this->is_kspace_built_ = true;
    global::profiler.stop("g0w0_build_sigc_KS");
}

void G0W0::build_sigc_matrix_KS_kgrid(const Atoms &geometry, const BlacsCtxtHandler &blacs_ctxt_h)
{
    librpa_int::global::mpi_comm_global_h.barrier();
    if (librpa_int::global::mpi_comm_global_h.myid == 0)
    {
        librpa_int::global::lib_printf("build_sigc_matrix_KS_kgrid: constructing self-energy matrix for SCF k-grid\n");
    }
    this->build_sigc_matrix_KS(this->mf.get_eigenvectors(), this->pbc.kfrac_list, geometry, blacs_ctxt_h);
}

void G0W0::build_sigc_matrix_KS_band(const std::map<int, std::map<int, ComplexMatrix>> &wfc,
                                     const std::vector<Vector3_Order<double>> &kfrac_band,
                                     const Atoms &geometry,
                                     const BlacsCtxtHandler &blacs_ctxt_h)
{
    librpa_int::global::mpi_comm_global_h.barrier();
    if (librpa_int::global::mpi_comm_global_h.myid == 0)
    {
        librpa_int::global::lib_printf("build_sigc_matrix_KS_kgrid: constructing self-energy matrix for band k-path\n");
    }
    this->build_sigc_matrix_KS(wfc, kfrac_band, geometry, blacs_ctxt_h);
}

} // namespace librpa_int
