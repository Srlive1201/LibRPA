#include <iostream>
#include <cstring>
#include <ctime>
#include "envs_io.h"
#include "parallel_mpi.h"
#include "envs_mpi.h"
#include "profiler.h"
#include "chi0.h"
#include "libri_utils.h"
#include "stl_io_helper.h"
#include "utils_io.h"
#include "utils_mem.h"
#include "pbc.h"
#include "ri.h"
#include <omp.h>
#include "matrix.h"
#include "complexmatrix.h"
#include "lapack_connector.h"
#include "pbc.h"
#include "constants.h"
#include "params.h"
#include "scalapack_connector.h"
#ifdef LIBRPA_USE_LIBRI
#include <RI/physics/RPA.h>
#endif
#include <array>
#include <map>

using LIBRPA::envs::mpi_comm_global_h;
using LIBRPA::ParallelRouting;
using LIBRPA::parallel_routing;
using LIBRPA::envs::ofs_myid;
using LIBRPA::utils::lib_printf;

void Chi0::build(const Cs_LRI &Cs,
                 const vector<Vector3_Order<int>> &Rlist,
                 const Vector3_Order<int> &R_period,
                 const vector<atpair_t> &atpairs_ABF,
                 const vector<Vector3_Order<double>> &qlist,
                 TFGrids::GRID_TYPES gt, bool use_space_time)
{
    gf_save = gf_discard = 0;
    // reset chi0_q in case the method was called before
    chi0_q.clear();

    // prepare the time-freq grids
    switch (gt)
    {
        case ( TFGrids::GRID_TYPES::Minimax):
        {
            double emin, emax;
            mf.get_E_min_max(emin, emax);
            tfg.generate_minimax(emin, emax);
            break;
        }
        case ( TFGrids::GRID_TYPES::EvenSpaced_TF ):
        {
            // WARN: only used to debug, adjust emin and interval manually
            tfg.generate_evenspaced_tf(0.005, 0.0, 0.005, 0.0);
            break;
        }
        default:
            throw invalid_argument("chi0 on requested grid is not implemented");
    }

    if (use_space_time && !tfg.has_time_grids())
        throw invalid_argument("current grid type does not have time grids");

    if (mpi_comm_global_h.is_root())
        tfg.show();
    mpi_comm_global_h.barrier();

    if (use_space_time)
    {
        // use space-time method
        for ( auto R: Rlist )
            Rlist_gf.push_back(R);
        for ( auto tau: tfg.get_time_nodes() )
        {
            for ( auto R: Rlist_gf)
            {
                build_gf_Rt(R, tau);
                build_gf_Rt(R, -tau);
            }
        }
        if (mpi_comm_global_h.is_root())
        {
            // cout << " Green threshold: " << gf_R_threshold << endl;
            LIBRPA::utils::lib_printf("Finished construction of R-space Green's function\n");
            LIBRPA::utils::lib_printf("| Saved Green's function:     %zu\n", gf_save);
            LIBRPA::utils::lib_printf("| Discarded Green's function: %zu\n\n", gf_discard);
        }
        build_chi0_q_space_time(Cs, R_period, atpairs_ABF, qlist);
    }
    else
    {
        // conventional method does not need to build Green's function explicitly
        build_chi0_q_conventional(Cs, R_period, atpairs_ABF, qlist);
    }
}

void Chi0::build_gf_Rt(Vector3_Order<int> R, double tau)
{
    Profiler::start("cal_Green_func", "space-time Green's function");
    const auto nkpts = mf.get_n_kpoints();
    const auto nspins = mf.get_n_spins();
    const auto nbands = mf.get_n_bands();
    const auto naos = mf.get_n_aos();
    assert (tau != 0);

    // temporary Green's function
    matrix gf_Rt_is_global(naos, naos);

    for (int is = 0; is != nspins; is++)
    {
        gf_Rt_is_global.zero_out();
        auto wg = mf.get_weight()[is];
        if ( tau > 0)
            for (int i = 0; i != wg.size; i++)
            {
                //wg.c[i] = 1.0 / nkpts *nspins - wg.c[i];
                wg.c[i] = 1.0 / nkpts  - wg.c[i]/2*nspins;//
                if (wg.c[i] < 0) wg.c[i] = 0;
            }
        else
        {
            wg *=0.5 *nspins;
        }
        matrix scale(nkpts, nbands);
        // tau-energy phase
        scale = - 0.5 * tau * (mf.get_eigenvals()[is] - mf.get_efermi());
        /* print_matrix("-(e-ef)*tau", scale); */
        for (int ie = 0; ie != scale.size; ie++)
        {
            // NOTE: enforce non-positive phase
            if ( scale.c[ie] > 0) scale.c[ie] = 0;
            scale.c[ie] = std::exp(scale.c[ie]) * wg.c[ie];
        }
        /* print_matrix("exp(-dE*tau)", scale); */
        for (int ik = 0; ik != nkpts; ik++)
        {
            double ang = - klist[ik] * (R * latvec) * TWO_PI;
            complex<double> kphase = complex<double>(cos(ang), sin(ang));
            /* LIBRPA::utils::lib_printf("kphase %f %fj\n", kphase.real(), kphase.imag()); */
            auto scaled_wfc_conj = conj(mf.get_eigenvectors()[is][ik]);
            for ( int ib = 0; ib != nbands; ib++)
                LapackConnector::scal(naos, scale(ik, ib), scaled_wfc_conj.c + naos * ib, 1);
            gf_Rt_is_global += (kphase * transpose(mf.get_eigenvectors()[is][ik], false) * scaled_wfc_conj).real();
        }
        if ( tau < 0 ) gf_Rt_is_global *= -1.;
        omp_lock_t gf_lock;
        omp_init_lock(&gf_lock);  
#pragma omp parallel for schedule(dynamic)
        for ( int I = 0; I != natom; I++)
        {
            const auto I_num = atom_nw[I];
            for (int J = 0; J != natom; J++)
            {
                const size_t J_num = atom_nw[J];
                matrix tmp_green(I_num, J_num);
                for (size_t i = 0; i != I_num; i++)
                {
                    size_t i_glo = atom_iw_loc2glo(I, i);
                    for (size_t j = 0; j != J_num; j++)
                    {
                        size_t j_glo = atom_iw_loc2glo(J, j);
                        tmp_green(i, j) = gf_Rt_is_global(i_glo, j_glo);
                    }
                }
                if (tmp_green.absmax() > gf_R_threshold)
                {
                    // cout<<" max_green_ele:  "<<tmp_green.absmax()<<endl;
                    omp_set_lock(&gf_lock);
                    gf_is_R_tau[is][I][J][R][tau] = std::move(tmp_green);
                    omp_unset_lock(&gf_lock);
                    // cout << is << " " << I << " " << J << " " << R << " " << tau << endl;
                    // print_matrix("gf_is_itau_R[is][I][J][R][itau]", gf_is_R_tau[is][I][J][R][tau]);
                    gf_save++;
                }
                else
                {
                    // LIBRPA::utils::lib_printf("Discarding GF of spin %1d, IJR {%d, %d, (%d, %d, %d)}, t %f\n",
                    //        is, I, J, R.x, R.y, R.z, tau);
                    gf_discard++;
                }
            }
        }
        omp_destroy_lock(&gf_lock);
#pragma omp barrier
    }
    Profiler::stop("cal_Green_func");
}

void Chi0::build_chi0_q_space_time(const Cs_LRI &Cs,
                                   const Vector3_Order<int> &R_period,
                                   const vector<atpair_t> &atpairs_ABF,
                                   const vector<Vector3_Order<double>> &qlist)
{
   // int R_tau_size = Rlist_gf.size() * tfg.size();
    if(parallel_routing == ParallelRouting::LIBRI)
    {
        if (mpi_comm_global_h.is_root())
            cout<<"Use LibRI for chi0"<<endl;
        build_chi0_q_space_time_LibRI_routing(Cs, R_period, atpairs_ABF, qlist);

    }
    else if (parallel_routing == ParallelRouting::R_TAU)
    {
        // if (para_mpi.is_master())
        //     cout << "R_tau_routing" << endl;
        build_chi0_q_space_time_R_tau_routing(Cs, R_period, atpairs_ABF, qlist);
    }
    else
    {
        // if (para_mpi.is_master())
        //     cout << "atom_pair_routing" << endl;
        build_chi0_q_space_time_atom_pair_routing(Cs, R_period, atpairs_ABF, qlist);
    }
}

void Chi0::build_chi0_q_space_time_LibRI_routing(const Cs_LRI &Cs,
                                                 const Vector3_Order<int> &R_period,
                                                 const vector<atpair_t> &atpairs_ABF,
                                                 const vector<Vector3_Order<double>> &qlist)
{
#ifndef LIBRPA_USE_LIBRI
    cout << "LibRI routing requested, but the executable is not compiled with LibRI" << endl;
    cout << "Please recompiler libRPA with -DUSE_LIBRI and configure include path" << endl;
    mpi_comm_global_h.barrier();
    throw std::logic_error("compilation");
#else
    Profiler::start("LibRI_routing", "Loop over LibRI");
    map<int,std::array<double,3>> atoms_pos;
    for(int i=0;i!=atom_mu.size();i++)
        atoms_pos.insert(pair<int,std::array<double,3>>{i,{0,0,0}});
    
    for ( auto freq: tfg.get_freq_nodes() )
        for ( auto q: qlist)
        {
            for ( auto atpair: atpairs_ABF)
            {
                auto Mu = atpair.first;
                auto Nu = atpair.second;
                chi0_q[freq][q][Mu][Nu].create(atom_mu[Mu], atom_mu[Nu]);
            }
        }
    std::array<double,3> xa{latvec.e11,latvec.e12,latvec.e13};
    std::array<double,3> ya{latvec.e21,latvec.e22,latvec.e23};
    std::array<double,3> za{latvec.e31,latvec.e32,latvec.e33};
    std::array<std::array<double,3>,3> lat_array{xa,ya,za};

    std::array<int,3> period_array{R_period.x,R_period.y,R_period.z};

    RI::RPA<int,int,3,double> rpa;
    rpa.set_parallel(mpi_comm_global_h.comm, atoms_pos,lat_array,period_array);
    rpa.set_csm_threshold(Params::libri_chi0_threshold_CSM);

    // local Rlist to collect after chi0s on each process
    auto Rlist_local = dispatch_vector(Rlist_gf, mpi_comm_global_h.myid, mpi_comm_global_h.nprocs, true);
    auto s0_s1 = get_s0_s1_for_comm_map2_first<atom_t, int>(atpairs_ABF);

    Profiler::start("chi0_libri_routing_set_cs", "Set Cs");
    // if (Params::debug)
    //     ofs_myid << Cs_libri;
    // cout << "Setting Cs for rpa object" << endl;
    // LIBRPA::utils::release_free_mem();
    // if(mpi_comm_global_h.is_root())
    // {
    //     printf("Begin set Cs !!! \n");
    //     LIBRPA::utils::display_free_mem();
    //     // printf("chi0_freq_q size: %d,  freq: %f, q:( %f, %f, %f )\n",chi0_wq.size(),freq, q.x,q.y,q.z );
    // }

    rpa.set_Cs(Cs.data_libri, Params::libri_chi0_threshold_C);

    // Cs_libri.clear();
    // LIBRPA::utils::release_free_mem();
    //
    // if(mpi_comm_global_h.is_root())
    // {
    //     printf("After set Cs !!! \n");
    //     LIBRPA::utils::display_free_mem();
    //     // printf("chi0_freq_q size: %d,  freq: %f, q:( %f, %f, %f )\n",chi0_wq.size(),freq, q.x,q.y,q.z );
    // }
    // cout << "Cs of rpa object set" << endl;
    Profiler::stop("chi0_libri_routing_set_cs");

    // dispatch GF accoding to atpair and R
    if (mpi_comm_global_h.is_root())
        LIBRPA::utils::lib_printf("Total count of GFs IJR: %zu\n", get_atom_pair(gf_is_R_tau[0]).size() * Rlist_gf.size());
    mpi_comm_global_h.barrier();
    const auto IJRs_gf_local = dispatch_vector_prod(get_atom_pair(gf_is_R_tau[0]),
                                                    Rlist_gf,
                                                    mpi_comm_global_h.myid, mpi_comm_global_h.nprocs, true, true);
    LIBRPA::utils::lib_printf("| Number of GFs IJR on Proc %4d: %zu\n", mpi_comm_global_h.myid, IJRs_gf_local.size());

    // omp_lock_t lock_chi0_fourier_cosine;
    // omp_init_lock(&lock_chi0_fourier_cosine);
    int count_gf=0;
    for (auto it = 0; it != tfg.size(); it++)
    {
        double tau = tfg.get_time_nodes()[it];
        // cout << tau << " ";
        for(const auto &isp: gf_is_R_tau)
        {
            const auto &gf_IJR_tau = isp.second;
            std::map<int, std::map<std::pair<int,std::array<int,3>>,RI::Tensor<double>>> gf_po_libri;
            std::map<int, std::map<std::pair<int,std::array<int,3>>,RI::Tensor<double>>> gf_ne_libri;
            std::clock_t cpu_clock_start_isp_tau = clock();
            double wtime_start_isp_tau = omp_get_wtime();
            for (const auto &IJR_gf: IJRs_gf_local)
            {
                const auto &I = IJR_gf.first.first;
                const auto &J = IJR_gf.first.second;
                const auto &R = IJR_gf.second;
                std::array<int,3> Ra{R.x,R.y,R.z};
                // cout << I << " " << J << " " << Ra << endl;
                if (gf_IJR_tau.count(I) && gf_IJR_tau.at(I).count(J) && gf_IJR_tau.at(I).at(J).count(R))
                {
                    // positive tau
                    const map<double, matrix> &gf_tau = gf_IJR_tau.at(I).at(J).at(R);
                    if (gf_tau.count(tau) * gf_tau.count(-tau) == 0)
                        continue;
                    // std::valarray<double> mat_po_array(gf_tau.at(tau).c, gf_tau.at(tau).size);
                    ofs_myid << "gf_libri 1,    tau:" << tau << "  I J R: "<<I<<"  "<<J<<"  "<<R<<" size:"<<gf_tau.at(tau).size<<"\n";
                    std::shared_ptr<std::valarray<double>> mat_po_ptr = std::make_shared<std::valarray<double>>(gf_tau.at(tau).c, gf_tau.at(tau).size);
                    // *mat_po_ptr=mat_po_array;
                    gf_po_libri[I][{J, Ra}] = RI::Tensor<double>({size_t(gf_tau.at(tau).nr), size_t(gf_tau.at(tau).nc)}, mat_po_ptr);
                    // negative tau
                    // std::valarray<double> mat_ne_array(gf_tau.at(-tau).c, gf_tau.at(-tau).size);
                    std::shared_ptr<std::valarray<double>> mat_ne_ptr = std::make_shared<std::valarray<double>>(gf_tau.at(-tau).c, gf_tau.at(-tau).size);
                    // *mat_ne_ptr=mat_ne_array;
                    gf_ne_libri[I][{J, Ra}] = RI::Tensor<double>({size_t(gf_tau.at(-tau).nr), size_t(gf_tau.at(-tau).nc)}, mat_ne_ptr);
                    // LIBRPA::fout_para << "gf_libri 2,    tau:" << tau << "  I J R: "<<I<<"  "<<J<<"  "<<R<<"\n";
                    // count_gf+=1;
        //             if(count_gf==500)
        //             {
        //                 printf("gf_libri myid: %d \n",mpi_comm_world_h.myid);
        //                 system("free -m");
        //                 count_gf=0;
        // // printf("chi0_freq_q size: %d,  freq: %f, q:( %f, %f, %f )\n",chi0_wq.size(),freq, q.x,q.y,q.z );
        //             }
                }
            }
            // ofs_myid << "gf_po_libri\n" << gf_po_libri << "\n";
            // ofs_myid << "gf_ne_libri\n" << gf_ne_libri << "\n";
            mpi_comm_global_h.barrier();
            // std::clock_t cpu_clock_done_init_gf = clock();
            ofs_myid << "rpa.cal_chi0s begin,    tau = " << tau << "\n";
            Profiler::start("chi0_libri_routing_cal_chi0s", "Call cal_chi0s");
            rpa.cal_chi0s(gf_po_libri, gf_ne_libri, Params::libri_chi0_threshold_G);
            Profiler::stop("chi0_libri_routing_cal_chi0s");
            ofs_myid << "rpa.cal_chi0s finished, tau = " << tau << "\n";

            // collect chi0 on selected atpairs of all R
            Profiler::start("chi0_libri_routing_collect_Rs", "Collect all R blocks");
            auto chi0s_IJR = RI::Communicate_Tensors_Map_Judge::comm_map2_first(mpi_comm_global_h.comm, rpa.chi0s, s0_s1.first, s0_s1.second);
            Profiler::stop("chi0_libri_routing_collect_Rs");
            std::clock_t cpu_clock_done_chi0s = clock();

            // parse back to chi0
            Profiler::start("chi0_libri_routing_ft_ct", "Fourier and Cosine transform");
            // a simple vector container for OpenMP parallel
            vector<pair<pair<Vector3_Order<double>, int>, atpair_t>> qifreq_atpair_all;
            for (auto q: qlist)
            {
                for (auto ifreq = 0; ifreq < tfg.size(); ifreq++)
                {
                    for (const auto &atpair: atpairs_ABF)
                    {
                        // pre-filter non-existing blocks
                        if (chi0s_IJR.count(atpair.first) == 0) continue;
                        qifreq_atpair_all.push_back({{q, ifreq}, atpair});
                    }
                }
            }
            cout<<"is: "<<isp.first<<" tau: "<<tau<<"  qifreq_atpair_all.size()"<<qifreq_atpair_all.size()<<endl;
            #pragma omp parallel for schedule(dynamic)
            for (int i=0; i<qifreq_atpair_all.size(); ++i)
            {
                auto qifreq_atpair=qifreq_atpair_all[i];
                const auto &q = qifreq_atpair.first.first;
                const auto &ifreq = qifreq_atpair.first.second;
                const double freq = tfg.get_freq_nodes()[ifreq];
                const double trans = tfg.get_costrans_t2f()(ifreq, it);
                const auto &Mu = qifreq_atpair.second.first;
                const auto &Nu = qifreq_atpair.second.second;
                for (const auto & JR_chi0: chi0s_IJR.at(Mu))
                {
                    if (Nu != JR_chi0.first.first) continue;
                    const auto &R = JR_chi0.first.second;
                    Vector3_Order<int> Rint(R[0], R[1], R[2]);
                    const auto &chi0_Rtau = JR_chi0.second;
                    ComplexMatrix cm_chi0(chi0_Rtau.shape[0], chi0_Rtau.shape[1]);
                    for (int i = 0; i < cm_chi0.size; i++)
                        cm_chi0.c[i] = *(chi0_Rtau.ptr()+i);

                    const double arg = q * (Rint * latvec) * TWO_PI;
                    const complex<double> kphase = complex<double>(cos(arg), sin(arg));
                    // if(freq==tfg.get_freq_nodes()[10] && tau == tfg.get_time_nodes()[10])
                    //     cout <<"freq:  "<<freq<<"   Mu Nu:"<<Mu<<", "<<Nu<<";  "<<complex<double>(cos(arg), sin(arg)) << " * " << trans <<"   chi0_tau: "<<cm_chi0(0,0)<<"   chi0_freq:"<<chi0_q[freq][q][Mu][Nu](0,0)<<endl;
                    cm_chi0 *= 2.0 / mf.get_n_spins() * (trans * kphase);
                    // omp_set_lock(&lock_chi0_fourier_cosine);
                    
                    chi0_q[freq][q][Mu][Nu] += cm_chi0;
                    // omp_unset_lock(&lock_chi0_fourier_cosine);
                }
            }
            Profiler::stop("chi0_libri_routing_ft_ct");

            std::clock_t cpu_clock_done_trans = clock();
            double wtime_end_isp_tau = omp_get_wtime();
            if (mpi_comm_global_h.myid == 0)
            {
                lib_printf("chi0s for time point %f, spin %d. CPU time: %f (tensor), %f (trans). Wall time %f\n",
                       tau, isp.first, cpu_time_from_clocks_diff(cpu_clock_start_isp_tau, cpu_clock_done_chi0s),
                       cpu_time_from_clocks_diff(cpu_clock_done_chi0s, cpu_clock_done_trans),
                       wtime_end_isp_tau - wtime_start_isp_tau);
            }
        }
    }
    // omp_destroy_lock(&lock_chi0_fourier_cosine);

    if (mpi_comm_global_h.is_root()) lib_printf("\n");
    mpi_comm_global_h.barrier();
    Profiler::stop("LibRI_routing");
#endif
}


void Chi0::build_chi0_q_space_time_R_tau_routing(const Cs_LRI &Cs,
                                                 const Vector3_Order<int> &R_period,
                                                 const vector<atpair_t> &atpairs_ABF,
                                                 const vector<Vector3_Order<double>> &qlist)
{
    assert (!Cs.use_libri);
    const auto &LRI_Cs = Cs.data_IJR;

    Profiler::start("R_tau_routing", "Loop over R-tau");
    // taus and Rs to compute on MPI task
    // tend to calculate more Rs on one process
    vector<pair<int, int>> itauiRs_local = dispatcher(0, tfg.size(), 0, Rlist_gf.size(),
                                                      mpi_comm_global_h.myid, mpi_comm_global_h.nprocs, true, false);
    map<Vector3_Order<double>,int> qlist2myid;
    auto loc_qlist = dispatch_vector(qlist , mpi_comm_global_h.myid, mpi_comm_global_h.nprocs, true);
    for(int id=0;id!=mpi_comm_global_h.nprocs;id++)
    {
        auto id_qlist=dispatch_vector(qlist , id, mpi_comm_global_h.nprocs, true);
        for(auto &id_q:id_qlist)
            qlist2myid.insert(std::make_pair(id_q,id));
    }


    map<double, map<Vector3_Order<double>, atom_mapping<ComplexMatrix>::pair_t_old>> chi0_q_tmp;
    for ( auto freq: tfg.get_freq_nodes() )
        for ( auto q: qlist)
        {
            for ( auto atpair: atpairs_ABF)
            {
                auto Mu = atpair.first;
                auto Nu = atpair.second;
                chi0_q_tmp[freq][q][Mu][Nu].create(atom_mu[Mu], atom_mu[Nu]);
            }
        }

    omp_lock_t chi0_lock;
    omp_init_lock(&chi0_lock);
    // double t_chi0_begin = omp_get_wtime();
    // double t_chi0_tot = 0;
#pragma omp parallel for schedule(dynamic)
    for ( int i=0;i!= itauiRs_local.size();i++)
    {
        auto itau = itauiRs_local[i].first;
        auto tau = tfg.get_time_nodes()[itau];
        auto iR = itauiRs_local[i].second;
        auto R = Rlist_gf[iR];
        double t_Rtau_begin = omp_get_wtime();
        for ( auto atpair: atpairs_ABF)
        {
            atom_t Mu = atpair.first;
            atom_t Nu = atpair.second;
            for ( int is = 0; is != mf.get_n_spins(); is++ )
            {
                // double chi0_ele_begin = omp_get_wtime();
                matrix chi0_tau;
                /* if (itau == 0) // debug first itau */
                    chi0_tau = 2.0 /mf.get_n_spins() * compute_chi0_s_munu_tau_R(LRI_Cs, R_period, is, Mu, Nu, tau, R);
                // print_matrix("", chi0_tau);
                /* else continue; // debug first itau */
                // double chi0_ele_t = omp_get_wtime() - chi0_ele_begin;
                omp_set_lock(&chi0_lock);
                for ( auto q: qlist )
                {
                    double arg = q * (R * latvec) * TWO_PI;
                    const complex<double> kphase = complex<double>(cos(arg), sin(arg));
                    for ( int ifreq = 0; ifreq != tfg.size(); ifreq++ )
                    {
                        double freq = tfg.get_freq_nodes()[ifreq];
                        double trans = tfg.get_costrans_t2f()(ifreq, itau);
                        const complex<double> weight = trans * kphase;
                        /* cout << weight << endl; */
                        // if(freq==tfg.get_freq_nodes()[10] && tau == tfg.get_time_nodes()[10])
                        //     cout <<"freq:  "<<freq<<"   Mu Nu:"<<Mu<<", "<<Nu<<";  "<<complex<double>(cos(arg), sin(arg)) << " * " << trans <<"   chi0_tau: "<<chi0_tau(0,0)<<"   chi0_freq:"<<chi0_q_tmp[freq][q][Mu][Nu](0,0)<<endl;
                        chi0_q_tmp[freq][q][Mu][Nu] += ComplexMatrix(chi0_tau) * weight;
                    }
                }
                omp_unset_lock(&chi0_lock);
            }
        }
        double t_Rtau_end = omp_get_wtime();
        double time_used = t_Rtau_end - t_Rtau_begin;
        // t_chi0_tot += time_used;
        lib_printf("CHI0 p_id: %3d, thread: %3d, R: ( %d,  %d,  %d ), tau: %f , TIME_USED: %f\n",
               mpi_comm_global_h.myid, omp_get_thread_num(), R.x, R.y, R.z, tau, time_used);
    }

    omp_destroy_lock(&chi0_lock);
    // Reduce to MPI
#pragma omp barrier
    for ( int ifreq = 0; ifreq < tfg.size(); ifreq++ )
    {
        double freq = tfg.get_freq_nodes()[ifreq];
        for ( auto q: qlist )
        {
            int id_contain_q = qlist2myid[q];
            for ( auto atpair: atpairs_ABF)
            {
                auto Mu = atpair.first;
                auto Nu = atpair.second;
                ComplexMatrix tmp_chi0_recv(atom_mu[Mu], atom_mu[Nu]);
                tmp_chi0_recv.zero_out();
                
                mpi_comm_global_h.barrier();
                /* cout << "nr/nc chi0_q_tmp: " << chi0_q_tmp[ifreq][iq][Mu][Nu].nr << ", "<< chi0_q_tmp[ifreq][iq][Mu][Nu].nc << endl; */
                /* cout << "nr/nc chi0_q: " << chi0_q[ifreq][iq][Mu][Nu].nr << ", "<< chi0_q[ifreq][iq][Mu][Nu].nc << endl; */
                mpi_comm_global_h.reduce_ComplexMatrix(chi0_q_tmp[freq][q][Mu][Nu], tmp_chi0_recv, id_contain_q);
                if(id_contain_q == mpi_comm_global_h.myid )
                    chi0_q[freq][q][Mu][Nu]=std::move(tmp_chi0_recv);
                /* if (LIBRPA::mpi_comm_global_h.myid==0 && Mu == 0 && Nu == 0 && ifreq == 0 && q == Vector3_Order<double>{0, 0, 0}) */
                /* if (LIBRPA::mpi_comm_global_h.myid==0 && Mu == 0 && Nu == 0 && ifreq == 0 ) */
                /* if (LIBRPA::mpi_comm_global_h.myid==0 && ifreq == 0 && q == Vector3_Order<double>{0, 0, 0}) */
                /* { */
                /*     cout <<  "freq: " << freq << ", q: " << q << endl; */
                /*     lib_printf("Mu %zu Nu %zu\n", Mu, Nu); */
                /*     print_complex_matrix("chi0_q ap, first freq first q", chi0_q[freq][q][Mu][Nu]); */
                /* } */
                chi0_q_tmp[freq][q][Mu].erase(Nu);
            }
        }
    }
    Profiler::stop("R_tau_routing");
}

void Chi0::build_chi0_q_space_time_atom_pair_routing(const Cs_LRI &Cs,
                                                     const Vector3_Order<int> &R_period,
                                                     const vector<atpair_t> &atpairs_ABF,
                                                     const vector<Vector3_Order<double>> &qlist)
{
    Profiler::start("atom_pair_routing", "Loop over atom pairs");
    //auto tot_pair = dispatch_vector(atpairs_ABF, LIBRPA::mpi_comm_global_h.myid, para_mpi.get_size(), false);
    lib_printf("Number of atom pairs on Proc %4d: %zu\n", mpi_comm_global_h.myid, atpairs_ABF.size());
    mpi_comm_global_h.barrier();
    omp_lock_t chi0_lock;
    omp_init_lock(&chi0_lock);
    double t_chi0_begin = omp_get_wtime();

    assert (!Cs.use_libri);
    const auto &LRI_Cs = Cs.data_IJR;

#pragma omp parallel
    {
#pragma omp for schedule(dynamic)
        for (size_t atom_pair = 0; atom_pair != atpairs_ABF.size(); atom_pair++)
        {
            auto Mu = atpairs_ABF[atom_pair].first;
            auto Nu = atpairs_ABF[atom_pair].second;
            // const size_t mu_num = atom_mu[Mu];
            // const size_t nu_num = atom_mu[Nu];
            double task_begin = omp_get_wtime();
            map<double, map<Vector3_Order<double>, ComplexMatrix>> tmp_chi0_freq_k;
            for ( auto freq: tfg.get_freq_nodes() )
                for ( auto q: qlist)
                {
                    tmp_chi0_freq_k[freq][q].create(atom_mu[Mu], atom_mu[Nu]);
                }

            for ( int is = 0; is != mf.get_n_spins(); is++ )
                for (auto &R : Rlist_gf)
                {
                    for (auto it = 0; it != tfg.size(); it++)
                    {
                        double tau = tfg.get_time_nodes()[it];
                        ComplexMatrix tmp_chi0_tau(2.0 /mf.get_n_spins() * ComplexMatrix(compute_chi0_s_munu_tau_R(LRI_Cs, R_period, is, Mu, Nu, tau, R)));
                        // print_complex_matrix("", tmp_chi0_tau);
                        for (auto &q : qlist)
                        {
                            const double arg = (q * (R * latvec)) * TWO_PI;
                            const complex<double> kphase = complex<double>(cos(arg), sin(arg));
                            for (auto ifreq = 0; ifreq != tfg.size(); ifreq++)
                            {
                                double freq = tfg.get_freq_nodes()[ifreq];
                                double trans = tfg.get_costrans_t2f()(ifreq, it);
                                const complex<double> weight = trans * kphase;
                                /* const complex<double> cos_weight_kpashe = kphase * tfg.get_costrans_t2f()[ifreq, it]; */
                                //  cout<<"  tmp_chi0  nr nc:  "<<tmp_chi0_freq_k[ifreq][ik_vec].nr<<"  "<<tmp_chi0_freq_k[ifreq][ik_vec].nc<<"    chi0_tau: "<<tmp_chi0_tau.nr<<"  "<<tmp_chi0_tau.nr<<endl;
                                tmp_chi0_freq_k[freq][q] += tmp_chi0_tau * weight;
                            }
                        }
                    }
                    // vector<matrix> chi0_freq_tmp(tmp_cosine_tran(freq_grid.size(),chi0_tau_tmp));
                }

            double task_end = omp_get_wtime();
            double time_used = task_end - task_begin;
            /* time_task_tot += time_used; */
            omp_set_lock(&chi0_lock);
            for (auto &freq : tfg.get_freq_nodes() )
            {
                for (auto &q : qlist)
                    chi0_q[freq][q][Mu][Nu] = std::move(tmp_chi0_freq_k[freq][q]);
            }
            double add_end = omp_get_wtime();
            double add_time = add_end - task_end;
            lib_printf("CHI0 p_id: %3d, thread: %3d, I: %zu, J: %zu, move time: %f  TIME_USED: %f\n", mpi_comm_global_h.myid, omp_get_thread_num(), Mu, Nu, add_time, time_used);
            omp_unset_lock(&chi0_lock);
        }
    }
    mpi_comm_global_h.barrier();
    double t_chi0_end= omp_get_wtime();
    Profiler::stop("atom_pair_routing");
    if(mpi_comm_global_h.is_root())
        lib_printf("| total chi0 time: %f\n",t_chi0_end-t_chi0_begin);
}

matrix Chi0::compute_chi0_s_munu_tau_R(const atpair_R_mat_t &Cs_IJR,
                                       const Vector3_Order<int> &R_period, 
                                       int spin_channel,
                                       atom_t Mu, atom_t Nu, double tau, Vector3_Order<int> R)
{
    Profiler::start("cal_chi0_element", "chi(tau,R,I,J)");
    /* lib_printf("     begin chi0  thread: %d,  I: %zu, J: %zu\n",omp_get_thread_num(), Mu, Nu); */

    assert ( tau > 0 );
    // Local RI requires
    const atom_t I_index = Mu;
    const atom_t J_index = Nu;

    const size_t i_num = atom_nw[Mu];
    const size_t mu_num = atom_mu[Mu];

    const size_t j_num = atom_nw[Nu];
    const size_t nu_num = atom_mu[Nu];

    int flag_G_IJRt = 0;
    int flag_G_IJRNt = 0;

    /* lib_printf("     check if already calculated\n"); */
    /* lib_printf("     size of Green_atom: %zu\n", Green_atom.size()); */
    const auto & gf_R_tau = gf_is_R_tau.at(spin_channel);
    if (gf_R_tau.at(I_index).count(J_index))
        if (gf_R_tau.at(I_index).at(J_index).count(R))
        {
            if (gf_R_tau.at(I_index).at(J_index).at(R).count(tau))
                flag_G_IJRt = 1;

            if (gf_R_tau.at(I_index).at(J_index).at(R).count(-tau))
                flag_G_IJRNt = 1;
        }

    matrix X_R2(i_num, j_num * nu_num);
    matrix X_conj_R2(i_num, j_num * nu_num);

    for (const auto &L_pair : Cs_IJR.at(J_index))
    {
        const auto L_index = L_pair.first;
        const size_t l_num = atom_nw[L_index];
        /* LIBRPA::utils::lib_printf("     begin is loop\n"); */
        if (gf_R_tau.at(I_index).count(L_index))
        {
            for (const auto &R2_index : L_pair.second)
            {
                const auto R2 = R2_index.first;
                const auto &Cs_mat2 = R2_index.second;

                Vector3_Order<int> R_temp_2(Vector3_Order<int>(R2 + R) % R_period);

                if (gf_R_tau.at(I_index).at(L_index).count(R_temp_2))
                {
                    assert(j_num * l_num == (*Cs_mat2).nr);
                    /* LIBRPA::utils::lib_printf("          thread: %d, X_R2 IJL:   %zu,%zu,%zu  R:(  %d,%d,%d  )  tau:%f\n",omp_get_thread_num(),I_index,J_index,L_index,R.x,R.y,R.z,time_tau); */
                    matrix Cs2_reshape(reshape_Cs(j_num, l_num, nu_num, Cs_mat2));
                    
                    if (gf_R_tau.at(I_index).at(L_index).at(R_temp_2).count(tau))
                    {
                        // cout<<"C";
                        // Profiler::start("X");
                        X_R2 += gf_R_tau.at(I_index).at(L_index).at(R_temp_2).at(tau) * Cs2_reshape;
                        // Profiler::stop("X");
                    }

                    if (gf_R_tau.at(I_index).at(L_index).at(R_temp_2).count(-tau))
                    {
                        // cout<<"D";
                        // Profiler::start("X");
                        X_conj_R2 += gf_R_tau.at(I_index).at(L_index).at(R_temp_2).at(-tau) * Cs2_reshape;
                        // Profiler::stop("X");
                    }
                    
                }
            }
        }
    }
    matrix X_R2_rs(reshape_mat(i_num, j_num, nu_num, X_R2));
    matrix X_conj_R2_rs(reshape_mat(i_num, j_num, nu_num, X_conj_R2));
    
    matrix O_sum(mu_num, nu_num);
    for (const auto &K_pair : Cs_IJR.at(I_index))
    {
        const auto K_index = K_pair.first;
        const size_t k_num = atom_nw[K_index];

        for (const auto &R1_index : K_pair.second)
        {
            const auto R1 = R1_index.first;
            const auto &Cs_mat1 = R1_index.second;
            // cout<<"R1:  begin  "<<R1<<endl;
                // matrix X_R2(i_num,j_num*nu_num);
                // matrix X_conj_R2(i_num,j_num*nu_num);
                matrix O(i_num, k_num * nu_num);
                matrix Z(k_num, i_num * nu_num);

                if (flag_G_IJRt || flag_G_IJRNt)
                {
                    matrix N_R2(k_num, j_num * nu_num);
                    matrix N_conj_R2(k_num, j_num * nu_num);
                    for (const auto &L_pair : Cs_IJR.at(J_index))
                    {
                        const auto L_index = L_pair.first;
                        const size_t l_num = atom_nw[L_index];
                        if (gf_R_tau.at(K_index).count(L_index))
                        {
                            for (const auto &R2_index : L_pair.second)
                            {
                                const auto R2 = R2_index.first;
                                const auto &Cs_mat2 = R2_index.second;
                                Vector3_Order<int> R_temp_1(Vector3_Order<int>(R + R2 - R1) % R_period);
                                // Vector3_Order<int> R_temp_2(R2+R);
                                if (gf_R_tau.at(K_index).at(L_index).count(R_temp_1))
                                {
                                    
                                    assert(j_num * l_num == (*Cs_mat2).nr);
                                    // LIBRPA::utils::lib_printf("          thread: %d, IJKL:   %d,%d,%d,%d  R:(  %d,%d,%d  )  tau:%f\n",omp_get_thread_num(),I_index,J_index,K_index,L_index,R.x,R.y,R.z,time_tau);
                                    matrix Cs2_reshape(reshape_Cs(j_num, l_num, nu_num, Cs_mat2));
                                    if (flag_G_IJRNt && gf_R_tau.at(K_index).at(L_index).at(R_temp_1).count(tau))
                                    {
                                        // cout<<"A";
                                        // Profiler::start("N");
                                        N_R2 += gf_R_tau.at(K_index).at(L_index).at(R_temp_1).at(tau) * Cs2_reshape;
                                        // Profiler::stop("N");
                                    }
                                    if (flag_G_IJRt && gf_R_tau.at(K_index).at(L_index).at(R_temp_1).count(-tau))
                                    {
                                        // cout<<"B";
                                        // Profiler::start("N");
                                        N_conj_R2 += gf_R_tau.at(K_index).at(L_index).at(R_temp_1).at(-tau) * Cs2_reshape;
                                        // Profiler::stop("N");
                                    }
                                    
                                }
                            }
                        }
                    }
                    if (flag_G_IJRt)
                    {
                        // Profiler::start("O");
                        matrix N_conj_R2_rs(reshape_mat(k_num, j_num, nu_num, N_conj_R2));
                        O += gf_R_tau.at(I_index).at(J_index).at(R).at(tau) * N_conj_R2_rs;
                        // Profiler::stop("O");
                    }
                    if (flag_G_IJRNt)
                    {
                        // Profiler::start("O");
                        matrix N_R2_rs(reshape_mat(k_num, j_num, nu_num, N_R2));
                        O += gf_R_tau.at(I_index).at(J_index).at(R).at(-tau) * N_R2_rs;
                        // Profiler::stop("O");
                    }
                }
                Vector3_Order<int> R_temp_3(Vector3_Order<int>(R - R1) % R_period);
                if (gf_R_tau.at(K_index).count(J_index))
                {
                    if (gf_R_tau.at(K_index).at(J_index).count(R_temp_3))
                    {
                        if (gf_R_tau.at(K_index).at(J_index).at(R_temp_3).count(-tau))
                        {
                            // Profiler::start("Z");
                            Z += gf_R_tau.at(K_index).at(J_index).at(R_temp_3).at(-tau) * X_R2_rs;
                            // Profiler::stop("Z");
                        }

                        if (gf_R_tau.at(K_index).at(J_index).at(R_temp_3).count(tau))
                        {
                            // Profiler::start("Z");
                            Z += gf_R_tau.at(K_index).at(J_index).at(R_temp_3).at(tau) * X_conj_R2_rs;
                            // Profiler::stop("Z");
                        }
                    }
                }

                matrix Z_rs(reshape_mat(k_num, i_num, nu_num, Z));

                O += Z_rs;
                matrix OZ(reshape_mat_21(i_num, k_num, nu_num, O));
                matrix Cs1_tran(transpose(*Cs_mat1));
                // Profiler::start("O");
                O_sum += Cs1_tran * OZ;
                // Profiler::stop("O");
                // cout<<"   K, R1:   "<<K_index<<"   "<<R1;
                // rt_m_max(O_sum);
        }
    }
    
    /* if ( LIBRPA::mpi_comm_global_h.myid==0 && Mu == 1 && Nu == 1 && tau == tfg.get_time_nodes()[0]) */
    /* { */
    /*     cout << R << " Mu=" << Mu << " Nu=" << Nu << " tau:" << tau << endl; */
    /*     print_matrix("space-time chi0", O_sum); */
    /* } */
    // if (Params::debug)
    // {
    //     cout << R << " Mu=" << Mu << " Nu=" << Nu << " tau=" << tau << " R=" << R << endl;
    //     print_matrix("space-time chi0", O_sum);
    // }
    Profiler::stop("cal_chi0_element");
    return O_sum;
}

void Chi0::build_chi0_q_conventional(const Cs_LRI &Cs,
                                     const Vector3_Order<int> &R_period,
                                     const vector<atpair_t> &atpair_ABF,
                                     const vector<Vector3_Order<double>> &qlist)
{
    // TODO: low implementation priority
    throw logic_error("Not implemented");
}

/*!
 * @warning ZMY: Since Green's function are currently stored in a real matrix,
 *          the complex conjugate in the equations are in fact not taken.
 *          But the result seems okay for the test cases, compared to FHI-aims.
 *          Definitely should be resolved in the future.
 */
map<size_t, matrix> compute_chi0_munu_tau_LRI_saveN_noreshape(const map<size_t, atom_mapping<matrix>::pair_t_old> &gf_occ_ab_t,
                                                              const map<size_t, atom_mapping<matrix>::pair_t_old> &gf_unocc_ab_t,
                                                              const atpair_R_mat_t &LRI_Cs,
                                                              const vector<Vector3_Order<int>> &Rlist, const Vector3_Order<int> &R_period,
                                                              const vector<int> iRs,
                                                              atom_t Mu, atom_t Nu)
{
    map<size_t, matrix> chi0_tau;

    // Store N and N* so one does not need to recompute for each R of chi
    // the extra memory consumption scales as Cs/nR/natoms,
    // which is a small amount compared to Cs for either small (small natoms, large nR) or large (small nR, large natoms) system
    // N(tau) at R and atom K, at(iR, iK)
    map<pair<Vector3_Order<int>, atom_t>, matrix> N;
    // N*(-tau) at R and atom K, at(iR, iK)
    map<pair<Vector3_Order<int>, atom_t>, matrix> Ncnt; // c -> conjugate, nt -> negative time

    const atom_t iI = Mu;
    const atom_t iJ = Nu;
    const size_t n_i = atom_nw[iI];
    const size_t n_mu = atom_mu[Mu];
    const size_t n_j = atom_nw[iJ];
    const size_t n_nu = atom_mu[Nu];

    // pre-compute N and N*
    // only need to calculate N on R-R1, with R in Green's function and R1 from Cs
    // TODO: may reuse the N code to compute N*, by parsing G(t) and G(-t) in same map as done by Shi Rong
    // TODO: decouple N and N*, reduce memory usage
    cout << "Mu: " << Mu << ", Nu: " << Nu << endl;
    for ( auto iR: iRs )
    {
        /* cout << "iR: " << iR << endl; */
        for ( auto const & iK_R1_Cs: LRI_Cs.at(Mu) )
        {
            auto iK = iK_R1_Cs.first;
            /* cout << "iK: " << iK << endl; */
            /* cout << "iR: " << iR << ", iK: " << iK << endl; */
            const size_t n_k = atom_nw[iK];
            for ( auto const & R1_Cs: iK_R1_Cs.second )
            {
                auto const & R1 = R1_Cs.first;
                auto const & RmR1 = Vector3_Order<int>(Rlist[iR] - R1);
                /* cout << "iRmR1: " << iRmR1 << endl; */
                /* cout << "Testing N.count(iRmR1, iK): " << N.count({iRmR1, iK}) << endl; */
                /* cout << "Testing gf_occ_ab_t.count(iR): " << gf_occ_ab_t.count(iR) << endl; */
                // N part, do not recompute and should have the GF counterpart
                if ( N.count({RmR1, iK}) == 0 && gf_occ_ab_t.count(iR) )
                {
                    N[{RmR1, iK}].create(n_nu, n_j*n_k);
                    matrix & N_iRiK = N[{RmR1, iK}];
                    for (auto const & iL_R2_Cs: LRI_Cs.at(Nu))
                    {
                        auto const & iL = iL_R2_Cs.first;
                        const size_t n_l = atom_nw[iL];
                        for ( auto const & R2_Cs: iL_R2_Cs.second )
                        {
                            auto const & R2 = R2_Cs.first;
                            auto mRpR1mR2 = Vector3_Order<int>(-RmR1-R2) % R_period;
                            int i_mRpR1mR2 = get_R_index(Rlist, mRpR1mR2);
                            if ( gf_unocc_ab_t.count(i_mRpR1mR2) )
                            {
                                if (gf_unocc_ab_t.at(i_mRpR1mR2).count(iL) == 0 ||
                                    gf_unocc_ab_t.at(i_mRpR1mR2).at(iL).count(iK) == 0) continue;
                                const matrix & gf_unocc = gf_unocc_ab_t.at(i_mRpR1mR2).at(iL).at(iK);
                                auto const & Cs_nu_jlR2 = R2_Cs.second;
                                assert( Cs_nu_jlR2->nr == n_j * n_l && Cs_nu_jlR2->nc == n_nu );
                                matrix tran_Cs_nu_jlR2 = transpose(*Cs_nu_jlR2);
                                assert( gf_unocc.nr == n_l && gf_unocc.nc == n_k );
                                for ( int i_nu = 0; i_nu != n_nu; i_nu++ )
                                    LapackConnector::gemm('N', 'N', n_j, n_k, n_l, 1.0,
                                                          tran_Cs_nu_jlR2.c+i_nu*n_j*n_l, n_l,
                                                          gf_unocc.c, gf_unocc.nc,
                                                          1.0, N_iRiK.c + i_nu*n_j*n_k, n_k);
                            }
                        }
                    }
                }
                // N* part
                if ( !Ncnt.count({RmR1, iK}) && gf_unocc_ab_t.count(iR) )
                {
                    Ncnt[{RmR1, iK}].create(n_nu, n_j*n_k);
                    matrix & Ncnt_iRiK = Ncnt.at({RmR1, iK});
                    for (auto const & iL_R2_Cs: LRI_Cs.at(Nu))
                    {
                        auto const & iL = iL_R2_Cs.first;
                        const size_t n_l = atom_nw[iL];
                        for ( auto const & R2_Cs: iL_R2_Cs.second )
                        {
                            auto const & R2 = R2_Cs.first;
                            auto mRpR1mR2 = Vector3_Order<int>(-RmR1-R2) % R_period;
                            int i_mRpR1mR2 = get_R_index(Rlist, mRpR1mR2);
                            if ( gf_occ_ab_t.count(i_mRpR1mR2) )
                            {
                                if (gf_occ_ab_t.at(i_mRpR1mR2).count(iL) == 0 ||
                                    gf_occ_ab_t.at(i_mRpR1mR2).at(iL).count(iK) == 0) continue;
                                const matrix & gf_occ = gf_occ_ab_t.at(i_mRpR1mR2).at(iL).at(iK);
                                auto const & Cs_nu_jlR2 = R2_Cs.second;
                                assert( Cs_nu_jlR2->nr == n_j * n_l && Cs_nu_jlR2->nc == n_nu );
                                matrix tran_Cs_nu_jlR2 = transpose(*Cs_nu_jlR2);
                                assert( gf_occ.nr == n_l && gf_occ.nc == n_k );
                                for ( int i_nu = 0; i_nu != n_nu; i_nu++ )
                                    LapackConnector::gemm('N', 'N', n_j, n_k, n_l, 1.0,
                                                          tran_Cs_nu_jlR2.c+i_nu*n_j*n_l, n_l,
                                                          gf_occ.c, n_k,
                                                          1.0, Ncnt_iRiK.c + i_nu*n_j*n_k, n_k);
                            }
                        }
                    }
                }
            }
        }
    }
    cout << "Size N: " << N.size() << endl;
    /* cout << "Precalcualted N, size (iR1 x iK): " << N.size() << endl; */

    for ( auto iR: iRs )
    {
        chi0_tau[iR].create(n_mu, n_nu); // TODO: selectively create by G

        for ( auto const & iK_R1_Cs: LRI_Cs.at(Mu) )
        {
            auto const iK = iK_R1_Cs.first;
            const size_t n_k = atom_nw[iK];
            for (auto const & R1_Cs: iK_R1_Cs.second)
            {
                const auto & Cs = R1_Cs.second;
                // a temporary matrix to store the sum of M(t) and M*(-t)
                matrix MpMc(n_nu, n_i*n_k);
                auto RmR1 = Vector3_Order<int>(Rlist[iR]-R1_Cs.first);
                // M(t)=G(t)N(t)
                if ( gf_occ_ab_t.count(iR) && N.count({RmR1, iK}))
                {
                    if ( gf_occ_ab_t.at(iR).count(iI) != 0 &&
                         gf_occ_ab_t.at(iR).at(iI).count(iJ) != 0)
                    {
                        /* cout << "Computing GN" << endl; */
                        const matrix & gf = gf_occ_ab_t.at(iR).at(iI).at(iJ);
                        const matrix & _N = N.at({RmR1, iK});
                        assert( _N.nr == n_nu && _N.nc == n_j * n_k);
                        for ( int i_nu = 0; i_nu != n_nu; i_nu++ )
                            LapackConnector::gemm('N', 'N', n_i, n_k, n_j, 1.0,
                                                  gf.c, n_j,
                                                  _N.c+i_nu*n_j*n_k, n_k,
                                                  1.0, MpMc.c+i_nu*n_i*n_k, n_k);
                    }
                }
                // M*(-t)=G*(-t)N*(-t)
                if ( gf_unocc_ab_t.count(iR) && N.count({RmR1, iK}))
                {
                    if ( gf_unocc_ab_t.at(iR).count(iI) != 0 &&
                         gf_unocc_ab_t.at(iR).at(iI).count(iJ) != 0)
                    {
                        /* cout << "Computing G*N*" << endl; */
                        const matrix & gf = gf_unocc_ab_t.at(iR).at(iI).at(iJ);
                        const matrix & _Ncnt = Ncnt.at({RmR1, iK});
                        assert( _Ncnt.nr == n_nu && _Ncnt.nc == n_j * n_k);
                        for ( int i_nu = 0; i_nu != n_nu; i_nu++ )
                            LapackConnector::gemm('N', 'N', n_i, n_k, n_j, 1.0,
                                                  gf.c, n_j,
                                                  _Ncnt.c+i_nu*n_j*n_k, n_k,
                                                  1.0, MpMc.c+i_nu*n_i*n_k, n_k);
                    }
                }
                LapackConnector::gemm('T', 'T', n_mu, n_nu, n_i*n_k, 1.0,
                        Cs->c, n_mu,
                        MpMc.c, n_i*n_k, 1.0, chi0_tau[iR].c, n_nu);
            }
        }
    }

    // clean up N to free memory, as they have finished their job
    N.clear();
    Ncnt.clear();
    
    for ( auto iR: iRs )
    {
        matrix X(n_nu, n_j*n_i), Ynt(n_nu, n_i*n_j);
        matrix Xcnt(n_nu, n_j*n_i), Yc(n_nu, n_i*n_j);
        // X
        for ( auto const & iL_R2_Cs: LRI_Cs.at(Nu) )
        {
            auto const iL = iL_R2_Cs.first;
            auto const n_l = atom_nw[iL];
            for ( auto & R2_Cs: iL_R2_Cs.second )
            {
                auto const R2 = R2_Cs.first;
                auto const & Cs = R2_Cs.second;
                auto mR2mR = Vector3_Order<int>(-Rlist[iR]-R2) % R_period;
                auto i_mR2mR = get_R_index(Rlist, mR2mR);
                if ( gf_unocc_ab_t.count(i_mR2mR) &&
                     gf_unocc_ab_t.at(i_mR2mR).count(iL) &&
                     gf_unocc_ab_t.at(i_mR2mR).at(iL).count(iI) )
                {
                    const auto tran_Cs_nu_jl = transpose(*Cs);
                    auto const & gfu = gf_unocc_ab_t.at(i_mR2mR).at(iL).at(iI);
                    for ( int i_nu = 0; i_nu != n_nu; i_nu++ )
                        LapackConnector::gemm('N', 'N', n_j, n_i, n_l, 1.0,
                                              tran_Cs_nu_jl.c+i_nu*n_j*n_l, n_l,
                                              gfu.c, n_i,
                                              1.0, X.c+i_nu*n_j*n_i, n_i);
                }
                if ( gf_occ_ab_t.count(i_mR2mR) &&
                     gf_occ_ab_t.at(i_mR2mR).count(iL) &&
                     gf_occ_ab_t.at(i_mR2mR).at(iL).count(iI) )
                {
                    const auto tran_Cs_nu_jl = transpose(*Cs);
                    // gf shall take conjugate here
                    auto const & gf = gf_occ_ab_t.at(i_mR2mR).at(iL).at(iI);
                    for ( int i_nu = 0; i_nu != n_nu; i_nu++ )
                        LapackConnector::gemm('N', 'N', n_j, n_i, n_l, 1.0,
                                              tran_Cs_nu_jl.c+i_nu*n_j*n_l, n_l,
                                              gf.c, n_i,
                                              1.0, Xcnt.c+i_nu*n_j*n_i, n_i);
                }
            }
        }
        // Y
        for ( auto const & iK_R1_Cs: LRI_Cs.at(Mu) )
        {
            auto const iK= iK_R1_Cs.first;
            auto const n_k = atom_nw[iK];
            for ( auto & R1_Cs: iK_R1_Cs.second )
            {
                auto const R1 = R1_Cs.first;
                auto const Cs = R1_Cs.second;
                auto RmR1 = Vector3_Order<int>(Rlist[iR]-R1) % R_period;
                auto i_RmR1 = get_R_index(Rlist, RmR1);
                if ( gf_occ_ab_t.count(i_RmR1) &&
                     gf_occ_ab_t.at(i_RmR1).count(iK) &&
                     gf_occ_ab_t.at(i_RmR1).at(iK).count(iJ) )
                {
                    const matrix tran_Cs_mu_ik = transpose(*Cs);
                    auto const & gf = gf_occ_ab_t.at(i_RmR1).at(iK).at(iJ);
                    for ( int i_mu = 0; i_mu != n_mu; i_mu++ )
                        // transpose so that result stored in ji order
                        LapackConnector::gemm('T', 'T', n_j, n_i, n_k, 1.0,
                                              gf.c, n_j,
                                              tran_Cs_mu_ik.c+i_mu*n_i*n_k, n_k,
                                              1.0, Ynt.c+i_mu*n_j*n_i, n_i);
                }
                if ( gf_unocc_ab_t.count(i_RmR1) &&
                     gf_unocc_ab_t.at(i_RmR1).count(iK) &&
                     gf_unocc_ab_t.at(i_RmR1).at(iK).count(iJ) )
                {
                    const matrix tran_Cs_mu_ik = transpose(*Cs);
                    auto const & gfu = gf_unocc_ab_t.at(i_RmR1).at(iK).at(iJ);
                    for ( int i_mu = 0; i_mu != n_mu; i_mu++ )
                        // transpose so that result stored in ji order
                        LapackConnector::gemm('T', 'T', n_j, n_i, n_k, 1.0,
                                              gfu.c, n_j,
                                              tran_Cs_mu_ik.c+i_mu*n_i*n_k, n_k,
                                              1.0, Yc.c+i_mu*n_j*n_i, n_i);
                }
            }
        }
        LapackConnector::gemm('N', 'T', n_mu, n_nu, n_i*n_j, 1.0,
                Ynt.c, n_i*n_j,
                X.c, n_i*n_j, 1.0, chi0_tau[iR].c, n_nu);
        LapackConnector::gemm('N', 'T', n_mu, n_nu, n_i*n_j, 1.0,
                Yc.c, n_i*n_j,
                Xcnt.c, n_i*n_j, 1.0, chi0_tau[iR].c, n_nu);
    }
    /* if (LIBRPA::mpi_comm_global_h.myid==0 && Mu == 0 && Nu == 0) */
    if (mpi_comm_global_h.myid==0 && Mu == 0 && Nu == 1)
    {
        /* for (auto iR_chi0_tau: chi0_tau) */
        /* { */
            /* if (iR_chi0_tau.first!=4) continue; */
            /* cout << Rlist[iR_chi0_tau.first] << " Mu=" << Mu << " Nu=" << Nu; */
            /* print_matrix("chi tauR at first tau and R", iR_chi0_tau.second); */
        /* } */
    }
    /* cout << "Done real-space chi0" << endl; */

    return chi0_tau;
}

void Chi0::free_chi0_q(const double freq, const Vector3_Order<double> q)
{
    auto &chi0_for_free = chi0_q.at(freq).at(q);
    chi0_for_free.clear();
    map<size_t, map<size_t,ComplexMatrix>>().swap(chi0_for_free); 
}
