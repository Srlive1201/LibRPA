#include "epsilon.h"

#include <math.h>
#include <omp.h>

#include <algorithm>
#include <array>
#include <valarray>
#include <set>
#include <stdexcept>

#include "lapack_connector.h"
#include "envs_io.h"
#include "envs_mpi.h"
#include "utils_io.h"
#include "stl_io_helper.h"
#include "libri_utils.h"
#include "matrix_m_parallel_utils.h"
#include "parallel_mpi.h"
#include "params.h"
#include "pbc.h"
#include "atoms.h"
#include "profiler.h"
#include "scalapack_connector.h"
#include <omp.h>
#if defined(__MACH__)
#include <malloc/malloc.h> // for malloc_zone_pressure_relief and malloc_default_zone
#else
#include <malloc.h>
#endif

#ifdef LIBRPA_USE_LIBRI
#include <RI/comm/mix/Communicate_Tensors_Map_Judge.h>
#include <RI/global/Tensor.h>
using RI::Tensor;
using RI::Communicate_Tensors_Map_Judge::comm_map2_first;
#endif

using LIBRPA::envs::mpi_comm_global_h;
using LIBRPA::envs::blacs_ctxt_global_h;
using LIBRPA::Array_Desc;
using LIBRPA::envs::ofs_myid;
using LIBRPA::utils::lib_printf;

CorrEnergy compute_RPA_correlation_blacs_2d_gamma_only( Chi0 &chi0, atpair_k_cplx_mat_t &coulmat)
{
    CorrEnergy corr;
    if (mpi_comm_global_h.myid == 0)
        lib_printf("Calculating EcRPA with BLACS/ScaLAPACK 2D gamma_only\n");
    // lib_printf("Calculating EcRPA with BLACS, pid:  %d\n", mpi_comm_global_h.myid);
    const auto & mf = chi0.mf;
    const double CONE=1.0;
    const int n_abf = LIBRPA::atomic_basis_abf.nb_total;
    const auto part_range = LIBRPA::atomic_basis_abf.get_part_range();

    mpi_comm_global_h.barrier();

    Array_Desc desc_nabf_nabf(blacs_ctxt_global_h);
    // use a square blocksize instead max block, otherwise heev and inversion will complain about illegal parameter
    desc_nabf_nabf.init_square_blk(n_abf, n_abf, 0, 0);
    const auto set_IJ_nabf_nabf = LIBRPA::get_necessary_IJ_from_block_2D_sy('U', LIBRPA::atomic_basis_abf, desc_nabf_nabf);
    const auto s0_s1 = get_s0_s1_for_comm_map2_first(set_IJ_nabf_nabf);
    auto chi0_block = init_local_mat<double>(desc_nabf_nabf, MAJOR::COL);
    auto coul_block = init_local_mat<double>(desc_nabf_nabf, MAJOR::COL);
    auto coul_chi0_block = init_local_mat<double>(desc_nabf_nabf, MAJOR::COL);

    vector<Vector3_Order<double>> qpts;
    for (const auto &qMuNuchi: chi0.get_chi0_q().at(chi0.tfg.get_freq_nodes()[0]))
        qpts.push_back(qMuNuchi.first);

    complex<double> tot_RPA_energy(0.0, 0.0);
    map<Vector3_Order<double>, complex<double>> cRPA_q;
    if(mpi_comm_global_h.is_root())
        lib_printf("Finish init RPA blacs 2d\n");
#ifdef LIBRPA_USE_LIBRI
    for (const auto &q: qpts)
    {
        coul_block.zero_out();

        int iq = std::distance(klist.begin(), std::find(klist.begin(), klist.end(), q));
        std::array<double, 3> qa = {q.x, q.y, q.z};
        // collect the block elements of coulomb matrices
        {
            double vq_begin = omp_get_wtime();
            // LibRI tensor for communication, release once done
            std::map<int, std::map<std::pair<int, std::array<double, 3>>, Tensor<double>>> coul_libri;

            for (const auto &Mu_Nu: local_atpair)
            {
                const auto Mu = Mu_Nu.first;
                const auto Nu = Mu_Nu.second;
                ofs_myid << "myid " << blacs_ctxt_global_h.myid << "Mu " << Mu << " Nu " << Nu << endl;
                if (coulmat.count(Mu) == 0 ||
                    coulmat.at(Mu).count(Nu) == 0 ||
                    coulmat.at(Mu).at(Nu).count(q) == 0) continue;
                const auto &Vq = coulmat.at(Mu).at(Nu).at(q);
                const auto n_mu = LIBRPA::atomic_basis_abf.get_atom_nb(Mu);
                const auto n_nu = LIBRPA::atomic_basis_abf.get_atom_nb(Nu);
                matrix tmp_vq_real=(*Vq).real();
                std::valarray<double> Vq_va(tmp_vq_real.c, Vq->size);
                auto pvq = std::make_shared<std::valarray<double>>();
                *pvq = Vq_va;
                coul_libri[Mu][{Nu, std::array<double, 3>{0,0,0}}] = Tensor<double>({n_mu, n_nu}, pvq);
                coulmat.at(Mu).at(Nu).at(q).reset();
            }
            #ifndef __MACH__
            malloc_trim(0);
            #else
            malloc_zone_pressure_relief(malloc_default_zone(), 0);
            #endif
            //printf("Finish RPA blacs 2d  vq arr\n");
            double arr_end = omp_get_wtime();
            mpi_comm_global_h.barrier();
            double comm_begin = omp_get_wtime();
            //printf("Begin comm_map2_first  myid: %d\n",mpi_comm_global_h.myid);
            const auto IJq_coul = comm_map2_first(mpi_comm_global_h.comm, coul_libri, s0_s1.first, s0_s1.second);
            double comm_end = omp_get_wtime();
            mpi_comm_global_h.barrier();

            double block_begin = omp_get_wtime();

            collect_block_from_ALL_IJ_Tensor(coul_block, desc_nabf_nabf, LIBRPA::atomic_basis_abf,
                        qa,true, CONE, IJq_coul, MAJOR::ROW);
            double block_end = omp_get_wtime();
            lib_printf("Vq Time  myid: %d  arr_time: %f  comm_time: %f   block_time: %f   pair_size: %d\n",mpi_comm_global_h.myid,arr_end-vq_begin, comm_end-comm_begin, block_end-block_begin,set_IJ_nabf_nabf.size());
            mpi_comm_global_h.barrier();
            double vq_end = omp_get_wtime();

            if(mpi_comm_global_h.myid == 0)
                lib_printf(" | Total vq time: %f  lri_coul: %f   comm_vq: %f   block_vq: %f\n",vq_end-vq_begin, comm_begin-vq_begin,block_begin-comm_begin,vq_end-block_begin);
        }

        double chi_arr_time=0.0;
        double chi_comm_time=0.0;
        double chi_2d_time=0.0;
        for (const auto &freq: chi0.tfg.get_freq_nodes())
        {
            const auto ifreq = chi0.tfg.get_freq_index(freq);
            const double freq_weight = chi0.tfg.find_freq_weight(freq);
            double pi_freq_begin = omp_get_wtime();
            chi0_block.zero_out();
            {
                double chi_begin_arr = omp_get_wtime();
                std::map<int, std::map<std::pair<int, std::array<double, 3>>, Tensor<double>>> chi0_libri;
                const auto &chi0_wq = chi0.get_chi0_q().at(freq).at(q);

                for (const auto &M_Nchi: chi0_wq)
                {
                    const auto &M = M_Nchi.first;
                    const auto n_mu = LIBRPA::atomic_basis_abf.get_atom_nb(M);
                    for (const auto &N_chi: M_Nchi.second)
                    {
                        const auto &N = N_chi.first;
                        const auto n_nu = LIBRPA::atomic_basis_abf.get_atom_nb(N);
                        const auto &chi = N_chi.second.real();
                        std::valarray<double> chi_va(chi.c, chi.size);
                        auto pchi = std::make_shared<std::valarray<double>>();
                        *pchi = chi_va;
                        chi0_libri[M][{N, std::array<double, 3>{0,0,0}}] = Tensor<double>({n_mu, n_nu}, pchi);
                    }
                }

                // if(mpi_comm_global_h.is_root())
                // {
                //     lib_printf("Begin to clean chi0 !!! \n");
                //     system("free -m");
                //     lib_printf("chi0_freq_q size: %d\n",chi0_wq.size());
                // }

                chi0.free_chi0_q(freq,q);
                #ifndef __MACH__
                malloc_trim(0);
                #else
                malloc_zone_pressure_relief(malloc_default_zone(), 0);
                #endif

                // if(mpi_comm_global_h.is_root())
                // {
                //     lib_printf("After clean chi0 !!! \n");
                //     system("free -m");
                //     lib_printf("chi0_freq_q size: %d\n",chi0_wq.size());
                // }

                mpi_comm_global_h.barrier();
                double chi_end_arr = omp_get_wtime();
                // ofs_myid << "chi0_libri" << endl << chi0_libri;

                const auto IJq_chi0 = comm_map2_first(mpi_comm_global_h.comm, chi0_libri, s0_s1.first, s0_s1.second);
                // ofs_myid << "IJq_chi0" << endl << IJq_chi0;
                double chi_end_comm = omp_get_wtime();

                collect_block_from_ALL_IJ_Tensor(chi0_block, desc_nabf_nabf, LIBRPA::atomic_basis_abf,
                        qa,true, CONE, IJq_chi0, MAJOR::ROW);
                //printf("End collect block myid: %d ifreq: %d   TIME_USED: %f\n",mpi_comm_global_h.myid,ifreq,chi_end_comm-chi_end_arr);
                mpi_comm_global_h.barrier();
                double chi_end_2d = omp_get_wtime();

                chi_arr_time=(chi_end_arr-chi_begin_arr);
                chi_comm_time=(chi_end_comm-chi_end_arr);
                chi_2d_time=(chi_end_2d-chi_end_comm);
            }

            double pi_begin = omp_get_wtime();
            ScalapackConnector::pgemm_f('N', 'N', n_abf, n_abf, n_abf, 1.0,
                    coul_block.ptr(), 1, 1, desc_nabf_nabf.desc,
                    chi0_block.ptr(), 1, 1, desc_nabf_nabf.desc, 0.0,
                    coul_chi0_block.ptr(), 1, 1, desc_nabf_nabf.desc);
            // char fnp[100];
            // sprintf(fnp, "pi_ifreq_%d_iq_%d.mtx", ifreq, iq);
            double pi_end = omp_get_wtime();
            //printf("End pgemm  myid: %d ifreq: %d \n",mpi_comm_global_h.myid,ifreq);
            double trace_pi=0.0;
            double trace_pi_loc=0.0;
            for(int i=0;i!= n_abf;i++)
            {
                const int ilo = desc_nabf_nabf.indx_g2l_r(i);
                const int jlo = desc_nabf_nabf.indx_g2l_c(i);
                if (ilo >= 0 && jlo >= 0)
                    trace_pi_loc+=coul_chi0_block(ilo,jlo);
            }

            coul_chi0_block*=-1.0;
            for(int i=0;i!= n_abf;i++)
            {
                const int ilo = desc_nabf_nabf.indx_g2l_r(i);
                const int jlo = desc_nabf_nabf.indx_g2l_c(i);
                if (ilo >= 0 && jlo >= 0)
                    coul_chi0_block(ilo,jlo)+=CONE;
            }

            int *ipiv = new int [desc_nabf_nabf.m_loc()*10];
            int info;
            //printf("begin det  myid: %d ifreq: %d \n",mpi_comm_global_h.myid,ifreq);
            double ln_det=compute_pi_det_blacs_2d_gamma_only(coul_chi0_block, desc_nabf_nabf, ipiv, info);
            //printf("End det  myid: %d ifreq: %d \n",mpi_comm_global_h.myid,ifreq);
            double det_end = omp_get_wtime();
            mpi_comm_global_h.barrier();
            MPI_Allreduce(&trace_pi_loc,&trace_pi,1,MPI_DOUBLE,MPI_SUM,mpi_comm_global_h.comm);
            double pi_freq_end = omp_get_wtime();

            if(mpi_comm_global_h.myid==0)
            {
                lib_printf("| TIME of DET-freq-q:  %f,  q: ( %f, %f, %f)  TOT: %f  CHI_arr: %f  CHI_comm: %f, CHI_2d: %f, Pi: %f, Det: %f\n",freq, q.x,q.y,q.z,pi_freq_end-pi_freq_begin, chi_arr_time,chi_comm_time,chi_2d_time,pi_end-pi_begin,det_end-pi_end);
                complex<double> rpa_for_omega_q=complex<double>(trace_pi+ln_det);
                cRPA_q[q] += rpa_for_omega_q * freq_weight * irk_weight[q] / TWO_PI;//!check
                tot_RPA_energy += rpa_for_omega_q * freq_weight * irk_weight[q] / TWO_PI;
            }
        }
    }
#else
    throw std::logic_error("need compilation with LibRI");
#endif
    if(mpi_comm_global_h.myid==0)
    {
        for (auto &q_crpa : cRPA_q)
        {
            corr.qcontrib[q_crpa.first] = q_crpa.second;
        }
    }
    mpi_comm_global_h.barrier();
    corr.value = tot_RPA_energy;

    corr.etype = CorrEnergy::type::RPA;
    return corr;
}


CorrEnergy compute_RPA_correlation_blacs_2d( Chi0 &chi0,  atpair_k_cplx_mat_t &coulmat)
{
    lib_printf("Begin to compute_RPA_correlation_blacs_2d  myid: %d\n",mpi_comm_global_h.myid );
    system("free -m");
    CorrEnergy corr;
    if (mpi_comm_global_h.myid == 0)
        lib_printf("Calculating EcRPA with BLACS/ScaLAPACK 2D\n");
    // lib_printf("Calculating EcRPA with BLACS, pid:  %d\n", mpi_comm_global_h.myid);
    const auto & mf = chi0.mf;
    const complex<double> CONE{1.0, 0.0};
    const int n_abf = LIBRPA::atomic_basis_abf.nb_total;
    const auto part_range = LIBRPA::atomic_basis_abf.get_part_range();

    mpi_comm_global_h.barrier();

    Array_Desc desc_nabf_nabf(blacs_ctxt_global_h);
    // use a square blocksize instead max block, otherwise heev and inversion will complain about illegal parameter
    desc_nabf_nabf.init_square_blk(n_abf, n_abf, 0, 0);
    const auto set_IJ_nabf_nabf = LIBRPA::get_necessary_IJ_from_block_2D_sy('U', LIBRPA::atomic_basis_abf, desc_nabf_nabf);
    const auto s0_s1 = get_s0_s1_for_comm_map2_first(set_IJ_nabf_nabf);
    auto chi0_block = init_local_mat<complex<double>>(desc_nabf_nabf, MAJOR::COL);
    auto coul_block = init_local_mat<complex<double>>(desc_nabf_nabf, MAJOR::COL);
    auto coul_chi0_block = init_local_mat<complex<double>>(desc_nabf_nabf, MAJOR::COL);
    // ofs_myid << "Iset Jset " << s0_s1 << endl;
    // ofs_myid << "atpair_unordered_local of myid " << blacs_ctxt_global_h.myid << " " << atpair_unordered_local << endl;


    vector<Vector3_Order<double>> qpts;
    for (const auto &qMuNuchi: chi0.get_chi0_q().at(chi0.tfg.get_freq_nodes()[0]))
        qpts.push_back(qMuNuchi.first);

    complex<double> tot_RPA_energy(0.0, 0.0);
    map<Vector3_Order<double>, complex<double>> cRPA_q;
    if(mpi_comm_global_h.is_root())
        lib_printf("Finish init RPA blacs 2d\n");
#ifdef LIBRPA_USE_LIBRI
    for (const auto &q: qpts)
    {
        coul_block.zero_out();

        int iq = std::distance(klist.begin(), std::find(klist.begin(), klist.end(), q));
        std::array<double, 3> qa = {q.x, q.y, q.z};
        // collect the block elements of coulomb matrices
        {
            double vq_begin = omp_get_wtime();
            // LibRI tensor for communication, release once done
            std::map<int, std::map<std::pair<int, std::array<double, 3>>, Tensor<complex<double>>>> coul_libri;
            coul_libri.clear();
            for (const auto &Mu_Nu: local_atpair)
            {
                const auto Mu = Mu_Nu.first;
                const auto Nu = Mu_Nu.second;
                ofs_myid << "myid " << blacs_ctxt_global_h.myid << "Mu " << Mu << " Nu " << Nu << endl;
                if (coulmat.count(Mu) == 0 ||
                    coulmat.at(Mu).count(Nu) == 0 ||
                    coulmat.at(Mu).at(Nu).count(q) == 0) continue;
                const auto &Vq = coulmat.at(Mu).at(Nu).at(q);
                const auto n_mu = LIBRPA::atomic_basis_abf.get_atom_nb(Mu);
                const auto n_nu = LIBRPA::atomic_basis_abf.get_atom_nb(Nu);
                std::valarray<complex<double>> Vq_va(Vq->c, Vq->size);
                auto pvq = std::make_shared<std::valarray<complex<double>>>();
                *pvq = Vq_va;
                coul_libri[Mu][{Nu, qa}] = Tensor<complex<double>>({n_mu, n_nu}, pvq);
            }
            //printf("Finish RPA blacs 2d  vq arr\n");
            double arr_end = omp_get_wtime();
            mpi_comm_global_h.barrier();
            double comm_begin = omp_get_wtime();
            //printf("Begin comm_map2_first  myid: %d\n",mpi_comm_global_h.myid);
            const auto IJq_coul = comm_map2_first(mpi_comm_global_h.comm, coul_libri, s0_s1.first, s0_s1.second);
            double comm_end = omp_get_wtime();
            mpi_comm_global_h.barrier();
            //printf("End vq comm_map2_first  myid: %d   TIME_USED: %f\n",mpi_comm_global_h.myid,comm_end-comm_begin);
            // ofs_myid << "IJq_coul" << endl << IJq_coul;
            //printf("Finish RPA blacs 2d  vq 2d\n");
            double block_begin = omp_get_wtime();
            // for (const auto &IJ: set_IJ_nabf_nabf)
            // {
            //     const auto &I = IJ.first;
            //     const auto &J = IJ.second;
            //     // cout << IJq_coul.at(I).at({J, qa});
            //     collect_block_from_IJ_storage_syhe(
            //         coul_block, desc_nabf_nabf, LIBRPA::atomic_basis_abf, IJ.first,
            //         IJ.second, true, CONE, IJq_coul.at(I).at({J, qa}).ptr(), MAJOR::ROW);
            //     // lib_printf("myid %d I %d J %d nr %d nc %d\n%s",
            //     //        blacs_ctxt_global_h.myid, I, J,
            //     //        coul_block.nr(), coul_block.nc(),
            //     //        str(coul_block).c_str());
            // }
            collect_block_from_ALL_IJ_Tensor(coul_block, desc_nabf_nabf, LIBRPA::atomic_basis_abf,
                                             qa, true, CONE, IJq_coul, MAJOR::ROW);
            double block_end = omp_get_wtime();
            lib_printf("Vq Time  myid: %d  arr_time: %f  comm_time: %f   block_time: %f   pair_size: %d\n",mpi_comm_global_h.myid,arr_end-vq_begin, comm_end-comm_begin, block_end-block_begin,set_IJ_nabf_nabf.size());
            mpi_comm_global_h.barrier();
            double vq_end = omp_get_wtime();

            if(mpi_comm_global_h.myid == 0)
                lib_printf(" | Total vq time: %f  lri_coul: %f   comm_vq: %f   block_vq: %f\n",vq_end-vq_begin, comm_begin-vq_begin,block_begin-comm_begin,vq_end-block_begin);
        }


        //if(mpi_comm_global_h.is_root())
        //printf("Finish RPA blacs 2d  vq comm\n");
        // char fn[100];
        // sprintf(fn, "coul_iq_%d.mtx", iq);
        // print_matrix_mm_file_parallel(fn, coul_block, desc_nabf_nabf);
        // ofs_myid << str(coul_block);
        // lib_printf("coul_block\n%s", str(coul_block).c_str());
        double chi_arr_time=0.0;
        double chi_comm_time=0.0;
        double chi_2d_time=0.0;
        for (const auto &freq: chi0.tfg.get_freq_nodes())
        {
            const auto ifreq = chi0.tfg.get_freq_index(freq);
            const double freq_weight = chi0.tfg.find_freq_weight(freq);
            double pi_freq_begin = omp_get_wtime();
            chi0_block.zero_out();
            {
                double chi_begin_arr = omp_get_wtime();
                std::map<int, std::map<std::pair<int, std::array<double, 3>>, Tensor<complex<double>>>> chi0_libri;
                const auto &chi0_wq = chi0.get_chi0_q().at(freq).at(q);
                chi0_libri.clear();
                for (const auto &M_Nchi: chi0_wq)
                {
                    const auto &M = M_Nchi.first;
                    const auto n_mu = LIBRPA::atomic_basis_abf.get_atom_nb(M);
                    for (const auto &N_chi: M_Nchi.second)
                    {
                        const auto &N = N_chi.first;
                        const auto n_nu = LIBRPA::atomic_basis_abf.get_atom_nb(N);
                        const auto &chi = N_chi.second;
                        std::valarray<complex<double>> chi_va(chi.c, chi.size);
                        auto pchi = std::make_shared<std::valarray<complex<double>>>();
                        *pchi = chi_va;
                        chi0_libri[M][{N,qa}] = Tensor<complex<double>>({n_mu, n_nu}, pchi);
                    }
                }
                if(mpi_comm_global_h.is_root())
                {
                    lib_printf("Begin to clean chi0 !!! \n");
                    system("free -m");
                    lib_printf("chi0_freq_q size: %d,  freq: %f, q:( %f, %f, %f )\n",chi0_wq.size(),freq, q.x,q.y,q.z );
                }
                chi0.free_chi0_q(freq,q);
                #ifndef __MACH__
                malloc_trim(0);
                #else
                malloc_zone_pressure_relief(malloc_default_zone(), 0);
                #endif
                // if(mpi_comm_global_h.is_root())
                // {
                //     lib_printf("After clean chi0 !!! \n");
                //     system("free -m");
                //     lib_printf("chi0_freq_q size: %d\n",chi0_wq.size());
                // }
                mpi_comm_global_h.barrier();
                double chi_end_arr = omp_get_wtime();
                // ofs_myid << "chi0_libri" << endl << chi0_libri;

                const auto IJq_chi0 = comm_map2_first(mpi_comm_global_h.comm, chi0_libri, s0_s1.first, s0_s1.second);
                // ofs_myid << "IJq_chi0" << endl << IJq_chi0;
                double chi_end_comm = omp_get_wtime();
                collect_block_from_ALL_IJ_Tensor(chi0_block, desc_nabf_nabf, LIBRPA::atomic_basis_abf,
                                                 qa, true, CONE, IJq_chi0, MAJOR::ROW);
                mpi_comm_global_h.barrier();
                double chi_end_2d = omp_get_wtime();

                chi_arr_time=(chi_end_arr-chi_begin_arr);
                chi_comm_time=(chi_end_comm-chi_end_arr);
                chi_2d_time=(chi_end_2d-chi_end_comm);
                // char fnc[100];
                // sprintf(fnc, "chi_ifreq_%d_iq_%d.mtx", ifreq, iq);
                // if( ifreq== 0)
                //     print_matrix_mm_file_parallel(fnc, chi0_block, desc_nabf_nabf);
            }

            double pi_begin = omp_get_wtime();
            ScalapackConnector::pgemm_f('N', 'N', n_abf, n_abf, n_abf, 1.0,
                    coul_block.ptr(), 1, 1, desc_nabf_nabf.desc,
                    chi0_block.ptr(), 1, 1, desc_nabf_nabf.desc, 0.0,
                    coul_chi0_block.ptr(), 1, 1, desc_nabf_nabf.desc);
            // char fnp[100];
            // sprintf(fnp, "pi_ifreq_%d_iq_%d.mtx", ifreq, iq);
            double pi_end = omp_get_wtime();

            complex<double> trace_pi(0.0,0.0);
            complex<double> trace_pi_loc(0.0,0.0);
            for(int i=0;i!= n_abf;i++)
            {
                const int ilo = desc_nabf_nabf.indx_g2l_r(i);
                const int jlo = desc_nabf_nabf.indx_g2l_c(i);
                if (ilo >= 0 && jlo >= 0)
                    trace_pi_loc+=coul_chi0_block(ilo,jlo);
            }

            coul_chi0_block*=-1.0;
            for(int i=0;i!= n_abf;i++)
            {
                const int ilo = desc_nabf_nabf.indx_g2l_r(i);
                const int jlo = desc_nabf_nabf.indx_g2l_c(i);
                if (ilo >= 0 && jlo >= 0)
                    coul_chi0_block(ilo,jlo)+=CONE;
            }
            // if( ifreq== 0 && mpi_comm_global_h.is_root() )
            //     print_whole_matrix("pi-2D-loc", coul_chi0_block);

            int *ipiv = new int [desc_nabf_nabf.m_loc()*10];
            int info;
            complex<double> ln_det=compute_pi_det_blacs_2d(coul_chi0_block, desc_nabf_nabf, ipiv, info);
            double det_end = omp_get_wtime();
            mpi_comm_global_h.barrier();
            MPI_Allreduce(&trace_pi_loc,&trace_pi,1,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_global_h.comm);
            double pi_freq_end = omp_get_wtime();
            //double task_end = omp_get_wtime();
            // if(mpi_comm_global_h.is_root())
            //     lib_printf("| After det for freq:  %f,  q: ( %f, %f, %f)   TIME_LOCMAT: %f   TIME_DET: %f  TIME_CAL_Pi: %f, TIME_TRAN_LOC: %f\n",ifreq, q.x,q.y,q.z,task_mid-task_begin,task_end-task_mid,pi_time,loc_tran_time);
            //para_mpi.mpi_barrier();

            if(mpi_comm_global_h.myid==0)
            {
                lib_printf("| TIME of DET-freq-q:  %f,  q: ( %f, %f, %f)  TOT: %f  CHI_arr: %f  CHI_comm: %f, CHI_2d: %f, Pi: %f, Det: %f\n",freq, q.x,q.y,q.z,pi_freq_end-pi_freq_begin, chi_arr_time,chi_comm_time,chi_2d_time,pi_end-pi_begin,det_end-pi_end);
                complex<double> rpa_for_omega_q=trace_pi+ln_det;
                //cout << " ifreq:" << freq << "      rpa_for_omega_k: " << rpa_for_omega_q << "      lnt_det: " << ln_det << "    trace_pi " << trace_pi << endl;
                cRPA_q[q] += rpa_for_omega_q * freq_weight * irk_weight[q] / TWO_PI;//!check
                tot_RPA_energy += rpa_for_omega_q * freq_weight * irk_weight[q] / TWO_PI;
            }
        }
    }
#else
    throw std::logic_error("need compilation with LibRI");
#endif
    if(mpi_comm_global_h.myid==0)
    {
        for (auto &q_crpa : cRPA_q)
        {
            corr.qcontrib[q_crpa.first] = q_crpa.second;
            //cout << q_crpa.first << q_crpa.second << endl;
        }
        //cout << "gx_num_" << chi0.tfg.size() << "  tot_RPA_energy:  " << setprecision(8)    <<tot_RPA_energy << endl;
    }
    mpi_comm_global_h.barrier();
    corr.value = tot_RPA_energy;

    corr.etype = CorrEnergy::type::RPA;
    return corr;
}
double compute_pi_det_blacs_2d_gamma_only(matrix_m<double> &loc_piT, const Array_Desc &arrdesc_pi, int *ipiv, int &info)
{
    int one=1;
    int range_all= N_all_mu;
    int DESCPI_T[9];

    double det_begin = omp_get_wtime();

    ScalapackConnector::pgetrf_f(range_all,range_all,loc_piT.ptr(),one,one,arrdesc_pi.desc,ipiv, info);
    double trf_end = omp_get_wtime();

    double ln_det_loc=0.0;
    double ln_det_all=0.0;

    for(int ig=0;ig!=range_all;ig++)
    {

        int locr = arrdesc_pi.indx_g2l_r(ig);
        int locc = arrdesc_pi.indx_g2l_c(ig);
        if (locr >= 0 && locc >= 0)
        {

            double tmp_ln_det;
            if(loc_piT(locr,locc)>0)
            {
                tmp_ln_det=std::log(loc_piT(locr,locc));
            }
            else
            {
                tmp_ln_det=std::log(-loc_piT(locr,locc));
            }
            ln_det_loc+=tmp_ln_det;
		}
    }
    double ln_end = omp_get_wtime();

    MPI_Allreduce(&ln_det_loc,&ln_det_all,1,MPI_DOUBLE,MPI_SUM, mpi_comm_global_h.comm);
    double det_end = omp_get_wtime();
    return ln_det_all;
}
complex<double> compute_pi_det_blacs_2d(matrix_m<complex<double>> &loc_piT, const Array_Desc &arrdesc_pi, int *ipiv, int &info)
{
    int one=1;
    int range_all= N_all_mu;
    int DESCPI_T[9];
    // if(out_pi)
    // {
    //     print_complex_real_matrix("first_pi",pi_freq_q.at(0).at(0));
    //     print_complex_real_matrix("first_loc_piT_mat",loc_piT);
    // }
    double det_begin = omp_get_wtime();
    //ScalapackConnector::transpose_desc(DESCPI_T, arrdesc_pi.desc);
    pzgetrf_(&range_all,&range_all,loc_piT.ptr(),&one,&one,arrdesc_pi.desc,ipiv, &info);
    double trf_end = omp_get_wtime();
    //ScalapackConnector::pgetrf_f(range_all,range_all,loc_piT.c,one,one,DESCPI_T,ipiv, info);
    //printf("   after LU myid: %d\n",mpi_comm_global_h.myid);
    //printf("desc myid: %d,  m n: %d,%d,  mb nb: %d, %d,  loc_m_n: %d, %d, myp: %d,%d, npr,npc: %d, %d\n",mpi_comm_global_h.myid, arrdesc_pi.m(),arrdesc_pi.n(), arrdesc_pi.mb(),arrdesc_pi.nb(), arrdesc_pi.m_loc(),arrdesc_pi.n_loc(),arrdesc_pi.myprow(),arrdesc_pi.mypcol(),arrdesc_pi.nprows(),arrdesc_pi.npcols());
    complex<double> ln_det_loc(0.0,0.0);
    complex<double> ln_det_all(0.0,0.0);
    // complex<double> det_loc(1.0,0.0);
    // complex<double> det_glo(0.0,0.0);
    // vector<complex<double>>  det_dig;
    // vector<complex<double>>  ln_det_dig;
    // vector<complex<double>>  det_dig_r;
    // vector<complex<double>>  det_dig_c;
    //printf(" myid: %d ig=25, locr,locc: %d, %d)\n",mpi_comm_global_h.myid,arrdesc_pi.indx_g2l_r(25),arrdesc_pi.indx_g2l_c(25));
    for(int ig=0;ig!=range_all;ig++)
    {
        // int locr=para_mpi.localIndex(ig,row_nblk,para_mpi.nprow,para_mpi.myprow);
		// int locc=para_mpi.localIndex(ig,col_nblk,para_mpi.npcol,para_mpi.mypcol);
        int locr = arrdesc_pi.indx_g2l_r(ig);
        int locc = arrdesc_pi.indx_g2l_c(ig);
        if (locr >= 0 && locc >= 0)
        {
            // if(ipiv[locr]!=(ig+1))
            // 	det_loc=-1*det_loc * loc_piT(locc,locr);
            // else
            // 	det_loc=det_loc * loc_piT(locc,locr);
            // det_dig.push_back(loc_piT(locr,locc));
            // det_dig_r.push_back(locr);
            // det_dig_c.push_back(locc);
            complex<double> tmp_ln_det;
            if(loc_piT(locr,locc).real()>0)
            {
                tmp_ln_det=std::log(loc_piT(locr,locc));
                //ln_det_dig.push_back(tmp_ln_det);
            }
            else
            {
                tmp_ln_det=std::log(-loc_piT(locr,locc));
                //ln_det_dig.push_back(tmp_ln_det);
            }
            ln_det_loc+=tmp_ln_det;
		}
    }
    double ln_end = omp_get_wtime();
//     ComplexMatrix det_mm(loc_piT.nr(),loc_piT.nc());
//     for(int i=0;i!=loc_piT.nr();i++)
//         for(int j=0;j!=loc_piT.nc();j++)
//             det_mm(i,j)=loc_piT(i,j);
//    // sort(det_dig.rbegin(),det_dig.rend());
//     ComplexMatrix det_dig_mm(det_dig.size(),4);
//     for(int i=0;i!=det_dig.size();i++)
//     {
//         det_dig_mm(i,0) =det_dig_r[i];
//         det_dig_mm(i,1) =det_dig_c[i];
//         det_dig_mm(i,2)=det_dig[i];
//         det_dig_mm(i,3)=ln_det_dig[i];
//     }
//     char fn[100];
//     sprintf(fn, "det_dig_myid_%d.mtx", mpi_comm_global_h.myid);
//     print_complex_matrix_file("det_dig_loc", det_dig_mm, fn, false);

//     sprintf(fn, "det_mat_myid_%d.mtx", mpi_comm_global_h.myid);
//     print_complex_matrix_file("det_mat_loc", det_mm, fn, false);


    MPI_Allreduce(&ln_det_loc,&ln_det_all,1,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_global_h.comm);
    double det_end = omp_get_wtime();
    // if(mpi_comm_global_h.myid == 0)
    //     lib_printf("    | Det time   trf: %f   ln: %f   allreduce: %f\n",trf_end-det_begin,ln_end-trf_end, det_end-ln_end);
    //MPI_Allreduce(&det_loc,&det_glo,1,MPI_DOUBLE_COMPLEX,MPI_PROD,mpi_comm_global_h.comm);
    //ln_det_all=std::log(det_glo);
    return ln_det_all;
}

complex<double> compute_pi_det_blacs(ComplexMatrix &loc_piT, const Array_Desc &arrdesc_pi, int *ipiv, int &info)
{
    // int range_all = atom_mu_part_range[natom-1]+atom_mu[natom-1];
    // int desc_pi[9];
    // int loc_row, loc_col, info;
    // int row_nblk=1;
    // int col_nblk=1;
    int one=1;
    int range_all= N_all_mu;
    // para_mpi.set_blacs_mat(desc_pi,loc_row,loc_col,range_all,range_all,row_nblk,col_nblk);
    // int *ipiv = new int [loc_row*10];
    // ComplexMatrix loc_piT(loc_col,loc_row);

    // for(int i=0;i!=loc_row;i++)
    // {
    //     int global_row = para_mpi.globalIndex(i,row_nblk,para_mpi.nprow,para_mpi.myprow);
    //     int mu;
    //     int I=atom_mu_glo2loc(global_row,mu);
    //     for(int j=0;j!=loc_col;j++)
    //     {
    //         int global_col = para_mpi.globalIndex(j,col_nblk,para_mpi.npcol,para_mpi.mypcol);
    //         int nu;
    //         int J=atom_mu_glo2loc(global_col,nu);

    //         if( global_col == global_row)
    //         {
    //             loc_piT(j,i)=complex<double>(1.0,0.0) - pi_freq_q.at(I).at(J)(mu,nu);
    //         }
    //         else
    //         {
    //             loc_piT(j,i)=-1*  pi_freq_q.at(I).at(J)(mu,nu);
    //         }

    //     }
    // }
    int DESCPI_T[9];
    // if(out_pi)
    // {
    //     print_complex_real_matrix("first_pi",pi_freq_q.at(0).at(0));
    //     print_complex_real_matrix("first_loc_piT_mat",loc_piT);
    // }

    ScalapackConnector::transpose_desc(DESCPI_T, arrdesc_pi.desc);

   // para_mpi.mpi_barrier();
    //printf("   before LU Myid: %d        Available DOS memory = %ld bytes\n",mpi_comm_global_h.myid, memavail());
    //printf("   before LU myid: %d  range_all: %d,  loc_mat.size: %d\n",mpi_comm_global_h.myid,range_all,loc_piT.size);
    pzgetrf_(&range_all,&range_all,loc_piT.c,&one,&one,DESCPI_T,ipiv, &info);
    //printf("   after LU myid: %d\n",mpi_comm_global_h.myid);
    complex<double> ln_det_loc(0.0,0.0);
    complex<double> ln_det_all(0.0,0.0);
    for(int ig=0;ig!=range_all;ig++)
    {
        // int locr=para_mpi.localIndex(ig,row_nblk,para_mpi.nprow,para_mpi.myprow);
		// int locc=para_mpi.localIndex(ig,col_nblk,para_mpi.npcol,para_mpi.mypcol);
        int locr = arrdesc_pi.indx_g2l_r(ig);
        int locc = arrdesc_pi.indx_g2l_c(ig);
        if (locr >= 0 && locc >= 0)
        {
            // if(ipiv[locr]!=(ig+1))
            // 	det_loc=-1*det_loc * loc_piT(locc,locr);
            // else
            // 	det_loc=det_loc * loc_piT(locc,locr);
            if(loc_piT(locc,locr).real()>0)
                ln_det_loc+=std::log(loc_piT(locc,locr));
            else
                ln_det_loc+=std::log(-loc_piT(locc,locr));
		}
    }
    MPI_Allreduce(&ln_det_loc,&ln_det_all,1,MPI_DOUBLE_COMPLEX,MPI_SUM,mpi_comm_global_h.comm);
    return ln_det_all;
}


CorrEnergy compute_RPA_correlation_blacs(const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat)
{
    CorrEnergy corr;
    if (mpi_comm_global_h.myid == 0)
        lib_printf("Calculating EcRPA with BLACS/ScaLAPACK row\n");

    const auto & mf = chi0.mf;
    const complex<double> CONE{1.0, 0.0};
    const int n_abf = LIBRPA::atomic_basis_abf.nb_total;
    const auto part_range = LIBRPA::atomic_basis_abf.get_part_range();

    mpi_comm_global_h.barrier();

    LIBRPA::Array_Desc arrdesc_pi(blacs_ctxt_global_h);
    arrdesc_pi.init_square_blk(n_abf, n_abf, 0, 0);
    int loc_row = arrdesc_pi.m_loc(), loc_col = arrdesc_pi.n_loc(), info;

    // para_mpi.set_blacs_mat(desc_pi,loc_row,loc_col,N_all_mu,N_all_mu,row_nblk,col_nblk);
    int *ipiv = new int [loc_row*10];
    // double vq_begin_m2t= omp_get_wtime();
    // std::map<int, std::map<std::pair<int, std::array<double, 3>>, Tensor<complex<double>>>> vq_libri;
    // for(auto &Ip:Vq)
    // {
    //     auto I=Ip.first;
    //     const auto n_mu = LIBRPA::atomic_basis_abf.get_atom_nb(I);
    //     for(auto &Jp:Ip.second)
    //     {
    //         auto J=Jp.first;
    //         const auto n_nu = LIBRPA::atomic_basis_abf.get_atom_nb(J);
    //         for(auto &qp:Jp.second)
    //         {
    //             auto q=qp.first;
    //             std::array<double, 3> qa = {q.x, q.y, q.z};
    //             const auto &vq_ptr=qp.second;
    //             std::valarray<complex<double>> Vq_va(vq_ptr->c, vq_ptr->size);
    //             auto pvq = std::make_shared<std::valarray<complex<double>>>();
    //             *pvq = Vq_va;
    //             vq_libri[I][{J, qa}] = Tensor<complex<double>>({n_mu, n_nu}, pvq);
    //             if(I!=J)
    //             {
    //                 auto vqT=transpose(*vq_ptr, 1);
    //                 std::valarray<complex<double>> VqT_va(vqT.c, vqT.size);
    //                 auto pvqT = std::make_shared<std::valarray<complex<double>>>();
    //                 *pvqT = VqT_va;
    //                 vq_libri[J][{I, qa}] = Tensor<complex<double>>({n_nu, n_mu}, pvqT);
    //             }
    //         }
    //     }
    // }
    // double vq_end_m2t = omp_get_wtime();
    // set<int> loc_atp_IJ;
    // for(auto &atp:local_atpair)
    // {
    //     loc_atp_IJ.insert(atp.first);
    //     loc_atp_IJ.insert(atp.second);
    // }
    // set<int> all_atom_set;
    // for(int I=0;I!=natom;I++)
    //     all_atom_set.insert(I);
    // const auto IJq_coul = Communicate_Tensors_Map_Judge::comm_map2_first(mpi_comm_global_h.comm, vq_libri, all_atom_set, loc_atp_IJ);
    // atpair_k_cplx_mat_t Vq_loc;
    // double vq_end_comm = omp_get_wtime();
    // for(auto Ip:IJq_coul)
    // {
    //     auto I=Ip.first;
    //     auto n_mu=atom_mu[I];
    //     for(auto &Jqp:Ip.second)
    //     {
    //         auto J=Jqp.first.first;
    //         auto n_nu=atom_mu[J];
    //         auto qa=Jqp.first.second;
    //         Vector3_Order<double> q{qa[0],qa[1],qa[2]};
    //         shared_ptr<ComplexMatrix> vq_ptr = make_shared<ComplexMatrix>();
    //         vq_ptr->create(n_mu, n_nu);
    //         const auto length=sizeof(complex<double>)* n_mu *n_nu;
    //         memcpy((*vq_ptr).c, Jqp.second.ptr(),length);
    //         Vq_loc[I][J][q]=vq_ptr;
    //         //printf("| process %d, I: %d  J: %d\n",mpi_comm_global_h.myid, I,J );
    //     }
    // }
    // double vq_end_t2m = omp_get_wtime();
    // mpi_comm_global_h.barrier();
    // if(mpi_comm_global_h.is_root())
    //     lib_printf("| Vq_time %f, TIME_m2t: %f   TIME_comm: %f  TIME_t2m: %f\n",vq_end_t2m-vq_begin_m2t,vq_end_m2t-vq_begin_m2t,vq_end_comm-vq_end_m2t,vq_end_t2m-vq_end_comm);
    map<double, map<Vector3_Order<double>, ComplexMatrix>> pi_freq_q;
    complex<double> tot_RPA_energy(0.0, 0.0);
    map<Vector3_Order<double>, complex<double>> cRPA_q;
    for (const auto &freq_q_MuNuchi0 : chi0.get_chi0_q())
    {
        const auto freq = freq_q_MuNuchi0.first;
        const double freq_weight = chi0.tfg.find_freq_weight(freq);
        for (const auto &q_MuNuchi0 : freq_q_MuNuchi0.second)
        {
            double task_begin = omp_get_wtime();
            const auto q = q_MuNuchi0.first;
            auto &MuNuchi0 = q_MuNuchi0.second;

            //ComplexMatrix loc_piT(loc_col,loc_row);
            auto loc_piT = init_local_mat<complex<double>>(arrdesc_pi, MAJOR::COL);
            complex<double> trace_pi(0.0,0.0);
            double vq_time=0.0;
            double pi_time=0.0;
            double loc_tran_time=0.0;
            for(int Mu=0;Mu!=natom;Mu++)
            {
                double Mu_begin=omp_get_wtime();
               // lib_printf(" |process %d,  Mu:  %d\n",mpi_comm_global_h.myid,Mu);
                const size_t n_mu = atom_mu[Mu];
                atom_mapping<ComplexMatrix>::pair_t_old Vq_row=gather_vq_row_q(Mu,coulmat,q);
                double Mu_after_vq=omp_get_wtime();
                // atom_mapping<ComplexMatrix>::pair_t_old Vq_row;
                // const auto IJq_coul = Communicate_Tensors_Map_Judge::comm_map2_first(mpi_comm_global_h.comm, vq_libri, {Mu}, loc_atp_atoms);
                // double Mu_vq_comm = omp_get_wtime();
                // for(auto Ip:IJq_coul)
                // {
                //     auto I=Ip.first;
                //     auto n_mu=atom_mu[I];
                //     for(auto &Jqp:Ip.second)
                //     {
                //         auto J=Jqp.first.first;
                //         auto n_nu=atom_mu[J];
                //         auto q=Jqp.first.second;
                //         Vq_row[I][J].create(n_mu,n_nu);
                //         const auto length=sizeof(complex<double>)* n_mu *n_nu;
                //         memcpy(Vq_row[I][J].c, Jqp.second.ptr(),length);
                //     }
                // }
                // double Mu_after_vq=omp_get_wtime();
                //printf("   |process %d, Mu: %d  vq_row.size: %d\n",para_mpi.get_myid(),Mu,Vq_row[Mu].size());
                //ComplexMatrix loc_pi_row=compute_Pi_freq_q_row(q,MuNuchi0,Vq_loc,Mu,q);
                ComplexMatrix loc_pi_row=compute_Pi_freq_q_row(q,MuNuchi0,Vq_row,Mu);
                //printf("   |process %d,   compute_pi\n",para_mpi.get_myid());
                ComplexMatrix glo_pi_row( n_mu,N_all_mu);
                mpi_comm_global_h.barrier();
                mpi_comm_global_h.allreduce_ComplexMatrix(loc_pi_row,glo_pi_row);
                double Mu_after_pi_loc=omp_get_wtime();
                //cout<<"  glo_pi_rowT nr,nc: "<<glo_pi_row.nr<<" "<<glo_pi_row.nc<<endl;

                for(int i_mu=0;i_mu!=n_mu;i_mu++)
                    trace_pi+=glo_pi_row(i_mu,atom_mu_part_range[Mu]+i_mu);
                //select glo_pi_rowT to pi_blacs
                for(int i=0;i!=loc_row;i++)
                {
                    // int global_row = para_mpi.globalIndex(i,row_nblk,para_mpi.nprow,para_mpi.myprow);
                    int global_row = arrdesc_pi.indx_l2g_r(i);
                    int mu_blacs;
                    int I_blacs=atom_mu_glo2loc(global_row,mu_blacs);
                    if(I_blacs== Mu)
                        for(int j=0;j!=loc_col;j++)
                        {
                            // int global_col = para_mpi.globalIndex(j,col_nblk,para_mpi.npcol,para_mpi.mypcol);
                            int global_col = arrdesc_pi.indx_l2g_c(j);
                            int nu_blacs;
                            int J_blacs=atom_mu_glo2loc(global_col,nu_blacs);
                            //cout<<" Mu: "<<Mu<<"  i,j: "<<i<<"  "<<j<<"    glo_row,col: "<<global_row<<"  "<<global_col<<"  J:"<<J_blacs<< "  index i,j: "<<atom_mu_part_range[J_blacs] + mu_blacs<<" "<<nu_blacs<<endl;
                            if( global_col == global_row)
                            {
                                loc_piT(i,j)=complex<double>(1.0,0.0) - glo_pi_row(mu_blacs, atom_mu_part_range[J_blacs]+nu_blacs);
                            }
                            else
                            {
                                loc_piT(i,j) = -glo_pi_row(mu_blacs, atom_mu_part_range[J_blacs]+nu_blacs);
                            }

                        }

                }
                double Mu_after_loc_tran=omp_get_wtime();
                vq_time+=(Mu_after_vq-Mu_begin);
                pi_time+=(Mu_after_pi_loc-Mu_after_vq);
                loc_tran_time+=(Mu_after_loc_tran-Mu_after_pi_loc);

            }
            // if(freq == chi0.tfg.get_freq_nodes()[0] && mpi_comm_global_h.is_root())
            //     print_complex_matrix(" loc_piT",loc_piT);
            double task_mid = omp_get_wtime();
            //printf("|process  %d, before det\n",mpi_comm_global_h.myid);
            complex<double> ln_det=compute_pi_det_blacs_2d(loc_piT, arrdesc_pi, ipiv, info);
            double task_end = omp_get_wtime();
            if(mpi_comm_global_h.is_root())
                lib_printf("| After det for freq:  %f,  q: ( %f, %f, %f)   TIME_Vq_COMM: %f   TIME_DET: %f  TIME_CAL_Pi: %f, TIME_TRAN_LOC: %f\n",freq, q.x,q.y,q.z,vq_time,task_end-task_mid,pi_time,loc_tran_time);
            //para_mpi.mpi_barrier();
            if(mpi_comm_global_h.myid==0)
            {
                complex<double> rpa_for_omega_q=trace_pi+ln_det;
                //cout << " ifreq:" << freq << "      rpa_for_omega_k: " << rpa_for_omega_q << "      lnt_det: " << ln_det << "    trace_pi " << trace_pi << endl;
                cRPA_q[q] += rpa_for_omega_q * freq_weight * irk_weight[q] / TWO_PI;//!check
                tot_RPA_energy += rpa_for_omega_q * freq_weight * irk_weight[q] / TWO_PI;
            }
        }
    }

    if(mpi_comm_global_h.myid==0)
    {
        for (auto &q_crpa : cRPA_q)
        {
            corr.qcontrib[q_crpa.first] = q_crpa.second;
            //cout << q_crpa.first << q_crpa.second << endl;
        }
        //cout << "gx_num_" << chi0.tfg.size() << "  tot_RPA_energy:  " << setprecision(8)    <<tot_RPA_energy << endl;
    }
    mpi_comm_global_h.barrier();
    corr.value = tot_RPA_energy;
    corr.etype = CorrEnergy::type::RPA;
    return corr;
}

CorrEnergy compute_RPA_correlation(const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat)
{
    CorrEnergy corr;
    if (mpi_comm_global_h.myid == 0)
        lib_printf("Calculating EcRPA without BLACS/ScaLAPACK\n");
    // lib_printf("Begin cal cRPA , pid:  %d\n", mpi_comm_global_h.myid);
    const auto & mf = chi0.mf;

    // freq, q
    map<double, map<Vector3_Order<double>, atom_mapping<ComplexMatrix>::pair_t_old>> pi_freq_q_Mu_Nu;
    if (LIBRPA::parallel_routing == LIBRPA::ParallelRouting::ATOM_PAIR || LIBRPA::parallel_routing == LIBRPA::ParallelRouting::LIBRI)
        pi_freq_q_Mu_Nu = compute_Pi_q_MPI(chi0, coulmat);
    else
        pi_freq_q_Mu_Nu = compute_Pi_q(chi0, coulmat);
    lib_printf("Finish Pi freq on Proc %4d, size %zu\n", mpi_comm_global_h.myid, pi_freq_q_Mu_Nu.size());
    //mpi_comm_global_h.barrier();

    int range_all = N_all_mu;


    vector<int> part_range;
    part_range.resize(atom_mu.size());
    part_range[0] = 0;
    int count_range = 0;
    for (int I = 0; I != atom_mu.size() - 1; I++)
    {
        count_range += atom_mu[I];
        part_range[I + 1] = count_range;
    }

    // cout << "part_range:" << endl;
    // for (int I = 0; I != atom_mu.size(); I++)
    // {
    //     cout << part_range[I] << endl;
    // }
    // cout << "part_range over" << endl;

    // pi_freq_q contains all atoms
    map<double, map<Vector3_Order<double>, ComplexMatrix>> pi_freq_q;

    for (const auto &freq_q_MuNupi : pi_freq_q_Mu_Nu)
    {
        const auto freq = freq_q_MuNupi.first;

        for (const auto &q_MuNupi : freq_q_MuNupi.second)
        {
            const auto q = q_MuNupi.first;
            const auto MuNupi = q_MuNupi.second;
            pi_freq_q[freq][q].create(range_all, range_all);

            ComplexMatrix pi_munu_tmp(range_all, range_all);
            pi_munu_tmp.zero_out();

            for (const auto &Mu_Nupi : MuNupi)
            {
                const auto Mu = Mu_Nupi.first;
                const auto Nupi = Mu_Nupi.second;
                const size_t n_mu = atom_mu[Mu];
                for (const auto &Nu_pi : Nupi)
                {
                    const auto Nu = Nu_pi.first;
                    const auto pimat = Nu_pi.second;
                    const size_t n_nu = atom_mu[Nu];

                    for (size_t mu = 0; mu != n_mu; ++mu)
                    {
                        for (size_t nu = 0; nu != n_nu; ++nu)
                        {
                            pi_munu_tmp(part_range[Mu] + mu, part_range[Nu] + nu) += pimat(mu, nu);
                        }
                    }
                }
            }
            if (LIBRPA::parallel_routing == LIBRPA::ParallelRouting::ATOM_PAIR || LIBRPA::parallel_routing == LIBRPA::ParallelRouting::LIBRI)
            {
                mpi_comm_global_h.reduce_ComplexMatrix(pi_munu_tmp, pi_freq_q.at(freq).at(q), 0);
            }
            else
            {
                pi_freq_q.at(freq).at(q) = std::move(pi_munu_tmp);
            }
        }
    }
    // lib_printf("Finish Pi communicate %4d, size %zu\n", mpi_comm_global_h.myid, pi_freq_q_Mu_Nu.size());
    mpi_comm_global_h.barrier();
    // if (mpi_comm_global_h.myid == 0)
    {
        complex<double> tot_RPA_energy(0.0, 0.0);
        map<Vector3_Order<double>, complex<double>> cRPA_q;
        for (const auto &freq_qpi : pi_freq_q)
        {
            const auto freq = freq_qpi.first;
            const double freq_weight = chi0.tfg.find_freq_weight(freq);
            for (const auto &q_pi : freq_qpi.second)
            {
                const auto q = q_pi.first;
                const auto pimat = q_pi.second;
                complex<double> rpa_for_omega_q(0.0, 0.0);
                ComplexMatrix identity(range_all, range_all);
                ComplexMatrix identity_minus_pi(range_all, range_all);
                identity.set_as_identity_matrix();
                identity_minus_pi = identity - pi_freq_q[freq][q];
                complex<double> det_for_rpa(1.0, 0.0);
                int info_LU = 0;
                int *ipiv = new int[range_all];
                LapackConnector::zgetrf(range_all, range_all, identity_minus_pi, range_all, ipiv, &info_LU);
                for (int ib = 0; ib != range_all; ib++)
                {
                    if (ipiv[ib] != (ib + 1))
                        det_for_rpa = -det_for_rpa * identity_minus_pi(ib, ib);
                    else
                        det_for_rpa = det_for_rpa * identity_minus_pi(ib, ib);
                }
                delete[] ipiv;

                complex<double> trace_pi;
                complex<double> ln_det;
                ln_det = std::log(det_for_rpa);
                trace_pi = trace(pi_freq_q.at(freq).at(q));
                // cout << "PI trace vector:" << endl;
                // cout << endl;
                rpa_for_omega_q = ln_det + trace_pi;
                // cout << " ifreq:" << freq << "      rpa_for_omega_k: " << rpa_for_omega_q << "      lnt_det: " << ln_det << "    trace_pi " << trace_pi << endl;
                cRPA_q[q] += rpa_for_omega_q * freq_weight * irk_weight[q] / TWO_PI;
                tot_RPA_energy += rpa_for_omega_q * freq_weight * irk_weight[q]  / TWO_PI;
            }
        }
        // lib_printf("Finish EcRPA %4d, size %zu\n", mpi_comm_global_h.myid, pi_freq_q_Mu_Nu.size());
        mpi_comm_global_h.barrier();
        map<Vector3_Order<double>, complex<double>> global_cRPA_q;
        for ( auto q_weight: irk_weight)
        {
            MPI_Reduce(&cRPA_q[q_weight.first],&global_cRPA_q[q_weight.first],1,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_global_h.comm);
        }

        for (auto &q_crpa : global_cRPA_q)
        {
            corr.qcontrib[q_crpa.first] = q_crpa.second;
        }
        complex<double> gather_tot_RPA_energy(0.0, 0.0);
        MPI_Reduce(&tot_RPA_energy,&gather_tot_RPA_energy,1,MPI_DOUBLE_COMPLEX,MPI_SUM,0,mpi_comm_global_h.comm);
        corr.value = gather_tot_RPA_energy;
    }
    corr.etype = CorrEnergy::type::RPA;
    return corr;
}

CorrEnergy compute_MP2_correlation(const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat)
{
    CorrEnergy corr;
    corr.etype = CorrEnergy::type::MP2;
    return corr;
}

map<double, map<Vector3_Order<double>, atom_mapping<ComplexMatrix>::pair_t_old>> compute_Pi_q(const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat)
{
    map<double, map<Vector3_Order<double>, atom_mapping<ComplexMatrix>::pair_t_old>> pi;
    lib_printf("Begin compute_Pi_q , pid:  %d\n", mpi_comm_global_h.myid);
    for (auto const & freq_qJQchi0: chi0.get_chi0_q())
    {
        const double freq = freq_qJQchi0.first;
        for (auto &q_JQchi0 : freq_qJQchi0.second)
        {
            Vector3_Order<double> q = q_JQchi0.first;
            for (auto &JQchi0 : q_JQchi0.second )
            {
                const size_t J = JQchi0.first;
                const size_t J_mu = atom_mu[J];
                for (auto &Qchi0 : JQchi0.second)
                {
                    const size_t Q = Qchi0.first;
                    const size_t Q_mu = atom_mu[Q];
                    // auto &chi0_mat = Qchi0.second;
                    for (int I=0;I!=natom;I++)
                    {
                        //const size_t I = I_p.first;
                        const size_t I_mu = atom_mu[I];
                        pi[freq][q][I][Q].create(I_mu, Q_mu);
                        if (J != Q)
                            pi[freq][q][I][J].create(I_mu, J_mu);
                    }
                }
            }
            // if(freq==chi0.tfg.get_freq_nodes()[0])
            //     for(auto &Ip:pi[freq][q])
            //         for(auto &Jp:Ip.second)
            //             lib_printf("  |process  %d, pi atpair: %d, %d \n",mpi_comm_global_h.myid,Ip.first,Jp.first);
        }

    }

    // ofstream fp;
    // std::stringstream ss;
    // ss<<"out_pi_rank_"<<mpi_comm_global_h.myid<<".txt";
    // fp.open(ss.str());
    for (auto &freq_p : chi0.get_chi0_q())
    {
        const double freq = freq_p.first;
        const auto chi0_freq = freq_p.second;
        for (auto &k_pair : chi0_freq)
        {
            Vector3_Order<double> ik_vec = k_pair.first;
            auto chi0_freq_k = k_pair.second;
            for (auto &J_p : chi0_freq_k)
            {
                const size_t J = J_p.first;
                for (auto &Q_p : J_p.second)
                {
                    const size_t Q = Q_p.first;
                    auto &chi0_mat = Q_p.second;
                    for (int I=0;I!=natom;I++)
                    {
                        //const size_t I = I_p.first;
                        //printf("cal_pi  pid: %d , IJQ:  %d  %d  %d\n", mpi_comm_global_h.myid, I, J, Q);
                        //  cout<<"         pi_IQ: "<<pi_k.at(freq).at(ik_vec).at(I).at(Q)(0,0)<<"   pi_IJ: "<<pi_k.at(freq).at(ik_vec).at(I).at(J)(0,0);
                        if (I <= J)
                        {
                            // if (freq == chi0.tfg.get_freq_nodes()[0])
                            //     lib_printf("cal_pi  pid: %d , IJQ:  %d  %d  %d   type: %d \n", mpi_comm_global_h.myid, I, J, Q,1);
                            //      << "  Vq: " << (*Vq.at(I).at(J).at(ik_vec))(0, 0) << endl;
                            pi.at(freq).at(ik_vec).at(I).at(Q) += (*Vq.at(I).at(J).at(ik_vec)) * chi0_mat;
                            // if (freq == chi0.tfg.get_freq_nodes()[0])
                            // {
                            //     std:stringstream sm;
                            //     complex<double> trace_pi;
                            //     trace_pi = trace(pi.at(freq).at(ik_vec).at(I).at(Q));
                            //     sm << " IJQ: " << I << " " << J << " " << Q << "  ik_vec: " << ik_vec << "  trace_pi:  " << trace_pi << endl;
                            //     print_complex_matrix_file(sm.str().c_str(), (*Vq.at(I).at(J).at(ik_vec)),fp,false);
                            //     print_complex_matrix_file("chi0:", chi0_mat,fp,false);
                            //     print_complex_matrix_file("pi_mat:", pi.at(freq).at(ik_vec).at(I).at(Q),fp,false);
                            // }
                        }
                        else
                        {
                            // if (freq == chi0.tfg.get_freq_nodes()[0])
                            //     lib_printf("cal_pi  pid: %d , IJQ:  %d  %d  %d   type: %d \n", mpi_comm_global_h.myid, I, J, Q,2);
                            //      << "  Vq: " << transpose(*Vq.at(J).at(I).at(ik_vec), 1)(0, 0) << endl;
                            pi.at(freq).at(ik_vec).at(I).at(Q) += transpose(*Vq.at(J).at(I).at(ik_vec), 1) * chi0_mat;
                        }

                        if (J != Q)
                        {
                            ComplexMatrix chi0_QJ = transpose(chi0_mat, 1);
                            if (I <= Q)
                            {
                                // if (freq == chi0.tfg.get_freq_nodes()[0])
                                //     lib_printf("cal_pi  pid: %d , IJQ:  %d  %d  %d   type: %d \n", mpi_comm_global_h.myid, I, J, Q,3);
                                //      << "  Vq: " << (*Vq.at(I).at(Q).at(ik_vec))(0, 0) << endl;
                                pi.at(freq).at(ik_vec).at(I).at(J) += (*Vq.at(I).at(Q).at(ik_vec)) * chi0_QJ;
                            }
                            else
                            {
                                // if (freq == chi0.tfg.get_freq_nodes()[0])
                                //     lib_printf("cal_pi  pid: %d , IJQ:  %d  %d  %d   type: %d \n", mpi_comm_global_h.myid, I, J, Q,4);
                                //      << "  Vq: " << transpose(*Vq.at(J).at(I).at(ik_vec), 1)(0, 0) << endl;
                                pi.at(freq).at(ik_vec).at(I).at(J) += transpose(*Vq.at(Q).at(I).at(ik_vec), 1) * chi0_QJ;
                            }
                        }
                    }
                }
            }
        }
    }
    // fp.close();
    //print_complex_matrix(" first_pi_mat:",pi.at(chi0.tfg.get_freq_nodes()[0]).at({0,0,0}).at(0).at(0));
    /* print_complex_matrix("  last_pi_mat:",pi.at(chi0.tfg.get_freq_nodes()[0]).at({0,0,0}).at(natom-1).at(natom-1)); */
    return pi;
}

map<double, map<Vector3_Order<double>, atom_mapping<ComplexMatrix>::pair_t_old>> compute_Pi_q_MPI(const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat)
{
    map<double, map<Vector3_Order<double>, atom_mapping<ComplexMatrix>::pair_t_old>> pi;
    lib_printf("Begin compute_Pi_q_MPI , pid:  %d\n", mpi_comm_global_h.myid);
    for (auto const & freq_qJQchi0: chi0.get_chi0_q())
    {
        const double freq = freq_qJQchi0.first;
        for (auto &q_JQchi0 : freq_qJQchi0.second)
        {
            Vector3_Order<double> q = q_JQchi0.first;
            for (auto &JQchi0 : q_JQchi0.second )
            {
                const size_t J = JQchi0.first;
                const size_t J_mu = atom_mu[J];
                for (auto &Qchi0 : JQchi0.second)
                {
                    const size_t Q = Qchi0.first;
                    const size_t Q_mu = atom_mu[Q];
                    // auto &chi0_mat = Qchi0.second;
                    for (int I=0;I!=natom;I++)
                    {
                        //const size_t I = I_p.first;
                        const size_t I_mu = atom_mu[I];
                        pi[freq][q][I][Q].create(I_mu, Q_mu);
                        if (J != Q)
                            pi[freq][q][I][J].create(I_mu, J_mu);
                    }
                }
            }
            // if(freq==chi0.tfg.get_freq_nodes()[0])
            //     for(auto &Ip:pi[freq][q])
            //         for(auto &Jp:Ip.second)
            //             lib_printf("  |process  %d, pi atpair: %d, %d \n",mpi_comm_global_h.myid,Ip.first,Jp.first);
        }

    }

    // ofstream fp;
    // std::stringstream ss;
    // ss<<"out_pi_rank_"<<mpi_comm_global_h.myid<<".txt";
    // fp.open(ss.str());
    for (auto &k_pair : irk_weight)
    {
        Vector3_Order<double> ik_vec = k_pair.first;
        for (int I=0;I!=natom;I++)
        {
            atom_mapping<ComplexMatrix>::pair_t_old Vq_row=gather_vq_row_q(I,coulmat,ik_vec);
            for (auto &freq_p : chi0.get_chi0_q())
            {
                const double freq = freq_p.first;
                const auto chi0_freq = freq_p.second;


                auto chi0_freq_k = freq_p.second.at(ik_vec);

                for (auto &J_p : chi0_freq_k)
                {
                    const size_t J = J_p.first;
                    for (auto &Q_p : J_p.second)
                    {
                        const size_t Q = Q_p.first;
                        auto &chi0_mat = Q_p.second;

                        //const size_t I = I_p.first;
                        //printf("cal_pi  pid: %d , IJQ:  %d  %d  %d\n", mpi_comm_global_h.myid, I, J, Q);
                        //  cout<<"         pi_IQ: "<<pi_k.at(freq).at(ik_vec).at(I).at(Q)(0,0)<<"   pi_IJ: "<<pi_k.at(freq).at(ik_vec).at(I).at(J)(0,0);

                        // if (freq == chi0.tfg.get_freq_nodes()[0])
                        //     lib_printf("cal_pi  pid: %d , IJQ:  %d  %d  %d   type: %d \n", mpi_comm_global_h.myid, I, J, Q,1);
                        //      << "  Vq: " << (*Vq.at(I).at(J).at(ik_vec))(0, 0) << endl;
                        pi.at(freq).at(ik_vec).at(I).at(Q) += Vq_row.at(I).at(J) * chi0_mat;
                        // if (freq == chi0.tfg.get_freq_nodes()[0])
                        // {
                        //     std:stringstream sm;
                        //     complex<double> trace_pi;
                        //     trace_pi = trace(pi.at(freq).at(ik_vec).at(I).at(Q));
                        //     sm << " IJQ: " << I << " " << J << " " << Q << "  ik_vec: " << ik_vec << "  trace_pi:  " << trace_pi << endl;
                        //     print_complex_matrix_file(sm.str().c_str(), Vq_row.at(I).at(J),fp,false);
                        //     print_complex_matrix_file("chi0:", chi0_mat,fp,false);
                        //     print_complex_matrix_file("pi_mat:", pi.at(freq).at(ik_vec).at(I).at(Q),fp,false);
                        // }

                        if (J != Q)
                        {
                            ComplexMatrix chi0_QJ = transpose(chi0_mat, 1);
                            // if (freq == chi0.tfg.get_freq_nodes()[0])
                            //     lib_printf("cal_pi  pid: %d , IJQ:  %d  %d  %d   type: %d \n", mpi_comm_global_h.myid, I, J,Q,3);
                            //      << "  Vq: " << (*Vq.at(I).at(Q).at(ik_vec))(0, 0) << endl;
                            pi.at(freq).at(ik_vec).at(I).at(J) += Vq_row.at(I).at(Q) * chi0_QJ;
                        }
                    }
                }
            }
        }
    }
    //fp.close();
    // print_complex_matrix(" first_pi_mat:",pi.at(chi0.tfg.get_freq_nodes()[0]).at({0,0,0}).at(0).at(0));
    /* print_complex_matrix("  last_pi_mat:",pi.at(chi0.tfg.get_freq_nodes()[0]).at({0,0,0}).at(natom-1).at(natom-1)); */
    lib_printf("End compute_Pi_q_MPI , pid:  %d\n", mpi_comm_global_h.myid);
    return pi;
}

ComplexMatrix compute_Pi_freq_q_row(const Vector3_Order<double> &ik_vec, const atom_mapping<ComplexMatrix>::pair_t_old &chi0_freq_q, const atom_mapping<ComplexMatrix>::pair_t_old &Vq_row, const int &I)
{
    map<size_t,ComplexMatrix> pi;
    // lib_printf("Begin cal_pi_k , pid:  %d\n", para_mpi.get_myid());
    auto I_mu=atom_mu[I];
    for(int J=0;J!=natom;J++)
        pi[J].create(I_mu,atom_mu[J]);

omp_lock_t pi_lock;
omp_init_lock(&pi_lock);
#pragma omp parallel for schedule(dynamic)
    for (int iap=0;iap!=local_atpair.size();iap++)
    {
        const size_t J = local_atpair[iap].first;
        const size_t Q = local_atpair[iap].second;
        auto &chi0_mat= chi0_freq_q.at(J).at(Q);
        auto tmp_pi_mat= Vq_row.at(I).at(J) * chi0_mat;
        ComplexMatrix chi0_QJ = transpose(chi0_mat, 1);
        auto tmp_pi_mat2= Vq_row.at(I).at(Q) * chi0_QJ;
        omp_set_lock(&pi_lock);
        pi.at(Q)+=tmp_pi_mat;
        if(J!=Q)
            {
                pi.at(J)+=tmp_pi_mat2;
            }
        omp_unset_lock(&pi_lock);
    }
    omp_destroy_lock(&pi_lock);
    // for (auto &J_p : chi0_freq_q)
    // {
    //     const size_t J = J_p.first;
    //     for (auto &Q_p : J_p.second)
    //     {
    //         const size_t Q = Q_p.first;
    //         auto &chi0_mat = Q_p.second;
    //         pi.at(Q) += Vq_row.at(I).at(J) * chi0_mat;
    //         if (J != Q)
    //         {
    //             ComplexMatrix chi0_QJ = transpose(chi0_mat, 1);
    //             pi.at(J) += Vq_row.at(I).at(Q) * chi0_QJ;
    //         }
    //     }
    // }
    //Pi_rowT
    // ComplexMatrix pi_row(N_all_mu,atom_mu[I]);
    // complex<double> *pi_row_ptr=pi_row.c;
    // for(auto &Jp:pi)
    // {
    //     auto J=Jp.first;
    //     auto J_mu=atom_mu[J];
    //     const auto length=sizeof(complex<double>)* I_mu *J_mu;
    //     memcpy(pi_row_ptr, pi.at(J).c,length);
    //     pi_row_ptr+=I_mu *J_mu;
    // }
    ComplexMatrix pi_row(atom_mu[I],N_all_mu);
    for(int i=0;i!=pi_row.nr;i++)
        for(int J=0;J!=natom;J++)
            for(int j=0;j!=atom_mu[J];j++)
                pi_row(i,atom_mu_part_range[J]+j)=pi.at(J)(i,j);
    return pi_row;
}

ComplexMatrix compute_Pi_freq_q_row_ri(const Vector3_Order<double> &ik_vec, const atom_mapping<ComplexMatrix>::pair_t_old &chi0_freq_q, const atpair_k_cplx_mat_t &Vq_loc, const int &I, const Vector3_Order<double> &q)
{
    map<size_t,ComplexMatrix> pi;
    // lib_printf("Begin cal_pi_k , pid:  %d\n", mpi_comm_global_h.myid);
    auto I_mu=atom_mu[I];
    for(int J=0;J!=natom;J++)
        pi[J].create(I_mu,atom_mu[J]);

omp_lock_t pi_lock;
omp_init_lock(&pi_lock);
#pragma omp parallel for schedule(dynamic)
    for (int iap=0;iap!=local_atpair.size();iap++)
    {
        const size_t J = local_atpair[iap].first;
        const size_t Q = local_atpair[iap].second;
        auto &chi0_mat= chi0_freq_q.at(J).at(Q);
        //printf("| IN cal Pi process %d, I: %d  J: %d  Q: %d\n",mpi_comm_global_h.myid, I,J,Q );
        auto tmp_pi_mat= *Vq_loc.at(I).at(J).at(q) * chi0_mat;
        ComplexMatrix chi0_QJ = transpose(chi0_mat, 1);
        auto tmp_pi_mat2= *Vq_loc.at(I).at(Q).at(q) * chi0_QJ;
        omp_set_lock(&pi_lock);
        pi.at(Q)+=tmp_pi_mat;
        if(J!=Q)
            {
                pi.at(J)+=tmp_pi_mat2;
            }
        omp_unset_lock(&pi_lock);
    }
    omp_destroy_lock(&pi_lock);
    // for (auto &J_p : chi0_freq_q)
    // {
    //     const size_t J = J_p.first;
    //     for (auto &Q_p : J_p.second)
    //     {
    //         const size_t Q = Q_p.first;
    //         auto &chi0_mat = Q_p.second;
    //         pi.at(Q) += Vq_row.at(I).at(J) * chi0_mat;
    //         if (J != Q)
    //         {
    //             ComplexMatrix chi0_QJ = transpose(chi0_mat, 1);
    //             pi.at(J) += Vq_row.at(I).at(Q) * chi0_QJ;
    //         }
    //     }
    // }
    //Pi_rowT
    // ComplexMatrix pi_row(N_all_mu,atom_mu[I]);
    // complex<double> *pi_row_ptr=pi_row.c;
    // for(auto &Jp:pi)
    // {
    //     auto J=Jp.first;
    //     auto J_mu=atom_mu[J];
    //     const auto length=sizeof(complex<double>)* I_mu *J_mu;
    //     memcpy(pi_row_ptr, pi.at(J).c,length);
    //     pi_row_ptr+=I_mu *J_mu;
    // }
    ComplexMatrix pi_row(atom_mu[I],N_all_mu);
    for(int i=0;i!=pi_row.nr;i++)
        for(int J=0;J!=natom;J++)
            for(int j=0;j!=atom_mu[J];j++)
                pi_row(i,atom_mu_part_range[J]+j)=pi.at(J)(i,j);
    return pi_row;
}

atom_mapping<ComplexMatrix>::pair_t_old gather_vq_row_q(const int &I, const atpair_k_cplx_mat_t &coulmat, const Vector3_Order<double> &ik_vec)
{
    auto I_mu=atom_mu[I];
    atom_mapping<ComplexMatrix>::pair_t_old Vq_row;
    for(int J_tmp=0;J_tmp!=natom;J_tmp++)
    {
        auto J_mu=atom_mu[J_tmp];
        ComplexMatrix loc_vq(atom_mu[I],atom_mu[J_tmp]);
        Vq_row[I][J_tmp].create(atom_mu[I],atom_mu[J_tmp]);
       // const auto length=sizeof(complex<double>)* I_mu *J_mu;
       // complex<double> *loc_vq_ptr=loc_vq.c;
        if(I<=J_tmp)
        {
            if(Vq.count(I))
                if(Vq.at(I).count(J_tmp))
                    loc_vq=*Vq.at(I).at(J_tmp).at(ik_vec);
        }
        else
        {
            if(Vq.count(J_tmp))
                if(Vq.at(J_tmp).count(I))
                    loc_vq=transpose(*Vq.at(J_tmp).at(I).at(ik_vec), 1);

        }
        mpi_comm_global_h.allreduce_ComplexMatrix(loc_vq,Vq_row[I][J_tmp]);
    }
    return Vq_row;
}

map<double, atom_mapping<std::map<Vector3_Order<double>, matrix_m<complex<double>>>>::pair_t_old>
compute_Wc_freq_q(const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat_eps, atpair_k_cplx_mat_t &coulmat_wc, const vector<std::complex<double>> &epsmac_LF_imagfreq)
{
    map<double, atom_mapping<std::map<Vector3_Order<double>, matrix_m<complex<double>>>>::pair_t_old> Wc_freq_q;
    const int range_all = LIBRPA::atomic_basis_abf.nb_total;
    const auto part_range = LIBRPA::atomic_basis_abf.get_part_range();

    if (mpi_comm_global_h.myid == 0)
    {
        cout << "Calculating Wc using LAPACK" << endl;
    }

    mpi_comm_global_h.barrier();
    // use q-points as the outmost loop, so that square root of Coulomb will not be recalculated at each frequency point
    vector<Vector3_Order<double>> qpts;
    for ( const auto &qMuNuchi: chi0.get_chi0_q().at(chi0.tfg.get_freq_nodes()[0]))
        qpts.push_back(qMuNuchi.first);

    for (const auto &q: qpts)
    {
        int iq = std::distance(klist.begin(), std::find(klist.begin(), klist.end(), q));
        char fn[80];

        ComplexMatrix Vq_all(range_all, range_all);
        for (const auto &Mu_NuqVq: coulmat_eps)
        {
            auto Mu = Mu_NuqVq.first;
            auto n_mu = atom_mu[Mu];
            for ( auto &Nu_qVq: Mu_NuqVq.second )
            {
                auto Nu = Nu_qVq.first;
                if ( 0 == Nu_qVq.second.count(q) ) continue;
                auto n_nu = atom_mu[Nu];
                for ( int i_mu = 0; i_mu != n_mu; i_mu++ )
                    for ( int i_nu = 0; i_nu != n_nu; i_nu++ )
                    {
                        Vq_all(part_range[Mu] + i_mu, part_range[Nu] + i_nu) = (*Nu_qVq.second.at(q))(i_mu, i_nu);
                        Vq_all(part_range[Nu] + i_nu, part_range[Mu] + i_mu) = conj((*Nu_qVq.second.at(q))(i_mu, i_nu));
                    }
            }
        }
        if (Params::debug)
        {
            sprintf(fn, "Vq_all_q_%d.mtx", iq);
            print_complex_matrix_mm(Vq_all, Params::output_dir + "/" + fn, 1e-15);
        }
        auto sqrtVq_all = power_hemat(Vq_all, 0.5, true, false, Params::sqrt_coulomb_threshold);
        // Vq_all is now eigenvectors of the original Coulomb matrix
        const auto& Vq_eigen = Vq_all;
        if (Params::debug)
        {
            sprintf(fn, "sqrtVq_all_q_%d.mtx", iq);
            print_complex_matrix_mm(sqrtVq_all, Params::output_dir + "/" + fn, 1e-15);
            // sprintf(fn, "rotated_sqrtVq_all_q_%d.mtx", iq);
            // print_complex_matrix_mm(Vq_all * sqrtVq_all * transpose(Vq_all, true), fn, 1e-15);
            // print_complex_matrix_mm(transpose(Vq_all, true) * sqrtVq_all * Vq_all, fn, 1e-15);
            sprintf(fn, "Vqeigenvec_q_%d.mtx", iq);
            print_complex_matrix_mm(Vq_eigen, Params::output_dir + "/" + fn, 1e-15);
        }

        // truncated (cutoff) Coulomb
        ComplexMatrix Vqcut_all(range_all, range_all);
        for ( auto &Mu_NuqVq: coulmat_wc )
        {
            auto Mu = Mu_NuqVq.first;
            auto n_mu = atom_mu[Mu];
            for ( auto &Nu_qVq: Mu_NuqVq.second )
            {
                auto Nu = Nu_qVq.first;
                if ( 0 == Nu_qVq.second.count(q) ) continue;
                auto n_nu = atom_mu[Nu];
                for ( int i_mu = 0; i_mu != n_mu; i_mu++ )
                    for ( int i_nu = 0; i_nu != n_nu; i_nu++ )
                    {
                        Vqcut_all(part_range[Mu] + i_mu, part_range[Nu] + i_nu) = (*Nu_qVq.second.at(q))(i_mu, i_nu);
                        Vqcut_all(part_range[Nu] + i_nu, part_range[Mu] + i_mu) = conj((*Nu_qVq.second.at(q))(i_mu, i_nu));
                    }
            }
        }
        auto sqrtVqcut_all = power_hemat(Vqcut_all, 0.5, false, true, Params::sqrt_coulomb_threshold);
        // sprintf(fn, "sqrtVqcut_all_q_%d.mtx", iq);
        // print_complex_matrix_mm(sqrtVqcut_all, fn, 1e-15);
        sprintf(fn, "Vqcut_all_filtered_q_%d.mtx", iq);
        // print_complex_matrix_mm(Vqcut_all, fn, 1e-15);
        // save the filtered truncated Coulomb back to the atom mapping object
        // TODO: revise the necessity
        for ( auto &Mu_NuqVq: coulmat_wc )
        {
            auto Mu = Mu_NuqVq.first;
            auto n_mu = atom_mu[Mu];
            for ( auto &Nu_qVq: Mu_NuqVq.second )
            {
                auto Nu = Nu_qVq.first;
                if ( 0 == Nu_qVq.second.count(q) ) continue;
                auto n_nu = atom_mu[Nu];
                for ( int i_mu = 0; i_mu != n_mu; i_mu++ )
                    for ( int i_nu = 0; i_nu != n_nu; i_nu++ )
                        (*Nu_qVq.second.at(q))(i_mu, i_nu) = Vqcut_all(part_range[Mu] + i_mu, part_range[Nu] + i_nu);
            }
        }

        ComplexMatrix chi0fq_all(range_all, range_all);
        for (const auto &freq_qMuNuchi: chi0.get_chi0_q())
        {
            auto freq = freq_qMuNuchi.first;
            auto ifreq = chi0.tfg.get_freq_index(freq);
            auto MuNuchi = freq_qMuNuchi.second.at(q);
            for (const auto &Mu_Nuchi: MuNuchi)
            {
                auto Mu = Mu_Nuchi.first;
                auto n_mu = atom_mu[Mu];
                for ( auto &Nu_chi: Mu_Nuchi.second )
                {
                    auto Nu = Nu_chi.first;
                    auto n_nu = atom_mu[Nu];
                    for ( int i_mu = 0; i_mu != n_mu; i_mu++ )
                        for ( int i_nu = 0; i_nu != n_nu; i_nu++ )
                        {
                            chi0fq_all(part_range[Mu] + i_mu, part_range[Nu] + i_nu) = Nu_chi.second(i_mu, i_nu);
                            chi0fq_all(part_range[Nu] + i_nu, part_range[Mu] + i_mu) = conj(Nu_chi.second(i_mu, i_nu));
                        }
                }
            }
            sprintf(fn, "chi0fq_all_q_%d_freq_%d.mtx", iq, ifreq);
            print_complex_matrix_mm(chi0fq_all, Params::output_dir + "/" + fn, 1e-15);

            ComplexMatrix identity(range_all, range_all);
            identity.set_as_identity_matrix();
            auto eps_fq = sqrtVq_all * chi0fq_all * sqrtVq_all;
            eps_fq = transpose(Vq_eigen, true) * eps_fq * Vq_eigen;
            if (!epsmac_LF_imagfreq.empty() && is_gamma_point(q))
            {
                // rotate to Coulomb-diagonal basis
                // lib_printf("Largest off-diagonal = %f\n", eps_fq.get_max_abs_offdiag());
                // print_matrix("rotated eps_fq: ", eps_fq.real());
                // replacing the element corresponding to largest Coulomb eigenvalue with dielectric function
                lib_printf("%22.12f %22.12f %22.12f %22.12f\n", freq, eps_fq(0, 0).real(), eps_fq(eps_fq.nr - 1, eps_fq.nc - 1).real(), epsmac_LF_imagfreq[ifreq].real());
                // eps_fq(eps_fq.nr - 1, eps_fq.nc - 1) = epsmac_LF_imagfreq[ifreq];
                eps_fq(0, 0) = 1.0 - epsmac_LF_imagfreq[ifreq];
            }
            if (Params::debug)
            {
                sprintf(fn, "rotated_vsxvs_q_%d_freq_%d.mtx", iq, ifreq);
                print_complex_matrix_mm(eps_fq, Params::output_dir + "/" + fn, 1e-10);
            }
            // rotate back to ABF
            eps_fq = Vq_eigen * eps_fq * transpose(Vq_eigen, true);
            eps_fq = identity - eps_fq;
            if (Params::debug)
            {
                sprintf(fn, "eps_q_%d_freq_%d.mtx", iq, ifreq);
                print_complex_matrix_mm(eps_fq, Params::output_dir + "/" + fn, 1e-10);
            }

            // invert the epsilon matrix
            power_hemat_onsite(eps_fq, -1);
            auto wc_all = sqrtVqcut_all * (eps_fq - identity) * sqrtVqcut_all;
            // sprintf(fn, "inveps_q_%d_freq_%d.mtx", iq, ifreq);
            // print_complex_matrix_mm(eps_fq, fn, 1e-15);
            // sprintf(fn, "wc_q_%d_freq_%d.mtx", iq, ifreq);
            // print_complex_matrix_mm(wc_all, fn, 1e-15);

            // save result to the atom mapping object
            for ( auto &Mu_Nuchi: MuNuchi )
            {
                auto Mu = Mu_Nuchi.first;
                auto n_mu = atom_mu[Mu];
                for ( auto &Nu_chi: Mu_Nuchi.second )
                {
                    auto Nu = Nu_chi.first;
                    auto n_nu = atom_mu[Nu];
                    shared_ptr<ComplexMatrix> wc_ptr = make_shared<ComplexMatrix>();
                    wc_ptr->create(n_mu, n_nu);
                    for ( int i_mu = 0; i_mu != n_mu; i_mu++ )
                        for ( int i_nu = 0; i_nu != n_nu; i_nu++ )
                        {
                            (*wc_ptr)(i_mu, i_nu) = wc_all(part_range[Mu] + i_mu, part_range[Nu] + i_nu);
                        }
                    Wc_freq_q[freq][Mu][Nu][q] = matrix_m<complex<double>>(n_mu, n_nu, wc_ptr->c, MAJOR::ROW, MAJOR::ROW);
                }
            }
        }
    }

    return Wc_freq_q;
}

map<double, atom_mapping<std::map<Vector3_Order<double>, matrix_m<complex<double>>>>::pair_t_old>
compute_Wc_freq_q_blacs(const Chi0 &chi0, const atpair_k_cplx_mat_t &coulmat_eps, atpair_k_cplx_mat_t &coulmat_wc, const vector<std::complex<double>> &epsmac_LF_imagfreq)
{
    map<double, atom_mapping<std::map<Vector3_Order<double>, matrix_m<complex<double>>>>::pair_t_old> Wc_freq_q;
    const complex<double> CONE{1.0, 0.0};
    const int n_abf = LIBRPA::atomic_basis_abf.nb_total;
    const auto part_range = LIBRPA::atomic_basis_abf.get_part_range();

    if (mpi_comm_global_h.myid == 0)
    {
        cout << "Calculating Wc using ScaLAPACK" << endl;
    }
    mpi_comm_global_h.barrier();

    Array_Desc desc_nabf_nabf(blacs_ctxt_global_h);
    // use a square blocksize instead max block, otherwise heev and inversion will complain about illegal parameter
    desc_nabf_nabf.init_square_blk(n_abf, n_abf, 0, 0);
    // obtain the indices of atom-pair block necessary to build 2D block of a Hermitian/symmetric matrix
    const auto set_IJ_nabf_nabf = LIBRPA::get_necessary_IJ_from_block_2D_sy('U', LIBRPA::atomic_basis_abf, desc_nabf_nabf);
    const auto s0_s1 = get_s0_s1_for_comm_map2_first(set_IJ_nabf_nabf);
    auto chi0_block = init_local_mat<complex<double>>(desc_nabf_nabf, MAJOR::COL);
    auto coul_block = init_local_mat<complex<double>>(desc_nabf_nabf, MAJOR::COL);
    auto coul_eigen_block = init_local_mat<complex<double>>(desc_nabf_nabf, MAJOR::COL);
    auto coul_chi0_block = init_local_mat<complex<double>>(desc_nabf_nabf, MAJOR::COL);
    auto coulwc_block = init_local_mat<complex<double>>(desc_nabf_nabf, MAJOR::COL);

    const auto atpair_local = dispatch_upper_trangular_tasks(
        natom, blacs_ctxt_global_h.myid, blacs_ctxt_global_h.nprows,
        blacs_ctxt_global_h.npcols, blacs_ctxt_global_h.myprow,
        blacs_ctxt_global_h.mypcol);
    ofs_myid << "atpair_local " << blacs_ctxt_global_h.myid << " " << atpair_local << endl;
    std::flush(ofs_myid);

    // IJ pair of Wc to be returned
    pair<set<int>, set<int>> Iset_Jset_Wc;
    for (const auto &ap: atpair_local)
    {
        Iset_Jset_Wc.first.insert(ap.first);
        Iset_Jset_Wc.second.insert(ap.second);
    }

    vector<Vector3_Order<double>> qpts;
    for (const auto &qMuNuchi: chi0.get_chi0_q().at(chi0.tfg.get_freq_nodes()[0]))
        qpts.push_back(qMuNuchi.first);

    vec<double> eigenvalues(n_abf);

#ifdef LIBRPA_USE_LIBRI
    for (const auto &q: qpts)
    {
        coul_block.zero_out();
        coulwc_block.zero_out();
        // lib_printf("coul_block\n%s", str(coul_block).c_str());

        int iq = std::distance(klist.begin(), std::find(klist.begin(), klist.end(), q));
        std::array<double, 3> qa = {q.x, q.y, q.z};

        // collect the block elements of truncated coulomb matrices first
        // as we reuse coul_eigen_block to reduce memory usage
        Profiler::start("epsilon_prepare_coulwc_sqrt", "Prepare sqrt of truncated Coulomb");
        {
            size_t n_singular_coulwc;
            // LibRI tensor for communication, release once done
            std::map<int, std::map<std::pair<int, std::array<double, 3>>, RI::Tensor<complex<double>>>> couleps_libri;
            Profiler::start("epsilon_prepare_coulwc_sqrt_1", "Setup libRI object");
            for (const auto &Mu_Nu: atpair_local)
            {
                const auto Mu = Mu_Nu.first;
                const auto Nu = Mu_Nu.second;
                ofs_myid << "Mu " << Mu << " Nu " << Nu << endl;
                if (coulmat_wc.count(Mu) == 0 ||
                    coulmat_wc.at(Mu).count(Nu) == 0 ||
                    coulmat_wc.at(Mu).at(Nu).count(q) == 0) continue;
                const auto &Vq = coulmat_wc.at(Mu).at(Nu).at(q);
                const auto n_mu = LIBRPA::atomic_basis_abf.get_atom_nb(Mu);
                const auto n_nu = LIBRPA::atomic_basis_abf.get_atom_nb(Nu);
                std::valarray<complex<double>> Vq_va(Vq->c, Vq->size);
                auto pvq = std::make_shared<std::valarray<complex<double>>>();
                *pvq = Vq_va;
                couleps_libri[Mu][{Nu, qa}] = RI::Tensor<complex<double>>({n_mu, n_nu}, pvq);
            }
            Profiler::stop("epsilon_prepare_coulwc_sqrt_1");

            Profiler::start("epsilon_prepare_coulwc_sqrt_2", "libRI Communicate");
            const auto IJq_coul = RI::Communicate_Tensors_Map_Judge::comm_map2_first(mpi_comm_global_h.comm, couleps_libri, s0_s1.first, s0_s1.second);
            Profiler::stop("epsilon_prepare_coulwc_sqrt_2");

            Profiler::start("epsilon_prepare_coulwc_sqrt_3", "Collect 2D-block from IJ");
            // for (const auto &IJ: set_IJ_nabf_nabf)
            // {
            //     const auto &I = IJ.first;
            //     const auto &J = IJ.second;
            //     collect_block_from_IJ_storage_syhe(
            //         coulwc_block, desc_nabf_nabf, LIBRPA::atomic_basis_abf, IJ.first,
            //         IJ.second, true, CONE, IJq_coul.at(I).at({J, qa}).ptr(), MAJOR::ROW);
            // }
            collect_block_from_ALL_IJ_Tensor(coulwc_block, desc_nabf_nabf, LIBRPA::atomic_basis_abf,
                                             qa, true, CONE, IJq_coul, MAJOR::ROW);
            Profiler::stop("epsilon_prepare_coulwc_sqrt_3");
            Profiler::start("epsilon_prepare_coulwc_sqrt_4", "Perform square root");
            power_hemat_blacs(coulwc_block, desc_nabf_nabf, coul_eigen_block, desc_nabf_nabf, n_singular_coulwc, eigenvalues.c, 0.5, Params::sqrt_coulomb_threshold);
            Profiler::stop("epsilon_prepare_coulwc_sqrt_4");
        }
        Profiler::stop("epsilon_prepare_coulwc_sqrt");
        lib_printf("Task %d: done coulwc sqrt\n", mpi_comm_global_h.myid);

        Profiler::start("epsilon_prepare_couleps_sqrt", "Prepare sqrt of bare Coulomb");
        // collect the block elements of coulomb matrices
        {
            // LibRI tensor for communication, release once done
            std::map<int, std::map<std::pair<int, std::array<double, 3>>, RI::Tensor<complex<double>>>> couleps_libri;
            ofs_myid << "Start build couleps_libri\n";
            for (const auto &Mu_Nu: atpair_local)
            {
                const auto Mu = Mu_Nu.first;
                const auto Nu = Mu_Nu.second;
                ofs_myid << "Mu " << Mu << " Nu " << Nu << endl;
                if (coulmat_eps.count(Mu) == 0 ||
                    coulmat_eps.at(Mu).count(Nu) == 0 ||
                    coulmat_eps.at(Mu).at(Nu).count(q) == 0) continue;
                const auto &Vq = coulmat_eps.at(Mu).at(Nu).at(q);
                const auto n_mu = LIBRPA::atomic_basis_abf.get_atom_nb(Mu);
                const auto n_nu = LIBRPA::atomic_basis_abf.get_atom_nb(Nu);
                std::valarray<complex<double>> Vq_va(Vq->c, Vq->size);
                auto pvq = std::make_shared<std::valarray<complex<double>>>();
                *pvq = Vq_va;
                couleps_libri[Mu][{Nu, qa}] = RI::Tensor<complex<double>>({n_mu, n_nu}, pvq);
            }
            ofs_myid << "Done build couleps_libri\n";
            // ofs_myid << "Couleps_libri" << endl << couleps_libri;
            // if (couleps_libri.size() == 0)
            //     throw std::logic_error("data at q-point not found in coulmat_eps");

            // perform communication
            ofs_myid << "Start collect couleps_libri, targets\n";
            ofs_myid << set_IJ_nabf_nabf << "\n";
            ofs_myid << "Extended blocks\n";
            ofs_myid << "atom 1: " << s0_s1.first << "\n";
            ofs_myid << "atom 2: " << s0_s1.second << "\n";
            ofs_myid << "Owned blocks\n";
            print_keys(ofs_myid, couleps_libri);
            std::flush(ofs_myid);
            mpi_comm_global_h.barrier();
            const auto IJq_coul = RI::Communicate_Tensors_Map_Judge::comm_map2_first(mpi_comm_global_h.comm, couleps_libri, s0_s1.first, s0_s1.second);
            ofs_myid << "Done collect couleps_libri, collected blocks\n";
            print_keys(ofs_myid, IJq_coul);
            std::flush(ofs_myid);
            mpi_comm_global_h.barrier();
            // ofs_myid << "IJq_coul" << endl << IJq_coul;

            ofs_myid << "Start construct couleps 2D block\n";
            std::flush(ofs_myid);
            mpi_comm_global_h.barrier();
            for (const auto &IJ: set_IJ_nabf_nabf)
            {
                const auto &I = IJ.first;
                const auto &J = IJ.second;
                try
                {
                    IJq_coul.at(I).at({J, qa});
                }
                catch (const std::out_of_range& e)
                {
                    ofs_myid << "Fail to find " << I << " " << J << " " << qa << ": " << e.what() << "\n";
                    std::flush(ofs_myid);
                }
                // collect_block_from_IJ_storage_syhe(
                //     coul_block, desc_nabf_nabf, LIBRPA::atomic_basis_abf, IJ.first,
                //     IJ.second, true, CONE, IJq_coul.at(I).at({J, qa}).ptr(), MAJOR::ROW);
                collect_block_from_ALL_IJ_Tensor(coul_block, desc_nabf_nabf, LIBRPA::atomic_basis_abf,
                                                 qa, true, CONE, IJq_coul, MAJOR::ROW);
                // lib_printf("myid %d I %d J %d nr %d nc %d\n%s",
                //        blacs_ctxt_global_h.myid, I, J,
                //        coul_block.nr(), coul_block.nc(),
                //        str(coul_block).c_str());
            }
            ofs_myid << "Done construct couleps 2D block\n";
            std::flush(ofs_myid);
            mpi_comm_global_h.barrier();
        }
        // char fn[100];
        // sprintf(fn, "couleps_iq_%d.mtx", iq);
        // print_matrix_mm_file_parallel(fn, coul_block, desc_nabf_nabf);
        // ofs_myid << str(coul_block);
        // lib_printf("coul_block\n%s", str(coul_block).c_str());

        size_t n_singular;
        lib_printf("Task %d: start power hemat couleps\n", mpi_comm_global_h.myid);
        auto sqrtveig_blacs = power_hemat_blacs(coul_block, desc_nabf_nabf, coul_eigen_block, desc_nabf_nabf, n_singular, eigenvalues.c, 0.5, Params::sqrt_coulomb_threshold);
        lib_printf("Task %d: done power hemat couleps\n", mpi_comm_global_h.myid);
        // lib_printf("nabf %d nsingu %lu\n", n_abf, n_singular);
        // release sqrtv when the q-point is not Gamma, or macroscopic dielectric constant at imaginary frequency is not prepared
        if (epsmac_LF_imagfreq.empty() || !is_gamma_point(q))
            sqrtveig_blacs.clear();
        const size_t n_nonsingular = n_abf - n_singular;
        Profiler::stop("epsilon_prepare_couleps_sqrt");
        lib_printf("Task %d: done couleps sqrt\n", mpi_comm_global_h.myid);

        for (const auto &freq: chi0.tfg.get_freq_nodes())
        {
            const auto ifreq = chi0.tfg.get_freq_index(freq);
            Profiler::start("epsilon_prepare_chi0_2d", "Prepare Chi0 2D block");
            chi0_block.zero_out();
            {
                std::map<int, std::map<std::pair<int, std::array<double, 3>>, RI::Tensor<complex<double>>>> chi0_libri;
                const auto &chi0_wq = chi0.get_chi0_q().at(freq).at(q);
                for (const auto &M_Nchi: chi0_wq)
                {
                    const auto &M = M_Nchi.first;
                    const auto n_mu = LIBRPA::atomic_basis_abf.get_atom_nb(M);
                    for (const auto &N_chi: M_Nchi.second)
                    {
                        const auto &N = N_chi.first;
                        const auto n_nu = LIBRPA::atomic_basis_abf.get_atom_nb(N);
                        const auto &chi = N_chi.second;
                        std::valarray<complex<double>> chi_va(chi.c, chi.size);
                        auto pchi = std::make_shared<std::valarray<complex<double>>>();
                        *pchi = chi_va;
                        chi0_libri[M][{N, qa}] = RI::Tensor<complex<double>>({n_mu, n_nu}, pchi);
                    }
                }
                // ofs_myid << "chi0_libri" << endl << chi0_libri;
                const auto IJq_chi0 = RI::Communicate_Tensors_Map_Judge::comm_map2_first(mpi_comm_global_h.comm, chi0_libri, s0_s1.first, s0_s1.second);
                // ofs_myid << "IJq_chi0" << endl << IJq_chi0;
                // for (const auto &IJ: set_IJ_nabf_nabf)
                // {
                //     const auto &I = IJ.first;
                //     const auto &J = IJ.second;
                //     collect_block_from_IJ_storage_syhe(
                //         chi0_block, desc_nabf_nabf, LIBRPA::atomic_basis_abf, IJ.first,
                //         IJ.second, true, CONE, IJq_chi0.at(I).at({J, qa}).ptr(), MAJOR::ROW);
                // }
                collect_block_from_ALL_IJ_Tensor(chi0_block, desc_nabf_nabf, LIBRPA::atomic_basis_abf,
                                                 qa, true, CONE, IJq_chi0, MAJOR::ROW);
                // sprintf(fn, "chi_ifreq_%d_iq_%d.mtx", ifreq, iq);
                // print_matrix_mm_file_parallel(fn, chi0_block, desc_nabf_nabf);
            }
            Profiler::stop("epsilon_prepare_chi0_2d");

            Profiler::start("epsilon_compute_eps", "Compute dielectric matrix");
            // for Gamma point, overwrite the head term
            if (epsmac_LF_imagfreq.size() > 0 && is_gamma_point(q))
            {
                ofs_myid << "Entering dielectric matrix head overwrite\n";
                if (Params::debug)
                {
                    lib_printf("Task %4d: Entering dielectric matrix head overwrite\n", mpi_comm_global_h.myid);
                }
                // rotate to Coulomb-eigenvector basis
                ScalapackConnector::pgemm_f('N', 'N', n_abf, n_nonsingular, n_abf, 1.0,
                        chi0_block.ptr(), 1, 1, desc_nabf_nabf.desc,
                        sqrtveig_blacs.ptr(), 1, n_singular + 1, desc_nabf_nabf.desc, 0.0,
                        coul_chi0_block.ptr(), 1, 1, desc_nabf_nabf.desc);
                ScalapackConnector::pgemm_f('C', 'N', n_nonsingular, n_nonsingular, n_abf, 1.0,
                        sqrtveig_blacs.ptr(), 1, n_singular + 1, desc_nabf_nabf.desc,
                        coul_chi0_block.ptr(), 1, 1, desc_nabf_nabf.desc, 0.0,
                        chi0_block.ptr(), 1, 1, desc_nabf_nabf.desc);
                const int ilo = desc_nabf_nabf.indx_g2l_r(n_nonsingular - 1);
                const int jlo = desc_nabf_nabf.indx_g2l_c(n_nonsingular - 1);
                if (ilo >= 0 && jlo >= 0)
                {
                    if (Params::debug)
                    {
                        lib_printf("Task %4d: perform the head element overwrite\n", mpi_comm_global_h.myid);
                    }
                    chi0_block(ilo, jlo) = 1.0 - epsmac_LF_imagfreq[ifreq];
                }
                // rotate back to ABF
                ScalapackConnector::pgemm_f('N', 'N', n_abf, n_nonsingular, n_nonsingular, 1.0,
                        coul_eigen_block.ptr(), 1, n_singular + 1, desc_nabf_nabf.desc,
                        chi0_block.ptr(), 1, 1, desc_nabf_nabf.desc, 0.0,
                        coul_chi0_block.ptr(), 1, 1, desc_nabf_nabf.desc);
                ScalapackConnector::pgemm_f('N', 'C', n_abf, n_abf, n_nonsingular, 1.0,
                        coul_chi0_block.ptr(), 1, 1, desc_nabf_nabf.desc,
                        coul_eigen_block.ptr(), 1, n_singular + 1, desc_nabf_nabf.desc, 0.0,
                        chi0_block.ptr(), 1, 1, desc_nabf_nabf.desc);
            }
            else
            {
                ScalapackConnector::pgemm_f('N', 'N', n_abf, n_abf, n_abf, 1.0,
                        coul_block.ptr(), 1, 1, desc_nabf_nabf.desc,
                        chi0_block.ptr(), 1, 1, desc_nabf_nabf.desc, 0.0,
                        coul_chi0_block.ptr(), 1, 1, desc_nabf_nabf.desc);
                ScalapackConnector::pgemm_f('N', 'N', n_abf, n_abf, n_abf, 1.0,
                        coul_chi0_block.ptr(), 1, 1, desc_nabf_nabf.desc,
                        coul_block.ptr(), 1, 1, desc_nabf_nabf.desc, 0.0,
                        chi0_block.ptr(), 1, 1, desc_nabf_nabf.desc);
            }
            // now chi0_block is actually v1/2 chi v1/2
            chi0_block *= -1.0;
            for (int i = 0; i != n_abf; i++)
            {
                const int ilo = desc_nabf_nabf.indx_g2l_r(i);
                if (ilo < 0) continue;
                const int jlo = desc_nabf_nabf.indx_g2l_c(i);
                if (jlo < 0) continue;
                chi0_block(ilo, jlo) += 1.0;
            }
            Profiler::stop("epsilon_compute_eps");
            // now chi0_block is actually the dielectric matrix
            // perform inversion
            Profiler::start("epsilon_invert_eps", "Invert dielectric matrix");
            invert_scalapack(chi0_block, desc_nabf_nabf);
            // subtract 1 from diagonal
            for (int i = 0; i != n_abf; i++)
            {
                const int ilo = desc_nabf_nabf.indx_g2l_r(i);
                if (ilo < 0) continue;
                const int jlo = desc_nabf_nabf.indx_g2l_c(i);
                if (jlo < 0) continue;
                chi0_block(ilo, jlo) -= 1.0;
            }
            Profiler::stop("epsilon_invert_eps");

            Profiler::start("epsilon_multiply_coulwc", "Multiply truncated Coulomb");
            ScalapackConnector::pgemm_f('N', 'N', n_abf, n_abf, n_abf, 1.0,
                    coulwc_block.ptr(), 1, 1, desc_nabf_nabf.desc,
                    chi0_block.ptr(), 1, 1, desc_nabf_nabf.desc, 0.0,
                    coul_chi0_block.ptr(), 1, 1, desc_nabf_nabf.desc);
            ScalapackConnector::pgemm_f('N', 'N', n_abf, n_abf, n_abf, 1.0,
                    coul_chi0_block.ptr(), 1, 1, desc_nabf_nabf.desc,
                    coulwc_block.ptr(), 1, 1, desc_nabf_nabf.desc, 0.0,
                    chi0_block.ptr(), 1, 1, desc_nabf_nabf.desc);
            Profiler::stop("epsilon_multiply_coulwc");
            // lib_printf("chi0_block\n%s", str(chi0_block).c_str());
            // now chi0_block is the screened Coulomb interaction Wc (i.e. W-V)

            Profiler::start("epsilon_convert_wc_2d_to_ij", "Convert Wc, 2D -> IJ");
            Profiler::start("epsilon_convert_wc_map_block", "Initialize Wc atom-pair map");
            map<int, map<int, matrix_m<complex<double>>>> Wc_MNmap;
            map_block_to_IJ_storage(Wc_MNmap, LIBRPA::atomic_basis_abf,
                                    LIBRPA::atomic_basis_abf, chi0_block,
                                    desc_nabf_nabf, MAJOR::ROW);
            Profiler::stop("epsilon_convert_wc_map_block");

            Profiler::start("epsilon_convert_wc_communicate", "Communicate");
            {
                std::map<int, std::map<std::pair<int, std::array<double, 3>>, RI::Tensor<complex<double>>>> Wc_libri;
                for (const auto &M_NWc: Wc_MNmap)
                {
                    const auto &M = M_NWc.first;
                    const auto n_mu = LIBRPA::atomic_basis_abf.get_atom_nb(M);
                    for (const auto &N_Wc: M_NWc.second)
                    {
                        const auto &N = N_Wc.first;
                        const auto n_nu = LIBRPA::atomic_basis_abf.get_atom_nb(N);
                        const auto &Wc = N_Wc.second;
                        // std::valarray<complex<double>> Wc_va(Wc.ptr(), Wc.size());
                        // auto pWc = std::make_shared<std::valarray<complex<double>>>();
                        // *pWc = Wc_va;
                        Wc_libri[M][{N, qa}] = RI::Tensor<complex<double>>({n_mu, n_nu}, Wc.sptr());
                    }
                }
                // cout << Wc_libri;
                const auto IJq_Wc = RI::Communicate_Tensors_Map_Judge::comm_map2_first(mpi_comm_global_h.comm, Wc_libri, Iset_Jset_Wc.first, Iset_Jset_Wc.second);
                // parse collected to
                for (const auto &MN: atpair_local)
                {
                    const auto &M = MN.first;
                    const auto &N = MN.second;
                    const auto n_mu = LIBRPA::atomic_basis_abf.get_atom_nb(M);
                    const auto n_nu = LIBRPA::atomic_basis_abf.get_atom_nb(N);
                    // Use row major for later usage in LibRI
                    Wc_freq_q[freq][M][N][q] = matrix_m<complex<double>>(
                        n_mu, n_nu, IJq_Wc.at(M).at({N, qa}).data, MAJOR::ROW);
                }
                // for ( int i_mu = 0; i_mu != n_mu; i_mu++ )
                //     for ( int i_nu = 0; i_nu != n_nu; i_nu++ )
                //     {
                //     }
            }
            Profiler::stop("epsilon_convert_wc_communicate");
            Profiler::stop("epsilon_convert_wc_2d_to_ij");
        }
    }
#else
    throw std::logic_error("need compilation with LibRI");
#endif
    return Wc_freq_q;
}

map<double, atom_mapping<std::map<Vector3_Order<int>, matrix_m<complex<double>>>>::pair_t_old>
FT_Wc_freq_q(const map<double, atom_mapping<std::map<Vector3_Order<double>, matrix_m<complex<double>>>>::pair_t_old> &Wc_freq_q,
             const TFGrids &tfg, vector<Vector3_Order<int>> Rlist)
{
    // major of Wc_freq_q input and Wc_tau_R output
    const MAJOR major_Wc = MAJOR::ROW;

    map<double, atom_mapping<std::map<Vector3_Order<int>, matrix_m<complex<double>>>>::pair_t_old> Wc_freq_R;
    const int ngrids = tfg.get_n_grids();
    if (Params::debug)
    {
        if (mpi_comm_global_h.is_root())
            lib_printf("Converting Wc q,w -> R,t\n");
        mpi_comm_global_h.barrier();
    }
    set<pair<atom_t, atom_t>> atpairs_unique;
    for (const auto &freq_MuNuqWc: Wc_freq_q)
    {
        for (const auto &Mu_NuqWc: freq_MuNuqWc.second)
        {
            const auto Mu = Mu_NuqWc.first;
            for (const auto &Nu_qWc: Mu_NuqWc.second)
            {
                const auto Nu = Nu_qWc.first;
                atpairs_unique.insert({Mu, Nu});
                for (const auto &q_Wc: Nu_qWc.second)
                {
                    assert(q_Wc.second.major() == major_Wc);
                }
            }
        }
    }

    vector<pair<pair<int, Vector3_Order<int>>, pair<atom_t, atom_t>>> ifreqR_atpair_all;
    // allocate space before hand
    for (auto R: Rlist)
    {
        for (int ifreq = 0; ifreq != ngrids; ifreq++)
        {
            auto freq = tfg.get_freq_nodes()[ifreq];
            for (auto atpair_unique: atpairs_unique)
            {
                const auto Mu = atpair_unique.first;
                const int n_mu = atom_mu[Mu];
                const auto Nu = atpair_unique.second;
                const int n_nu = atom_mu[Nu];
                Wc_freq_R[freq][Mu][Nu][R] = matrix_m<complex<double>>(n_mu, n_nu, major_Wc);
                ifreqR_atpair_all.push_back({{ifreq, R}, atpair_unique});
            }
        }
    }

    #pragma omp parallel for schedule(dynamic)
    for (auto ifreqR_atpair: ifreqR_atpair_all)
    {
        const auto ifreq = ifreqR_atpair.first.first;
        const auto freq = tfg.get_freq_nodes()[ifreq];
        const auto R = ifreqR_atpair.first.second;
        const auto Mu = ifreqR_atpair.second.first;
        const auto Nu = ifreqR_atpair.second.second;
        const int n_mu = atom_mu[Mu];
        const int n_nu = atom_mu[Nu];

        // thread local temporary matrix
        matrix_m<complex<double>> WfR_temp(n_mu, n_nu, major_Wc);

        if (Wc_freq_q.count(freq) == 0) continue;
        if (Wc_freq_q.at(freq).count(Mu) == 0) continue;
        if (Wc_freq_q.at(freq).at(Mu).count(Nu) == 0) continue;

        for (auto &Wc_q: Wc_freq_q.at(freq).at(Mu).at(Nu))
        {
            const auto q = Wc_q.first;
            const auto &Wc = Wc_q.second;
            for (auto q_bz: map_irk_ks[q])
            {
                double ang = - q_bz * (R * latvec) * TWO_PI;
                complex<double> kphase = complex<double>(cos(ang), sin(ang));
                if (q == q_bz)
                    WfR_temp += Wc * kphase;
                else
                    WfR_temp += conj(Wc) * kphase;
            }
        }
        // omp_set_lock(&lock_Wc);
        Wc_freq_R[freq][Mu][Nu][R] += WfR_temp;
        // omp_unset_lock(&lock_Wc);
    }

    if (mpi_comm_global_h.is_root())
        lib_printf("Done converting Wc(q,w) -> Wc(R,w)\n");
    mpi_comm_global_h.barrier();

    return Wc_freq_R;
}


map<double, atom_mapping<std::map<Vector3_Order<int>, matrix_m<complex<double>>>>::pair_t_old>
CT_FT_Wc_freq_q(const map<double, atom_mapping<std::map<Vector3_Order<double>, matrix_m<complex<double>>>>::pair_t_old> &Wc_freq_q,
                const TFGrids &tfg, vector<Vector3_Order<int>> Rlist)
{
    // major of Wc_freq_q input and Wc_tau_R output
    const MAJOR major_Wc = MAJOR::ROW;

    map<double, atom_mapping<std::map<Vector3_Order<int>, matrix_m<complex<double>>>>::pair_t_old> Wc_tau_R;
    if (!tfg.has_time_grids())
        throw logic_error("TFGrids object does not have time grids");
    const int ngrids = tfg.get_n_grids();

    if (mpi_comm_global_h.is_root())
        lib_printf("Converting Wc(q,w) -> W(R,t)\n");
    mpi_comm_global_h.barrier();

    set<pair<atom_t, atom_t>> atpairs_unique;
    for (const auto &freq_MuNuqWc: Wc_freq_q)
    {
        for (const auto &Mu_NuqWc: freq_MuNuqWc.second)
        {
            const auto Mu = Mu_NuqWc.first;
            for (const auto &Nu_qWc: Mu_NuqWc.second)
            {
                const auto Nu = Nu_qWc.first;
                atpairs_unique.insert({Mu, Nu});
                for (const auto &q_Wc: Nu_qWc.second)
                {
                    assert(q_Wc.second.major() == major_Wc);
                }
            }
        }
    }

    vector<pair<pair<int, Vector3_Order<int>>, pair<atom_t, atom_t>>> itauR_atpair_all;
    // allocate space before hand
    for (auto R: Rlist)
    {
        for (int itau = 0; itau != ngrids; itau++)
        {
            auto tau = tfg.get_time_nodes()[itau];
            for (auto atpair_unique: atpairs_unique)
            {
                const auto Mu = atpair_unique.first;
                const int n_mu = atom_mu[Mu];
                const auto Nu = atpair_unique.second;
                const int n_nu = atom_mu[Nu];
                Wc_tau_R[tau][Mu][Nu][R] = matrix_m<complex<double>>(n_mu, n_nu, major_Wc);
                itauR_atpair_all.push_back({{itau, R}, atpair_unique});
            }
        }
    }

    #pragma omp parallel for schedule(dynamic)
    for (auto itauR_atpair: itauR_atpair_all)
    {
        const auto itau = itauR_atpair.first.first;
        const auto tau = tfg.get_time_nodes()[itau];
        const auto R = itauR_atpair.first.second;
        const auto Mu = itauR_atpair.second.first;
        const auto Nu = itauR_atpair.second.second;
        const int n_mu = atom_mu[Mu];
        const int n_nu = atom_mu[Nu];

        // thread local temporary matrix
        matrix_m<complex<double>> WtR_temp(n_mu, n_nu, major_Wc);

        for (int ifreq = 0; ifreq < ngrids; ifreq++)
        {
            const auto freq = tfg.get_freq_nodes()[ifreq];
            const auto f2t = tfg.get_costrans_f2t()(itau, ifreq);
            // ofs_myid << "f2t cos eff for freq " << freq << " -> tau " << tau  << ": " << f2t << "\n";
            if (Wc_freq_q.count(freq) == 0) continue;
            if (Wc_freq_q.at(freq).count(Mu) == 0) continue;
            if (Wc_freq_q.at(freq).at(Mu).count(Nu) == 0) continue;
            // cout << "freq: " << freq << "\n";

            const auto &Wc_q_all = Wc_freq_q.at(freq).at(Mu).at(Nu);
            for (auto &Wc_q: Wc_q_all)
            {
                const auto q = Wc_q.first;
                const auto &Wc = Wc_q.second;
                for (auto q_bz: map_irk_ks[q])
                {
                    double ang = - q_bz * (R * latvec) * TWO_PI;
                    complex<double> kphase = complex<double>(cos(ang), sin(ang));
                    complex<double> weight = kphase * f2t;
                    // ofs_myid << q << " " << q_bz << " weight = " << weight << "\n";
                    // ofs_myid << q_Wc.second;
                    if (q == q_bz)
                        WtR_temp += Wc * weight;
                    else
                        WtR_temp += conj(Wc) * weight;
                }
            }
        }
        // omp_set_lock(&lock_Wc);
        Wc_tau_R[tau][Mu][Nu][R] += WtR_temp;
        // omp_unset_lock(&lock_Wc);
    }

    if (mpi_comm_global_h.is_root())
        lib_printf("Done converting Wc q,w -> R,t\n");
    mpi_comm_global_h.barrier();

    // myz debug: check the imaginary part of the matrix
    // NOTE: if G(R) is real, is W(R) real as well?
    // if (Params::debug)
    // {
    //     for (const auto & tau_MuNuRWc: Wc_tau_R)
    //     {
    //         char fn[80];
    //         auto tau = tau_MuNuRWc.first;
    //         auto itau = tfg.get_time_index(tau);
    //         for (const auto & Mu_NuRWc: tau_MuNuRWc.second)
    //         {
    //             auto Mu = Mu_NuRWc.first;
    //             // const int n_mu = atom_mu[Mu];
    //             for (const auto & Nu_RWc: Mu_NuRWc.second)
    //             {
    //                 auto Nu = Nu_RWc.first;
    //                 // const int n_nu = atom_mu[Nu];
    //                 for (const auto & R_Wc: Nu_RWc.second)
    //                 {
    //                     auto R = R_Wc.first;
    //                     auto Wc = R_Wc.second;
    //                     auto iteR = std::find(Rlist.cbegin(), Rlist.cend(), R);
    //                     auto iR = std::distance(Rlist.cbegin(), iteR);
    //                     sprintf(fn, "Wc_Mu_%zu_Nu_%zu_iR_%zu_itau_%d_id_%d.mtx", Mu, Nu, iR, itau, mpi_comm_global_h.myid);
    //                     print_matrix_mm_file(Wc, Params::output_dir + "/" + fn, 1e-10);
    //                 }
    //             }
    //         }
    //     }
    // }
    // end myz debug
    return Wc_tau_R;
}

void test_libcomm_for_system(const atpair_k_cplx_mat_t &coulmat)
{
    if (mpi_comm_global_h.myid == 0)
        lib_printf("test_libcomm_for_system Coulumb\n");
    // lib_printf("Calculating EcRPA with BLACS, pid:  %d\n", mpi_comm_global_h.myid);
    const complex<double> CONE{1.0, 0.0};
    const int n_abf = LIBRPA::atomic_basis_abf.nb_total;
    const auto part_range = LIBRPA::atomic_basis_abf.get_part_range();

    mpi_comm_global_h.barrier();

    Array_Desc desc_nabf_nabf(blacs_ctxt_global_h);
    // use a square blocksize instead max block, otherwise heev and inversion will complain about illegal parameter
    desc_nabf_nabf.init_square_blk(n_abf, n_abf, 0, 0);
    const auto set_IJ_nabf_nabf = LIBRPA::get_necessary_IJ_from_block_2D_sy('U', LIBRPA::atomic_basis_abf, desc_nabf_nabf);
    const auto s0_s1 = get_s0_s1_for_comm_map2_first(set_IJ_nabf_nabf);
    
    auto coul_block = init_local_mat<complex<double>>(desc_nabf_nabf, MAJOR::COL);


    vector<Vector3_Order<double>> qpts;
    for (const auto &qMuNuchi: irk_weight)
        qpts.push_back(qMuNuchi.first);

#ifdef LIBRPA_USE_LIBRI
    for (const auto &q: qpts)
    {
        coul_block.zero_out();

        int iq = std::distance(klist.begin(), std::find(klist.begin(), klist.end(), q));
        std::array<double, 3> qa = {q.x, q.y, q.z};
        // collect the block elements of coulomb matrices
        {
            double vq_begin = omp_get_wtime();
            // LibRI tensor for communication, release once done
            std::map<int, std::map<std::pair<int, std::array<double, 3>>, Tensor<complex<double>>>> coul_libri;
            coul_libri.clear();
            int count_coul=0;
            for (const auto &Mu_Nu: local_atpair)
            {
                const auto Mu = Mu_Nu.first;
                const auto Nu = Mu_Nu.second;
                // ofs_myid << "myid " << blacs_ctxt_global_h.myid << "Mu " << Mu << " Nu " << Nu << endl;
                if (coulmat.count(Mu) == 0 ||
                    coulmat.at(Mu).count(Nu) == 0 ||
                    coulmat.at(Mu).at(Nu).count(q) == 0) continue;
                const auto &Vq = coulmat.at(Mu).at(Nu).at(q);
                const auto n_mu = LIBRPA::atomic_basis_abf.get_atom_nb(Mu);
                const auto n_nu = LIBRPA::atomic_basis_abf.get_atom_nb(Nu);
                std::valarray<complex<double>> Vq_va(Vq->c, Vq->size);
                auto pvq = std::make_shared<std::valarray<complex<double>>>();
                *pvq = Vq_va;
                coul_libri[Mu][{Nu,qa}] = Tensor<complex<double>>({n_mu, n_nu}, pvq);
                count_coul+=1;
            }
            int count_pair=0;
            for(auto &Mu:coul_libri)
            {
                for(auto &nu_q:Mu.second)
                {
                    count_pair+=1;
                }
            }
            //printf("Finish RPA blacs 2d  vq arr\n");
            double arr_end = omp_get_wtime();
            mpi_comm_global_h.barrier();
            double comm_begin = omp_get_wtime();
            lib_printf("Begin comm_map2_first  myid: %d  q:(%f, %f, %f)  count_coul: %d  count_pair: %d\n",mpi_comm_global_h.myid,q.x,q.y,q.z,count_coul,count_pair);
            const auto IJq_coul = comm_map2_first(mpi_comm_global_h.comm, coul_libri, s0_s1.first, s0_s1.second);
            double comm_end = omp_get_wtime();
            mpi_comm_global_h.barrier();
            //printf("End vq comm_map2_first  myid: %d   TIME_USED: %f\n",mpi_comm_global_h.myid,comm_end-comm_begin);
            // ofs_myid << "IJq_coul" << endl << IJq_coul;
            //printf("Finish RPA blacs 2d  vq 2d\n");
            double block_begin = omp_get_wtime();
            // for (const auto &IJ: set_IJ_nabf_nabf)
            // {
            //     const auto &I = IJ.first;
            //     const auto &J = IJ.second;
            //     // cout << IJq_coul.at(I).at({J, qa});
            //     collect_block_from_IJ_storage_syhe(
            //         coul_block, desc_nabf_nabf, LIBRPA::atomic_basis_abf, IJ.first,
            //         IJ.second, true, CONE, IJq_coul.at(I).at({J, qa}).ptr(), MAJOR::ROW);
            //     // lib_printf("myid %d I %d J %d nr %d nc %d\n%s",
            //     //        blacs_ctxt_global_h.myid, I, J,
            //     //        coul_block.nr(), coul_block.nc(),
            //     //        str(coul_block).c_str());
            // }
            collect_block_from_ALL_IJ_Tensor(coul_block, desc_nabf_nabf, LIBRPA::atomic_basis_abf,
                                             qa, true, CONE, IJq_coul, MAJOR::ROW);
            double block_end = omp_get_wtime();
            // lib_printf("Vq Time  myid: %d  arr_time: %f  comm_time: %f   block_time: %f   pair_size: %d\n",mpi_comm_global_h.myid,arr_end-vq_begin, comm_end-comm_begin, block_end-block_begin,set_IJ_nabf_nabf.size());
            mpi_comm_global_h.barrier();
            double vq_end = omp_get_wtime();

            if(mpi_comm_global_h.myid == 0)
                lib_printf(" | Total vq time: %f  lri_coul: %f   comm_vq: %f   block_vq: %f\n",vq_end-vq_begin, comm_begin-vq_begin,block_begin-comm_begin,vq_end-block_begin);
        }
    }
    lib_printf("Success test_libcomm_for_system\n");
#endif
}
