#include "cal_periodic_chi0.h"
#include "lapack_connector.h"
#include <omp.h>
#include <iomanip>
#include "parallel_mpi.h"
#include "profiler.h"
#include "meanfield.h"
#include "timefreq.h"
#include "constants.h"
#include "ri.h"
#include "input.h"
#include <iostream>
/*
#include "scalapack_connector.h"
#include "src_lcao/global_fp.h"
#include "src_pw/global.h"
#include "src_lcao/exx_lcao.h"
#include "src_lcao/abfs.h"
#include "src_external/src_test/test_function.h"
#include "Gauss_Quadrature.h"
*/

using namespace std;

void Cal_Periodic_Chi0::chi0_main(const char *Input_ngrid, const char *Input_green_threshold)
{
    if (meanfield.get_n_spins() != 1)
        throw invalid_argument("Unsupported nspins");
    /* print_matrix("wg", meanfield.get_weight()[0]); */
    /* return; */
    prof.start("chi0_main");
    grid_N = stoi(Input_ngrid);
    Green_threshold = stod(Input_green_threshold);
    double emax, emin;
    double gap = meanfield.get_E_min_max(emin, emax);
    /* cout << emin << " " << " " << emax << " " << gap << endl; */
    // TODO: zero division check?
    init(emax/emin);
    cout << "grid_N: " << grid_N << endl;
    double t_begin = omp_get_wtime();
    // check eigenvectors
    // for (int is = 0; is != meanfield.get_n_spins(); is++)
    //     for (int ik = 0; ik != meanfield.get_n_kpoints(); ik++)
    //         print_complex_matrix("wfc:", meanfield.get_eigenvectors()[is][ik]);

    p_id_local = para_mpi.get_myid();
    // LOC.wfc_dm_2d.cal_dm();
    cout << "finish init" << endl;
    R_grid = construct_R_grid();
    cout << "R_grid_num= " << R_grid.size() << endl;
    // cout<< minimax_grid_path+"/time_grid.txt"<<endl;
    time_grid = read_local_grid(grid_N, "local_" + to_string(grid_N) + "_time_points.dat", 'T', gap);
    freq_grid = read_local_grid(grid_N, "local_" + to_string(grid_N) + "_freq_points.dat", 'F', gap);
    tran_gamma = read_cosine_trans_grid(to_string(grid_N) + "_time2freq_grid_cos.txt");
    // time_grid=read_file_grid( minimax_grid_path+"/time_grid.txt",'T');
    // freq_grid=read_file_grid( minimax_grid_path+"/freq_grid.txt",'F');
    // tran_gamma=read_cosine_trans_grid(minimax_grid_path+"/time2freq_transform_grid.txt");
    first_freq = freq_grid.begin()->first;
    first_tau = time_grid.begin()->first;
    n_tau = time_grid.size();

    tau_vec.resize(n_tau);
    int itt = 0;
    for (auto tp : time_grid)
    {
        tau_vec[itt] = tp.first;
        cout << "   itt: " << itt << "   tau_vec: " << tau_vec[itt] << endl;
        itt++;
    }
    freq_vec.resize(n_tau);
    int itf = 0;
    for (auto fp : freq_grid)
    {
        freq_vec[itf] = fp.first;
        itf++;
    }
    // cout << " out_wfc_k" << endl;
    //  for (int ik = 0; ik != n_kpoints; ik++)
    //      print_complex_matrix("wfc_k", wfc_k[ik]);

    for (const auto &time_tau : time_grid)
    {
        for (const Vector3_Order<int> &R : R_grid)
        {
            cal_Green_func_R_tau(time_tau.first, R, kvec_c);
            cal_Green_func_R_tau(-1 * time_tau.first, R, kvec_c);
        }
    }
    cout << "Size of Green_atom map: " << Green_atom.size() << endl;

    /* printf("spin 0, first tau GF: %f\n", first_tau); */
    /* for (auto &R : R_grid) */
    /* { */
    /*     for (int I = 0; I != natom; I++) */
    /*         for (int J = 0; J != natom; J++) */
    /*         { */
    /*             if ( I == J && I == 0 && R == Vector3_Order<int>{-1, -1, -1}) */
    /*             { */
    /*                 cout << "Green_atom  IJ R:  " << I << "  " << J << "    " << R << endl; */
    /*                 print_matrix("unocc Green-atom", Green_atom[0][I][J][R][first_tau]); */
    /*                 print_matrix("occ Green-atom", Green_atom[0][I][J][R][-first_tau]); */
    /*             } */
    /*         } */
    /* } */

    cout << " FINISH GREEN FUN" << endl;
    cout << " Green threshold: " << Green_threshold << endl;
    cout << " Green save num: " << green_save << endl;
    cout << " Green discard num: " << green_discard << endl;
    /* return; */

    // for (auto &tau_p : time_grid)
    // {
    //     for (auto &R : R_grid)
    //     {
    //         for (auto &I_p : Vq)
    //         {
    //             const size_t I = I_p.first;
    //             const size_t mu_num = atom_mu[I];
    //             for (auto &J_p : Vq)
    //             {
    //                 const size_t J = J_p.first;
    //                 const size_t nu_num = atom_mu[J];
    //                 chi0[tau_p.first][R][I][J].create(mu_num, nu_num);
    //             }
    //         }
    //     }
    // }

    vector<pair<atom_t, atom_t>> tot_pair_all(get_atom_pair(Vq));
    int ntask_ap = tot_pair_all.size();
    int ntask_R_tau = R_grid.size() * time_grid.size();
    cout << "  ntask_ap: " << ntask_ap << "    ntask_R_tau: " << ntask_R_tau << endl;
    if (ntask_R_tau > ntask_ap)
    {
        parall_type = "R_tau";
        R_tau_routing();
    }
    else
        atom_pair_routing();
    // atom_pair_routing();

    double t_end = omp_get_wtime();
    cout << "TOTAL TIME USED :  " << t_end - t_begin << endl;
    prof.stop("chi0_main");
}

void Cal_Periodic_Chi0::R_tau_routing()
{
    prof.start("R_tau_routing");
    cout << "Go R_tau_routing!!!    tau_vec: " << tau_vec[0] << endl;
    vector<vector<pair<int, Vector3_Order<int>>>> p_task(process_task_tau_R(time_grid, R_grid));
    map<double, map<Vector3_Order<double>, map<size_t, map<size_t, ComplexMatrix>>>> tmp_chi0_freq_k;
    for (auto &freq_p : freq_grid)
    {
        for (auto &k_p : irk_weight)
        {
            for (auto &I_p : Vq)
            {
                const size_t I = I_p.first;
                const size_t mu_num = atom_mu[I];
                for (auto &J_p : I_p.second)
                {
                    const size_t J = J_p.first;
                    const size_t nu_num = atom_mu[J];
                    tmp_chi0_freq_k[freq_p.first][k_p.first][I][J].create(mu_num, nu_num);
                }
            }
        }
    }
    /* cout << "tmp_chi0_freq_k built" << endl; */

    omp_lock_t chi0_lock;
    omp_init_lock(&chi0_lock);
    double t_chi0_begin = omp_get_wtime();
    double t_chi0_tot = 0;
#pragma omp parallel for schedule(dynamic)
    for (size_t tau_R_index = 0; tau_R_index != p_task[p_id_local].size(); tau_R_index++)
    {
        auto itt = p_task[p_id_local][tau_R_index].first;
        const double tau_loc = tau_vec[itt];
        auto R = p_task[p_id_local][tau_R_index].second;
        printf("   itt:  %d,  tau_loc:  %f , R: ( %d, %d , %d  )\n",itt,tau_loc,R.x,R.y,R.z);
        double t_Rtau_begin = omp_get_wtime();
        for (auto &I_p : Vq)
        {
            const size_t I = I_p.first;
            for (auto &J_p : I_p.second)
            {
                const size_t J = J_p.first;
                double chi0_ele_begin = omp_get_wtime();
                ComplexMatrix tmp_chi0_tau(ComplexMatrix(cal_chi0_element(tau_loc, R, I, J)));
                // chi0[tau_loc][R][I][J] = tmp_chi0_tau.real();
                double chi0_ele_t = omp_get_wtime() - chi0_ele_begin;

                /* printf("      thread: %d, tau: %f,   I: %d  J: %d   , TIME_part: %f \n",omp_get_thread_num(),tau_loc, I,J,chi0_ele_t); */
                omp_set_lock(&chi0_lock);
                for (auto &k_pair : irk_weight)
                {
                    auto ik_vec = k_pair.first;
                    const double arg = (ik_vec * (R * latvec)) * TWO_PI;
                    const complex<double> kphase = complex<double>(cos(arg), sin(arg));
                    for (auto ifreq = 0; ifreq != freq_grid.size(); ifreq++)
                    {
                        const complex<double> cos_weight_kpashe = kphase * tran_gamma[ifreq * n_tau + itt];
                        //  cout<<"  tmp_chi0  nr nc:  "<<tmp_chi0_freq_k[ifreq][ik_vec].nr<<"  "<<tmp_chi0_freq_k[ifreq][ik_vec].nc<<"    chi0_tau: "<<tmp_chi0_tau.nr<<"  "<<tmp_chi0_tau.nr<<endl;
                        tmp_chi0_freq_k[freq_vec[ifreq]][ik_vec][I][J] += tmp_chi0_tau * cos_weight_kpashe;
                    }
                }
                omp_unset_lock(&chi0_lock);
            }
        }
        double t_Rtau_end = omp_get_wtime();
        double time_used = t_Rtau_end - t_Rtau_begin;
        t_chi0_tot += time_used;
        printf("CHI0 p_id: %d,   thread: %d,     R:  ( %d,  %d,  %d ),    tau: %f ,   TIME_USED: %f \n ", p_id_local, omp_get_thread_num(), R.x, R.y, R.z, tau_loc, time_used);
    }
    omp_destroy_lock(&chi0_lock);
    double t_chi0_end = omp_get_wtime();
#pragma omp barrier
    //#pragma omp single
    for (auto &freq_p : freq_grid)
    {
        for (auto &k_p : irk_weight)
        {
            for (auto &I_p : Vq)
            {
                const size_t I = I_p.first;
                const size_t mu_num = atom_mu[I];
                for (auto &J_p : I_p.second)
                {
                    const size_t J = J_p.first;
                    const size_t nu_num = atom_mu[J];
                    chi0_k[freq_p.first][k_p.first][I][J].create(mu_num, nu_num);
                    // ComplexMatrix tmp_chi0=std::move(tmp_chi0_freq_k[freq_vec[ifreq]][ik_vec][I][J]);
                    para_mpi.reduce_ComplexMatrix(tmp_chi0_freq_k[freq_p.first][k_p.first][I][J], chi0_k[freq_p.first][k_p.first][I][J]);
                    tmp_chi0_freq_k[freq_p.first][k_p.first][I].erase(I); // should be erase(J)?
                }
            }
        }
    }
    // Cosine_to_chi0_freq(time_grid,freq_grid);
    // FT_R_to_K_chi0();

    cal_pi_k_use_aims_vq();
    RPA_correlation_energy(freq_grid);
    prof.stop("R_tau_routing");
    cout << "ONLY CHI0 TIME USED:  " << t_chi0_end - t_chi0_begin << endl;
    cout << "  Average per_task time:" << t_chi0_tot / (p_task[p_id_local].size());
}

void Cal_Periodic_Chi0::atom_pair_routing()
{
    prof.start("atom_pair_routing");
    cout << "Go atom_pair_routing!!!" << endl;
    vector<pair<size_t, size_t>> tot_pair(process_task_ap_local());
    cout << "  tot_pair.size :  " << tot_pair.size() << endl;
    omp_lock_t chi0_lock;
    omp_init_lock(&chi0_lock);
    double time_task_tot = 0;
    double t_chi0_begin = omp_get_wtime();
#pragma omp parallel
    {
#pragma omp for schedule(dynamic)
        for (size_t atom_pair = 0; atom_pair != tot_pair.size(); atom_pair++)
        {
            auto I = tot_pair[atom_pair].first;
            auto J = tot_pair[atom_pair].second;
            const size_t mu_num = atom_mu[I];
            const size_t nu_num = atom_mu[J];
            double task_begin = omp_get_wtime();
            map<double, map<Vector3_Order<double>, ComplexMatrix>> tmp_chi0_freq_k;
            for (auto &freq_p : freq_grid)
            {
                for (auto &k_p : irk_weight)
                    tmp_chi0_freq_k[freq_p.first][k_p.first].create(mu_num, nu_num);
            }

            for (auto &R : R_grid)
            {
                for (auto it = 0; it != tau_vec.size(); it++)
                {
                    ComplexMatrix tmp_chi0_tau(ComplexMatrix(cal_chi0_element(tau_vec[it], R, I, J)));
                    // chi0[tau_vec[it]][R][I][J] = tmp_chi0_tau.real();
                    //  chi0_tau_tmp[it]=std::move(tmp_chi0);
                    for (auto &k_pair : irk_weight)
                    {
                        auto ik_vec = k_pair.first;
                        const double arg = (ik_vec * (R * latvec)) * TWO_PI;
                        const complex<double> kphase = complex<double>(cos(arg), sin(arg));
                        for (auto ifreq = 0; ifreq != freq_grid.size(); ifreq++)
                        {
                            const complex<double> cos_weight_kpashe = kphase * tran_gamma[ifreq * n_tau + it];
                            //  cout<<"  tmp_chi0  nr nc:  "<<tmp_chi0_freq_k[ifreq][ik_vec].nr<<"  "<<tmp_chi0_freq_k[ifreq][ik_vec].nc<<"    chi0_tau: "<<tmp_chi0_tau.nr<<"  "<<tmp_chi0_tau.nr<<endl;
                            tmp_chi0_freq_k[freq_vec[ifreq]][ik_vec] += tmp_chi0_tau * cos_weight_kpashe;
                        }
                    }
                }
                // vector<matrix> chi0_freq_tmp(tmp_cosine_tran(freq_grid.size(),chi0_tau_tmp));
            }

            double task_end = omp_get_wtime();
            double time_used = task_end - task_begin;
            time_task_tot += time_used;
            omp_set_lock(&chi0_lock);
            for (auto &freq_p : freq_grid)
            {
                for (auto &k_p : irk_weight)
                    chi0_k[freq_p.first][k_p.first][I][J] = std::move(tmp_chi0_freq_k[freq_p.first][k_p.first]);
            }
            double add_end = omp_get_wtime();
            double add_time = add_end - task_end;
            printf("CHI0 p_id: %d,   thread: %d,     I: %zu,   J: %zu,   add time: %f  TIME_USED: %f \n ", p_id_local, omp_get_thread_num(), I, J, add_time, time_used);
            omp_unset_lock(&chi0_lock);
        }
    }
    omp_destroy_lock(&chi0_lock);
#pragma omp barrier
    double t_chi0_end = omp_get_wtime();
#pragma omp single
    cal_pi_k_use_aims_vq();

    RPA_correlation_energy(freq_grid);

    double t_end = omp_get_wtime();
    cout << "ONLY CHI0 TIME USED:  " << t_chi0_end - t_chi0_begin << endl;
    cout << "  Average per_task time:" << time_task_tot / (tot_pair.size());
    prof.stop("atom_pair_routing");
}

vector<vector<pair<int, Vector3_Order<int>>>> Cal_Periodic_Chi0::process_task_tau_R(const map<double, double> &time_grid, const set<Vector3_Order<int>> &R_grid)
{
    int tot_num = time_grid.size() * R_grid.size();
    int p_num = para_mpi.get_size();
    vector<vector<pair<int, Vector3_Order<int>>>> p_task(p_num);
    int itR = 0;
    for (auto &R : R_grid)
    {
        int pid = itR % p_num;
        int itt = 0;
        for (auto &tau : time_grid)
        {
            p_task[pid].push_back(pair<int, Vector3_Order<int>>(itt, R));
            itt++;
        }
        itR++;
    }
    return p_task;
}

vector<pair<size_t, size_t>> Cal_Periodic_Chi0::process_task_ap_local()
{
    vector<pair<size_t, size_t>> tot_pair(get_atom_pair(Vq));
    int tot_num = tot_pair.size();
    int p_num = para_mpi.get_size();
    int p_id = para_mpi.get_myid();
    vector<pair<size_t, size_t>> p_task_ap;
    int itp = 0;
    for (const auto &ap : tot_pair)
    {
        int p_task = itp % p_num;
        if (p_task == p_id)
        {
            // cout<<"  p_id:  "<<p_id<<"   I J :  "<<ap.first<<"   "<<ap.second<<endl;
            p_task_ap.push_back(ap);
        }
        itp++;
    }
    // for(auto &ap:p_task_ap)
    //     printf("task p_id: %d,   I J:   %d   %d \n",p_id,ap.first,ap.second);
    return p_task_ap;
}

vector<pair<size_t, size_t>> Cal_Periodic_Chi0::get_atom_pair(const map<size_t, map<size_t, map<Vector3_Order<double>, std::shared_ptr<ComplexMatrix>>>> &m)
{
    vector<pair<size_t, size_t>> atom_pairs;
    for (const auto &mA : m)
        for (const auto &mB : mA.second)
            atom_pairs.push_back({mA.first, mB.first});
    // printf("  tot atom_pairs tasks: %d\n",atom_pairs.size());
    return atom_pairs;
}

void Cal_Periodic_Chi0::init(double erange)
{
    cout << "Green threshold:  " << Green_threshold << endl;
    // para_mpi.mpi_init(argc,argv);
    int temp_Cs_count = 0;
    cout << "  Tot  Cs_dim   " << Cs.size() << "    " << Cs[0].size() << "   " << Cs[0][0].size() << endl;
    cout << "  Tot  Vq_dim   " << Vq.size() << "    " << Vq[0].size() << "   " << Vq[0][0].size() << endl;
    for (auto &Ip : Cs)
    {
        size_t I = Ip.first;
        for (auto &Jp : Ip.second)
        {
            size_t J = Jp.first;
            cout << "  I J  Cs_R.num " << I << "  " << J << "   " << Cs[I][J].size() << endl;
        }
    }
    cout << "Cs[0][0] R :" << endl;
    for (auto &R : Cs[0][0])
        cout << R.first << endl;
    // char shell_commond[50];
    string tmps;
    tmps = "python " + GX_path + " " + to_string(grid_N) + " " + to_string(erange);
    cout << tmps.c_str() << endl;
    system(tmps.c_str());

    // cout<<" knpoints_num= "<<n_kpoints<<endl;
    // green_k.resize(n_kpoints);
    /*
        wfc_k.resize(n_kpoints);
        //cout<<"begin to init wfc_k"<<endl;
        for(int i=0;i!=wfc_k.size();++i)
        {
            wfc_k[i].create(LOC.wfc_dm_2d.wfc_k[i].nr,LOC.wfc_dm_2d.wfc_k[i].nc);
            //cout<<"nc="<<LOC.wfc_dm_2d.wfc_k[i].nr<<"   nc="<<LOC.wfc_dm_2d.wfc_k[i].nc<<endl;
            //cout<<"after create "<<endl;
            //cout<<"wfc_k[i].size="<<wfc_k[i].size<<"   LOC.wfc_dm_2d.wfc_k[i].size="<<LOC.wfc_dm_2d.wfc_k[i].size<<endl;
            for(int j=0;j!=LOC.wfc_dm_2d.wfc_k[i].size;++j)
            {
                wfc_k[i].c[j]=LOC.wfc_dm_2d.wfc_k[i].c[j];
            }
            cout<<"   ik="<<i<<"  ";
            //print_complex_matrix("  wfc_k_mat ",LOC.wfc_dm_2d.wfc_k[i].nc,LOC.wfc_dm_2d.wfc_k[i].nr,LOC.wfc_dm_2d.wfc_k[i].c,LOC.wfc_dm_2d.wfc_k[i].nc);
        }
        //cout<<"end init wfc_k"<<endl;

        */
}

set<Vector3_Order<int>> Cal_Periodic_Chi0::construct_R_grid()
{
    // cout<<" begin to construct_R_grid"<<endl;
    set<Vector3_Order<int>> R_grid;
    R_grid.clear();

    R_grid.insert({0, 0, 0});
    for (int x = -(kv_nmp[0]) / 2; x <= (kv_nmp[0] - 1) / 2; ++x)
        for (int y = -(kv_nmp[1]) / 2; y <= (kv_nmp[1] - 1) / 2; ++y)
            for (int z = -(kv_nmp[2]) / 2; z <= (kv_nmp[2] - 1) / 2; ++z)
                R_grid.insert({x, y, z});

    R_periodic = Vector3<int>{kv_nmp[0], kv_nmp[1], kv_nmp[2]};
    for (const auto &R : R_grid)
        cout << R << endl;
    cout << " R_periodic: " << R_periodic << endl;
    return R_grid;
}

vector<vector<ComplexMatrix>> Cal_Periodic_Chi0::cal_Green_func_element_k(const vector<matrix> &wg, const double time_tau)
{
    // cout<<wg.nr<<"   "<<wg.nc<<"  "<<time_tau<<endl;
    //  dm = wfc.T * wg * wfc.conj()
    //  dm[ik](iw1,iw2) = \sum_{ib} wfc[ik](ib,iw1).T * wg(ik,ib) * wfc[ik](ib,iw2).conj()
    // assert(wg.nc<=NLOCAL);
    // assert(wg.nr==n_kpoints);
    vector<vector<ComplexMatrix>> green_sk(meanfield.get_n_spins());
    for (int is = 0; is != meanfield.get_n_spins(); is++)
        green_sk[is].resize(meanfield.get_n_kpoints());

    auto efermi = meanfield.get_efermi();
    
    for (int is = 0; is != meanfield.get_n_spins(); is++)
    {
        auto const & wfc_s = meanfield.get_eigenvectors()[is];
        auto const & ekb = meanfield.get_eigenvals()[is];

        for (int ik = 0; ik != meanfield.get_n_kpoints(); ++ik)
        {
            // cout<<" ik= "<<ik<<endl;
            std::vector<double> wg_local(wg[is].nc, 0.0);
            for (int ib_global = 0; ib_global != wg[is].nc; ++ib_global)
            {
                // const int ib_local = ParaO.trace_loc_col[ib_global];
                double temp_e_up = -1 * 0.5 * (ekb(ik, ib_global) - efermi) * time_tau;
                double ephase;
                if (temp_e_up >= 0)
                {
                    ephase = 1.0;
                }
                else
                {
                    ephase = std::exp(temp_e_up);
                }
                wg_local[ib_global] = wg[is](ik, ib_global) * ephase;
            }
            // wg_wfc(ib,iw) = wg[ib] * wfc(ib,iw).conj();
            // cout<<wfc_k.size()<<endl;
            // cout<<"wfc_k  dim "<<wfc_k[ik].nr<<wfc_k[ik].nc<<endl;
            ComplexMatrix wg_wfc = conj(wfc_s[ik]);
            // cout<<" wg_wfc "<<wg_wfc.nr<<"  "<<wg_wfc.nc<<endl;
            for (int ir = 0; ir != wg_wfc.nr; ++ir)
                LapackConnector::scal(wg_wfc.nc, wg_local[ir], wg_wfc.c + ir * wg_wfc.nc, 1);
            // C++: dm(iw1,iw2) = wfc(ib,iw1).T * wg_wfc(ib,iw2)
            green_sk[is][ik].create(wfc_s[ik].nr, wfc_s[ik].nc);
            green_sk[is][ik] = transpose(wfc_s[ik], 0) * wg_wfc;
            /* printf("green_sk is=%d, ik=%d, tau=%f\n", is, ik, time_tau); */
            /* print_complex_matrix("green_sk:", green_sk[is][ik]); */
        }
    }
    return green_sk;
}

void Cal_Periodic_Chi0::cal_Green_func_R_tau(const double &time_tau, const Vector3_Order<int> &R, const Vector3<double> *kvec_c)
{
    prof.start("cal_Green_func_R_tau");
    auto nspins = meanfield.get_n_spins();
    auto nkpts = meanfield.get_n_kpoints();
    auto nbands = meanfield.get_n_bands();
    auto const & wg = meanfield.get_weight();

    Vector3_Order<double> R_2(R.x, R.y, R.z);

    Vector3_Order<int> R_0(0, 0, 0);
    // cout<<"R_0 init"<<endl;
    // cout<<"R= "<<R<<"     time_tau= "<<time_tau<<endl;

    // cout<<"   zero_out";
    // cout<<"Green_function_zero_out"<<endl;
    vector<vector<ComplexMatrix>> green_sk;
    if (time_tau > 0)
    {
        vector<matrix> wg_unocc(nspins);
        for (int is = 0; is != nspins; is++)
        {
            wg_unocc[is].create(nkpts, nbands);
            for (int i = 0; i != wg_unocc[is].nr * wg_unocc[is].nc; ++i)
            {
                // FIXME: check if imat should be 2.0 / nkpts / nspins
                wg_unocc[is].c[i] = 1.0 / nkpts * nspins - wg[is].c[i];
                if (wg_unocc[is].c[i] < 0) wg_unocc[is].c[i] = 0;
            }
        }

        // print_matrix("wg_occ_mat:",wg);
        /* print_matrix("wg_unocc_mat:",wg_unocc[0]); */
        // cout<<"cal green ele"<<endl;
        green_sk = cal_Green_func_element_k(wg_unocc, time_tau);
    }
    else
    {
        green_sk = cal_Green_func_element_k(wg, time_tau);
    }
    /* if ( R == Vector3_Order<int>{-1, -1, -1} && time_tau < 0.012 ) */
    /* { */
    /*     cout << R << " " << time_tau << endl; */
    /*     print_complex_matrix("green_sk[0][0]:", green_sk[0][0]); */
    /* } */
    /* return; */

    for (int is = 0; is != meanfield.get_n_spins(); is++)
    {
        ComplexMatrix Green_glo_tmp(green_sk[is][0].nr, green_sk[is][0].nc);
        Green_glo_tmp.zero_out();
        for (int ik = 0; ik != meanfield.get_n_kpoints(); ++ik)
        {
            /* cout << "ik: " << ik << endl; */
            const double arg = -1 * (kvec_c[ik] * (R_2 * latvec)) * TWO_PI; // latvec
            const complex<double> kphase = complex<double>(cos(arg), sin(arg));
            Green_glo_tmp += green_sk[is][ik] * kphase;
        }
        /* if (time_tau == first_tau) */
        /*     print_complex_matrix("Green_glo_tmp:", Green_glo_tmp); */

        if (time_tau <= 0)
        {
            Green_glo_tmp *= -1;
        }
        for (int I = 0; I != natom; I++)
        {
            const size_t I_num = atom_nw[I];
            for (int J = 0; J != natom; J++)
            {
                const size_t J_num = atom_nw[J];
                // Green_atom[is][I][J][time_tau][R].create(I_num,J_num);
                matrix tmp_green(I_num, J_num);
                for (size_t i = 0; i != I_num; i++)
                {
                    size_t i_glo = atom_iw_loc2glo(I, i);
                    for (size_t j = 0; j != J_num; j++)
                    {
                        size_t j_glo = atom_iw_loc2glo(J, j);
                        // Green_atom[is][I][J][time_tau][R](i,j)=Green_glo_tmp(i_glo,j_glo).real();
                        tmp_green(i, j) = Green_glo_tmp(i_glo, j_glo).real();
                    }
                }
                if (tmp_green.absmax() > Green_threshold)
                {
                    // cout<<" max_green_ele:  "<<tmp_green.absmax()<<endl;
                    Green_atom[is][I][J][R][time_tau] = std::move(tmp_green);
                    green_save++;
                }
                else
                {
                    green_discard++;
                }
            }
        }
    }
    /* cout << "Done Green R" << endl; */
    prof.stop("cal_Green_func_R_tau");
}

void Cal_Periodic_Chi0::Cosine_to_chi0_freq(map<double, double> &time_grid, map<double, double> &freq_grid)
{
    // map<double,double> freq_grid(construct_gauss_grid(120));
    int n_freq = freq_grid.size();
    int n_tau = time_grid.size();
    assert(n_freq == n_tau);
    vector<double> tran_gamma = read_cosine_trans_grid("/time2freq_transform_grid.txt");
    for (const auto &time_pair : chi0)
    {
        // const double time=time_pair.first;
        const auto chi0_time = time_pair.second;
        for (const auto &R_pair : chi0_time)
        {
            const auto R = R_pair.first;
            const auto chi0_time_R = R_pair.second;
            for (const auto &I_pair : chi0_time_R)
            {
                const size_t I = I_pair.first;
                const auto chi0_time_R_I = I_pair.second;
                for (const auto &J_pair : chi0_time_R_I)
                {
                    const size_t J = J_pair.first;
                    const auto mat = J_pair.second;
                    for (const auto &freq_pair : freq_grid)
                    {
                        const double freq = freq_pair.first;
                        const double freq_weight = freq_pair.second;
                        chi0_freq[freq][R][I][J].create(mat.nr, mat.nc);
                    }
                }
            }
        }
        break;
    }

    int iter_tau = 0;
    for (const auto &time_pair : chi0)
    {
        const double time = time_pair.first;
        const auto chi0_time = time_pair.second;
        for (const auto &R_pair : chi0_time)
        {
            const auto R = R_pair.first;
            const auto chi0_time_R = R_pair.second;
            for (const auto &I_pair : chi0_time_R)
            {
                const size_t I = I_pair.first;
                const auto chi0_time_R_I = I_pair.second;
                for (const auto &J_pair : chi0_time_R_I)
                {
                    const size_t J = J_pair.first;
                    const auto mat = J_pair.second;
                    int iter_freq = 0;
                    for (const auto &freq_pair : freq_grid)
                    {

                        const double freq = freq_pair.first;
                        const double cos_weight = tran_gamma[iter_freq * n_tau + iter_tau];
                        chi0_freq.at(freq).at(R).at(I).at(J) += mat * cos_weight;
                        //cout<<time<<"   "<<time_weight<<"    freq:"<<freq<<"   cos_weight"<<cos_weight  \
                            <<"      mat:  "<<mat(0,0)<<"    "<<chi0_freq[freq][R][I][J](0,0)<<endl;
                        iter_freq++;
                    }
                }
            }
        }
        iter_tau++;
    }
    cout << "finish cosine tran chi0 !" << endl;
}

// vector<matrix> tmp_cosine_tran(double n_freq, vector<matrix> &tau_mat_vec)
// {
//     int n_tau=tau_mat_vec.size();
//     vector<matrix> tmp_chi0_freq_vec(freq_grid.size())
//     for(auto it=0;it!=n_tau;it++)
//     {
//         for(auto ifreq=0;ifreq!=tmp_chi0_freq_vec.size();ifreq++)
//         {
//             const double cos_weight=tran_gamma[ifreq*n_tau+it];
//             tmp_chi0_freq_vec[ifreq]+=tau_mat_vec[it]*cos_weight;
//         }
//     }
//     return tmp_chi0_freq_vec;
// }

matrix Cal_Periodic_Chi0::cal_chi0_element(const double &time_tau, const Vector3_Order<int> &R, const size_t &I_index, const size_t &J_index)
{
    prof.start("cal_chi0_element");
    /* printf("     begin chi0  thread: %d,  I: %zu, J: %zu\n",omp_get_thread_num(),I_index,J_index); */
    // for(const auto &K_pair:Cs[I_index])
    // Vector3_Order<int> R_0(0,0,0);
    int flag_G_IJRt = 0;
    int flag_G_IJRNt = 0;

    /* printf("     check if already calculated\n"); */
    /* printf("     size of Green_atom: %zu\n", Green_atom.size()); */
    Green_atom.at(0);
    if (Green_atom.at(0).at(I_index).count(J_index))
        if (Green_atom.at(0).at(I_index).at(J_index).count(R))
        {
            if (Green_atom.at(0).at(I_index).at(J_index).at(R).count(time_tau))
                flag_G_IJRt = 1;

            if (Green_atom.at(0).at(I_index).at(J_index).at(R).count(-time_tau))
                flag_G_IJRNt = 1;
        }
    /* printf("     checked\n"); */

    const size_t i_num = atom_nw[I_index];
    const size_t mu_num = atom_mu[I_index];

    const size_t j_num = atom_nw[J_index];
    const size_t nu_num = atom_mu[J_index];

    // chi0[time_tau][R][I_index][J_index].create(mu_num,nu_num);
    matrix X_R2(i_num, j_num * nu_num);
    matrix X_conj_R2(i_num, j_num * nu_num);
    
    /* printf("     begin L_pair loop\n"); */
    for (const auto &L_pair : Cs[J_index])
    {
        const auto L_index = L_pair.first;
        const size_t l_num = atom_nw[L_index];
        /* printf("     begin is loop\n"); */
        for (int is = 0; is != meanfield.get_n_spins(); is++)
        {
            if (Green_atom.at(is).at(I_index).count(L_index))
            {
                for (const auto &R2_index : L_pair.second)
                {
                    const auto R2 = R2_index.first;
                    const auto &Cs_mat2 = R2_index.second;

                    Vector3_Order<int> R_temp_2(Vector3_Order<int>(R2 + R) % R_periodic);

                    if (Green_atom.at(is).at(I_index).at(L_index).count(R_temp_2))
                    {
                        assert(j_num * l_num == (*Cs_mat2).nr);
                        /* printf("          thread: %d, X_R2 IJL:   %zu,%zu,%zu  R:(  %d,%d,%d  )  tau:%f\n",omp_get_thread_num(),I_index,J_index,L_index,R.x,R.y,R.z,time_tau); */
                        matrix Cs2_reshape(reshape_Cs(j_num, l_num, nu_num, Cs_mat2));
                        
                        if (Green_atom.at(is).at(I_index).at(L_index).at(R_temp_2).count(time_tau))
                        {
                            // cout<<"C";
                            prof.start("X");
                            X_R2 += Green_atom.at(is).at(I_index).at(L_index).at(R_temp_2).at(time_tau) * Cs2_reshape;
                            prof.stop("X");
                        }

                        if (Green_atom.at(is).at(I_index).at(L_index).at(R_temp_2).count(-time_tau))
                        {
                            // cout<<"D";
                            prof.start("X");
                            X_conj_R2 += Green_atom.at(is).at(I_index).at(L_index).at(R_temp_2).at(-time_tau) * Cs2_reshape;
                            prof.stop("X");
                        }
                        
                    }
                }
            }
        }
    }
    matrix X_R2_rs(reshape_mat(i_num, j_num, nu_num, X_R2));
    matrix X_conj_R2_rs(reshape_mat(i_num, j_num, nu_num, X_conj_R2));
    

    matrix O_sum(mu_num, nu_num);
    for (const auto &K_pair : Cs[I_index])
    {
        const auto K_index = K_pair.first;
        const size_t k_num = atom_nw[K_index];

        for (const auto &R1_index : K_pair.second)
        {
            const auto R1 = R1_index.first;
            const auto &Cs_mat1 = R1_index.second;
            // cout<<"R1:  begin  "<<R1<<endl;
            for (int is = 0; is != meanfield.get_n_spins(); is++)
            {
                // matrix X_R2(i_num,j_num*nu_num);
                // matrix X_conj_R2(i_num,j_num*nu_num);
                matrix O(i_num, k_num * nu_num);
                matrix Z(k_num, i_num * nu_num);

                if (flag_G_IJRt || flag_G_IJRNt)
                {
                    matrix N_R2(k_num, j_num * nu_num);
                    matrix N_conj_R2(k_num, j_num * nu_num);
                    for (const auto &L_pair : Cs[J_index])
                    {
                        const auto L_index = L_pair.first;
                        const size_t l_num = atom_nw[L_index];
                        if (Green_atom.at(is).at(K_index).count(L_index))
                        {
                            for (const auto &R2_index : L_pair.second)
                            {
                                const auto R2 = R2_index.first;
                                const auto &Cs_mat2 = R2_index.second;
                                Vector3_Order<int> R_temp_1(Vector3_Order<int>(R + R2 - R1) % R_periodic);
                                // Vector3_Order<int> R_temp_2(R2+R);
                                if (Green_atom.at(is).at(K_index).at(L_index).count(R_temp_1))
                                {
                                    
                                    assert(j_num * l_num == (*Cs_mat2).nr);
                                    // printf("          thread: %d, IJKL:   %d,%d,%d,%d  R:(  %d,%d,%d  )  tau:%f\n",omp_get_thread_num(),I_index,J_index,K_index,L_index,R.x,R.y,R.z,time_tau);
                                    matrix Cs2_reshape(reshape_Cs(j_num, l_num, nu_num, Cs_mat2));
                                    if (flag_G_IJRNt && Green_atom.at(is).at(K_index).at(L_index).at(R_temp_1).count(time_tau))
                                    {
                                        // cout<<"A";
                                        prof.start("N");
                                        N_R2 += Green_atom.at(is).at(K_index).at(L_index).at(R_temp_1).at(time_tau) * Cs2_reshape;
                                        prof.stop("N");
                                    }
                                    if (flag_G_IJRt && Green_atom.at(is).at(K_index).at(L_index).at(R_temp_1).count(-time_tau))
                                    {
                                        // cout<<"B";
                                        prof.start("N");
                                        N_conj_R2 += Green_atom.at(is).at(K_index).at(L_index).at(R_temp_1).at(-time_tau) * Cs2_reshape;
                                        prof.stop("N");
                                    }
                                    
                                }
                            }
                        }
                    }
                    if (flag_G_IJRt)
                    {
                        prof.start("O");
                        matrix N_conj_R2_rs(reshape_mat(k_num, j_num, nu_num, N_conj_R2));
                        O += Green_atom.at(is).at(I_index).at(J_index).at(R).at(time_tau) * N_conj_R2_rs;
                        prof.stop("O");
                    }
                    if (flag_G_IJRNt)
                    {
                        prof.start("O");
                        matrix N_R2_rs(reshape_mat(k_num, j_num, nu_num, N_R2));
                        O += Green_atom.at(is).at(I_index).at(J_index).at(R).at(-time_tau) * N_R2_rs;
                        prof.stop("O");
                    }
                }
                // printf("          thread: %d, 2\n",omp_get_thread_num());

                // matrix X_R2_rs(reshape_mat(i_num,j_num,nu_num,X_R2));
                // matrix X_conj_R2_rs(reshape_mat(i_num,j_num,nu_num,X_conj_R2));

                // cal and sum O

                // printf("          thread: %d, 3\n",omp_get_thread_num());
                //  if(Green_atom.at(is).at(I_index).count(J_index))
                //  {
                //      if(Green_atom.at(is).at(I_index).at(J_index).count(-time_tau))
                //      {
                //          if(Green_atom.at(is).at(I_index).at(J_index).at(-time_tau).count(R))
                //          {
                //              O+=Green_atom.at(is).at(I_index).at(J_index).at(-time_tau).at(R)*N_R2_rs;
                //          }
                //      }
                //      if(Green_atom.at(is).at(I_index).at(J_index).count(time_tau))
                //      {
                //          if(Green_atom.at(is).at(I_index).at(J_index).at(time_tau).count(R))
                //          {
                //              O+=Green_atom.at(is).at(I_index).at(J_index).at(time_tau).at(R)*N_conj_R2_rs;
                //          }
                //      }
                //  }
                Vector3_Order<int> R_temp_3(Vector3_Order<int>(R - R1) % R_periodic);
                if (Green_atom.at(is).at(K_index).count(J_index))
                {
                    if (Green_atom.at(is).at(K_index).at(J_index).count(R_temp_3))
                    {
                        if (Green_atom.at(is).at(K_index).at(J_index).at(R_temp_3).count(-time_tau))
                        {
                            prof.start("Z");
                            Z += Green_atom.at(is).at(K_index).at(J_index).at(R_temp_3).at(-time_tau) * X_R2_rs;
                            prof.stop("Z");
                        }

                        if (Green_atom.at(is).at(K_index).at(J_index).at(R_temp_3).count(time_tau))
                        {
                            prof.start("Z");
                            Z += Green_atom.at(is).at(K_index).at(J_index).at(R_temp_3).at(time_tau) * X_conj_R2_rs;
                            prof.stop("Z");
                        }
                    }
                }

                matrix Z_rs(reshape_mat(k_num, i_num, nu_num, Z));

                O += Z_rs;
                matrix OZ(reshape_mat_21(i_num, k_num, nu_num, O));
                matrix Cs1_tran(transpose(*Cs_mat1));
                prof.start("O");
                O_sum += Cs1_tran * OZ;
                prof.stop("O");
                // cout<<"   K, R1:   "<<K_index<<"   "<<R1;
                // rt_m_max(O_sum);
            }
        }
    }
    prof.stop("cal_chi0_element");
    // chi0[time_tau][R][I_index][J_index]=O_sum;
    // printf("     finish chi0  thread: %d,  I: %d, J: %d",omp_get_thread_num(),I_index,J_index);
    return O_sum;
}

matrix Cal_Periodic_Chi0::reshape_Cs(const size_t n1, const size_t n2, const size_t n3, const shared_ptr<matrix> &Cs) //(n1*n2,n3) -> (n2,n1*n3)
{
    prof.start("reshape_Cs");
    const auto length = sizeof(double) * n3;
    const auto n13 = n1 * n3;
    const double *m_ptr = (*Cs).c;
    matrix m_new(n2, n1 * n3, false);
    for (size_t i1 = 0; i1 != n1; ++i1)
    {
        double *m_new_ptr = m_new.c + i1 * n3;
        for (size_t i2 = 0; i2 != n2; ++i2, m_ptr += n3, m_new_ptr += n13)
            memcpy(m_new_ptr, m_ptr, length);
    }
    prof.stop("reshape_Cs");
    return m_new;
}

matrix Cal_Periodic_Chi0::reshape_dim_Cs(const size_t n1, const size_t n2, const size_t n3, const shared_ptr<matrix> &Cs) //(n1*n2,n3) -> (n1,n2*n3)
{
    const auto length = sizeof(double) * n1 * n2 * n3;
    const auto n13 = n1 * n3;
    const double *m_ptr = (*Cs).c;
    matrix m_new(n1, n2 * n3, false);
    double *m_new_ptr = m_new.c;
    memcpy(m_new_ptr, m_ptr, length);
    return m_new;
}

matrix Cal_Periodic_Chi0::reshape_mat(const size_t n1, const size_t n2, const size_t n3, const matrix &mat) //(n1,n2*n3) -> (n2,n1*n3)
{
    prof.start("reshape_mat");
    const auto length = sizeof(double) * n3;
    const auto n13 = n1 * n3;
    const double *m_ptr = mat.c;
    matrix m_new(n2, n1 * n3, false);
    for (size_t i1 = 0; i1 != n1; ++i1)
    {
        double *m_new_ptr = m_new.c + i1 * n3;
        for (size_t i2 = 0; i2 != n2; ++i2, m_ptr += n3, m_new_ptr += n13)
            memcpy(m_new_ptr, m_ptr, length);
    }
    prof.stop("reshape_mat");
    return m_new;
}
matrix Cal_Periodic_Chi0::reshape_mat_21(const size_t n1, const size_t n2, const size_t n3, const matrix &mat) //(n1,n2*n3) -> (n1*n2,n3)
{
    const auto length = sizeof(double) * n1 * n2 * n3;
    const auto n13 = n1 * n2 * n3;
    const double *m_ptr = mat.c;
    matrix m_new(n1 * n2, n3, false);
    double *m_new_ptr = m_new.c;
    memcpy(m_new_ptr, m_ptr, length);
    return m_new;
}

std::map<double, map<Vector3_Order<double>, map<size_t, map<size_t, ComplexMatrix>>>> Cal_Periodic_Chi0::cal_pi_k_use_aims_vq()
{
    // cout << "out_chi0_R" << endl;
    // for (auto &Rp : chi0.at(first_tau))
    // {
    //     auto R = Rp.first;
    //     for (auto &Ip : Rp.second)
    //     {
    //         auto I = Ip.first;
    //         for (auto &Jp : Ip.second)
    //         {
    //             auto J = Jp.first;
    //             cout << "  I J R:  " << I << "  " << J << "   " << R << endl;
    //             print_matrix("chi0_R", Jp.second);
    //         }
    //     }
    // }
    printf("Begin cal_pi_k , pid:  %d\n", para_mpi.get_myid());
    for (auto &freq_p : chi0_k)
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
                const size_t J_mu = atom_mu[J];
                for (auto &Q_p : J_p.second)
                {
                    const size_t Q = Q_p.first;
                    const size_t Q_mu = atom_mu[Q];
                    auto &chi0_mat = Q_p.second;
                    for (auto &I_p : Vq)
                    {
                        const size_t I = I_p.first;
                        const size_t I_mu = atom_mu[I];
                        pi_k[freq][ik_vec][I][Q].create(I_mu, Q_mu);
                        if (J != Q)
                            pi_k[freq][ik_vec][I][J].create(I_mu, J_mu);
                    }
                }
            }
        }
    }

    for (auto &freq_p : chi0_k)
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
                    for (auto &I_p : Vq)
                    {
                        const size_t I = I_p.first;
                        // printf("cal_pi  pid: %d , IJQ:  %d  %d  %d\n", para_mpi.get_myid(), I, J, Q);
                        //  cout<<"         pi_IQ: "<<pi_k.at(freq).at(ik_vec).at(I).at(Q)(0,0)<<"   pi_IJ: "<<pi_k.at(freq).at(ik_vec).at(I).at(J)(0,0);
                        if (I <= J)
                        {
                            // cout << "   1"
                            //      << "  Vq: " << (*Vq.at(I).at(J).at(ik_vec))(0, 0) << endl;
                            pi_k.at(freq).at(ik_vec).at(I).at(Q) += (*Vq.at(I).at(J).at(ik_vec)) * chi0_mat;
                            // if (freq == first_freq)
                            // {
                            //     complex<double> trace_pi;
                            //     trace_pi = trace(pi_k.at(freq).at(ik_vec).at(I).at(Q));
                            //     cout << " IJQ: " << I << " " << J << " " << Q << "  ik_vec: " << ik_vec << "  trace_pi:  " << trace_pi << endl;
                            //     print_complex_matrix("vq:", (*Vq.at(I).at(J).at(ik_vec)));
                            //     print_complex_matrix("chi0:", chi0_mat);
                            //     print_complex_matrix("pi_mat:", pi_k.at(freq).at(ik_vec).at(I).at(Q));
                            // }
                        }
                        else
                        {
                            // cout << "   2"
                            //      << "  Vq: " << transpose(*Vq.at(J).at(I).at(ik_vec), 1)(0, 0) << endl;
                            pi_k.at(freq).at(ik_vec).at(I).at(Q) += transpose(*Vq.at(J).at(I).at(ik_vec), 1) * chi0_mat;
                        }

                        if (J != Q)
                        {
                            ComplexMatrix chi0_QJ = transpose(chi0_mat, 1);
                            if (I <= Q)
                            {
                                // cout << "   3"
                                //      << "  Vq: " << (*Vq.at(I).at(Q).at(ik_vec))(0, 0) << endl;
                                pi_k.at(freq).at(ik_vec).at(I).at(J) += (*Vq.at(I).at(Q).at(ik_vec)) * chi0_QJ;
                            }
                            else
                            {
                                // cout << "   4"
                                //      << "  Vq: " << transpose(*Vq.at(J).at(I).at(ik_vec), 1)(0, 0) << endl;
                                pi_k.at(freq).at(ik_vec).at(I).at(J) += transpose(*Vq.at(Q).at(I).at(ik_vec), 1) * chi0_QJ;
                            }
                        }
                    }
                }
            }
        }
    }

    //         for(auto &I_p:Vq)
    //         {
    //             const auto I=I_p.first;
    //             for(auto &J_p:I_p.second)
    //             {
    //                 const auto J=J_p.first;
    //                 if(chi0_freq_k.count(J))
    //                 {
    //                     for(auto &Q_p:chi0_freq_k.at(J))
    //                     {
    //                         const auto Q=Q_p.first;
    //                         auto Vq_mat=*Q_p.second.at(ik_vec);
    //                         if(freq==first_freq)
    //                         {
    //                             //cout<<"p_id:  "<<para_mpi.get_myid()<<"  IQJ:  "<<I<<"   "<<Q<<"  "<<J<<"   "<<ik_vec;
    //                             printf("pid:  %d    IQJ:  %d  %d  %d\n",para_mpi.get_myid(),I,J,Q);
    //                             //rt_cm_max(Vq_mat);
    //                             // print_complex_matrix("vq",Vq_mat);
    //                             // if(J<=Q)
    //                             // {
    //                             //     print_complex_matrix("chi0_mat",chi0_k[freq][ik_vec][J][Q]);
    //                             //   //  print_complex_matrix("pi_mat",Vq_mat*chi0_k[freq][ik_vec][J][Q]);
    //                             // }
    //                             // else
    //                             // {
    //                             //     print_complex_matrix("chi0_mat",transpose(chi0_k[freq][ik_vec][Q][J],1));
    //                             //     //print_complex_matrix("pi_mat",Vq_mat*transpose(chi0_k[freq][ik_vec][Q][J],1));
    //                             // }
    //                         }
    //                         if(J<=Q)
    //                         {
    //                             //cout<<"   1";
    //                             pi_k.at(freq).at(ik_vec).at(I).at(Q)+=Vq_mat*chi0_k.at(freq).at(ik_vec).at(J).at(Q);
    //                             // if(freq==first_freq)
    //                             // {
    //                             //     // cout<<" pi ";
    //                             //     // rt_cm_max(pi_k.at(freq).at(ik_vec).at(I).at(Q));
    //                             //     // cout<<" chi0 ";
    //                             //     // rt_cm_max(chi0_k.at(freq).at(ik_vec).at(J).at(Q));
    //                             //     // cout<<" Vq ";
    //                             //     // rt_cm_max(Vq_mat);
    //                             // }
    //                             // if(freq==first_freq && I==natom-1 && Q==natom-1)
    //                             // {
    //                             //     cout<<"IQJ: "<<I<<"   "<<Q<<"   "<<J<<"    pi_ele: "<<pi_k.at(freq).at(ik_vec).at(I).at(Q)(0,0)<<endl;
    //                             //     rt_cm_max(chi0_k.at(freq).at(ik_vec).at(J).at(Q));
    //                             //     rt_cm_max(Vq_mat);
    //                             // }

    //                         }
    //                         else
    //                         {
    //                             //cout<<"   2";
    //                             pi_k.at(freq).at(ik_vec).at(I).at(Q)+=Vq_mat*transpose(chi0_k.at(freq).at(ik_vec).at(Q).at(J),1);
    //                             // if(freq==first_freq && I==natom-1 && Q==natom-1)
    //                             // {
    //                             //     cout<<"IQJ: "<<I<<"   "<<Q<<"   "<<J<<"    pi_ele: "<<pi_k.at(freq).at(ik_vec).at(I).at(Q)(0,0)<<endl;
    //                             //     rt_cm_max(chi0_k.at(freq).at(ik_vec).at(Q).at(J));
    //                             //     rt_cm_max(Vq_mat);
    //                             // }
    //                         }

    //                         if(I!=J)
    //                         {
    //                             if(I<=Q)
    //                             {
    //                                // cout<<"   3"<<"  pi: "<<pi_k[freq][ik_vec][J][Q].nr<<pi_k[freq][ik_vec][J][Q].nc<<"  vq nr nc: "<<(transpose(Vq_mat,1)).nr<<"  "<<(transpose(Vq_mat,1)).nc<<"   chi0:"<<chi0_k[freq][ik_vec][I][Q].nr<<"  "<<chi0_k[freq][ik_vec][I][Q].nc<<endl;
    //                                 pi_k.at(freq).at(ik_vec).at(J).at(Q)+=transpose(Vq_mat,1)*chi0_k.at(freq).at(ik_vec).at(I).at(Q);
    //                                 // if(freq==first_freq && J==natom-1 && Q==natom-1)
    //                                 // {
    //                                 //     cout<<"IQJ: "<<I<<"   "<<Q<<"   "<<J<<"    pi_ele: "<<pi_k.at(freq).at(ik_vec).at(J).at(Q)(0,0)<<endl;
    //                                 //     rt_cm_max(chi0_k.at(freq).at(ik_vec).at(I).at(Q));
    //                                 //     rt_cm_max(Vq_mat);
    //                                 // }
    //                             }
    //                             else
    //                             {
    //                                 //cout<<"   4";
    //                                 pi_k.at(freq).at(ik_vec).at(J).at(Q)+=transpose(Vq_mat,1)*transpose(chi0_k.at(freq).at(ik_vec).at(Q).at(I),1);
    //                                 // if(freq==first_freq && J==natom-1 && Q==natom-1)
    //                                 // {
    //                                 //     cout<<"IQJ: "<<I<<"   "<<Q<<"   "<<J<<"    pi_ele: "<<pi_k.at(freq).at(ik_vec).at(J).at(Q)(0,0)<<endl;
    //                                 //     rt_cm_max(chi0_k.at(freq).at(ik_vec).at(Q).at(I));
    //                                 //     rt_cm_max(Vq_mat);
    //                                 // }
    //                             }
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }
    //  print_complex_matrix("  last_pi_mat:",pi_k.at(first_freq).at({0,0,0}).at(natom-1).at(natom-1));
    printf("finish cal_pi_k , pid:  %d\n", para_mpi.get_myid());
    return pi_k;
}

void Cal_Periodic_Chi0::FT_R_to_K_chi0()
{
    cout << " bigin to FT_R_to_K_chi0" << endl;
    for (auto &freq_pair : chi0_freq)
    {
        const double freq = freq_pair.first;
        const auto chi0_freq = freq_pair.second;
        for (auto &R_pair : chi0_freq)
        {
            const Vector3_Order<int> R = R_pair.first;
            auto chi0_freq_R = R_pair.second;
            for (auto &I_pair : chi0_freq_R)
            {
                const size_t I = I_pair.first;
                auto chi0_freq_R_I = I_pair.second;
                // const size_t mu_num=exx_lcao.index_abfs[ucell.iat2it[I]].count_size;
                for (auto &J_pair : chi0_freq_R_I)
                {
                    size_t J = J_pair.first;
                    // const auto chi0_mat=J_pair.second;
                    //(*Vps[I][J][R]).c
                    // cout<<"2"<<endl;
                    // const size_t nu_num=exx_lcao.index_abfs[ucell.iat2it[J]].count_size;
                    for (auto &q_pair : irk_weight)
                    {
                        // cout<<"3"<<endl;
                        Vector3_Order<double> ik_vec = q_pair.first;
                        // auto &vq_mat=*q_pair.second;
                        chi0_k[freq][ik_vec][I][J].create(atom_mu[I], atom_mu[J]);
                        // cout<<freq<<"    "<<I<<"   "<<J<<"    ik:"<<ik<<"   ik_vec "<<ik_vec<<endl;
                    }
                }
            }
            break;
        }
    }
    cout << "  finish init chi0_k" << endl;
    temp_chi0_k_count = 0;
    for (auto &freq_pair : chi0_freq)
    {
        const double freq = freq_pair.first;
        const auto chi0_freq = freq_pair.second;

        for (auto &R_pair : chi0_freq)
        {
            const Vector3_Order<int> R = R_pair.first;
            auto chi0_freq_R = R_pair.second;
            for (auto &I_pair : chi0_freq_R)
            {
                const size_t I = I_pair.first;
                auto chi0_freq_R_I = I_pair.second;
                const size_t mu_num = atom_mu[I];
                for (auto &J_pair : chi0_freq_R_I)
                {
                    size_t J = J_pair.first;
                    const size_t nu_num = atom_mu[J];
                    auto &chi0_mat = J_pair.second;
                    size_t mat_size = chi0_mat.nr * chi0_mat.nc;
                    for (auto &ikp : chi0_k[freq])
                    {
                        auto ik_vec = ikp.first;
                        const double arg = (ik_vec * (R * latvec)) * TWO_PI;
                        const complex<double> kphase = complex<double>(cos(arg), sin(arg));
                        // cout<<"  I J:  "<<I<<"  "<<J<<endl;
                        // cout<<chi0_k[freq][ikp.first][I][J].nr<<"  "<<chi0_k[freq][ikp.first][I][J].nc<<"  "<<chi0_mat.nr<<"  "<<chi0_mat.nc<<endl;
                        for (size_t ci = 0; ci != mat_size; ci++)
                            chi0_k.at(freq).at(ikp.first).at(I).at(J).c[ci] += kphase * chi0_mat.c[ci];
                    }
                }
            }
        }
        temp_chi0_k_count++;
    }
}

void Cal_Periodic_Chi0::cal_MP2_energy_pi_k(map<double, double> &freq_grid)
{
    // std::map<double,map<Abfs::Vector3_Order<double>,map<size_t,map<size_t,ComplexMatrix>>>> pi_k(cal_pi_k());
    std::map<double, map<Vector3_Order<double>, map<size_t, map<size_t, ComplexMatrix>>>> pi_k_square;
    for (auto &freq_pair : pi_k)
    {
        const double freq = freq_pair.first;
        const auto pi_omega = freq_pair.second;
        for (auto &k_pair : pi_omega)
        {
            const auto k = k_pair.first;
            auto pi_omega_k = k_pair.second;
            for (auto &I_pair : pi_omega_k)
            {
                const size_t I = I_pair.first;
                auto pi_omega_k_I = I_pair.second;
                const size_t mu_num = atom_mu[I];
                pi_k_square[freq][k][I][I].create(mu_num, mu_num);
                for (auto &J_pair : pi_omega_k_I)
                {
                    const size_t J = J_pair.first;
                    const auto pi_mat = J_pair.second;
                    const auto pi2_mat = pi_omega_k[J][I];
                    const size_t nu_num = atom_mu[J];
                    // cout<<freq<<"    "<<R<<"   "<<I<<"    "<<J<<"    mu_num="<<mu_num<<"    nu_num="<<nu_num<<endl;

                    const size_t xi_num = atom_mu[J];

                    LapackConnector::gemm(
                        'N', 'N',
                        mu_num, mu_num, xi_num,
                        1,
                        pi_mat.c, xi_num,
                        pi2_mat.c, mu_num,
                        1,
                        pi_k_square[freq][k][I][I].c,
                        mu_num);
                }
                // print_complex_matrix("pi_square_mat",pi_square[freq][R][I][I].nc,pi_square[freq][R][I][I].nr,pi_square[freq][R][I][I].c,pi_square[freq][R][I][I].nc);
            }
        }
    }

    complex<double> tot_trace_pi(0.0, 0.0);
    for (auto &freq_pair : pi_k_square)
    {
        auto freq = freq_pair.first;
        auto pi_omega = freq_pair.second;
        complex<double> trace_pi(0.0, 0.0);
        for (auto &k_pair : pi_omega)
        {
            auto pi_omega_k = k_pair.second;
            for (auto &I_pair : pi_omega_k)
            {
                const size_t I = I_pair.first;
                auto pi_omega_k_I = I_pair.second;
                trace_pi += trace(pi_omega_k_I[I]);
            }
        }
        tot_trace_pi += trace_pi * freq_grid[freq] * (double(meanfield.get_n_spins()) / meanfield.get_n_kpoints() / TWO_PI / 2 * (-1));
        cout << " pi_k   freq:" << freq << "   trace:" << trace_pi << endl;
    }
    cout << " tot_MP2_energy_pi_k: " << tot_trace_pi;
}

void Cal_Periodic_Chi0::RPA_correlation_energy(map<double, double> &freq_grid)
{
    printf("Begin cal cRPA , pid:  %d\n", para_mpi.get_myid());
    int range_all = 0;

    for (auto &iat : atom_mu)
    {
        range_all += iat.second;
    }

    vector<int> part_range;
    part_range.resize(atom_mu.size());
    part_range[0] = 0;
    int count_range = 0;
    for (int I = 0; I != atom_mu.size() - 1; I++)
    {
        count_range += atom_mu[I];
        part_range[I + 1] = count_range;
    }

    cout << "part_range:" << endl;
    for (int I = 0; I != atom_mu.size(); I++)
    {
        cout << part_range[I] << endl;
    }
    cout << "part_range over" << endl;
    // std::map<double,map<Abfs::Vector3_Order<double>,map<size_t,map<size_t,ComplexMatrix>>>> pi_k(cal_pi_k());
    std::map<double, map<Vector3_Order<double>, ComplexMatrix>> pi_freq_2D;

    std::map<double, map<Vector3_Order<double>, ComplexMatrix>> chi0_2D;
    std::map<Vector3_Order<double>, ComplexMatrix> Vq_2D;
    for (const auto &freq_pair : pi_k)
    {
        const auto freq = freq_pair.first;
        const auto pi_freq = freq_pair.second;

        for (const auto &k_pair : pi_freq)
        {
            const auto kvec_c = k_pair.first;
            const auto pi_freq_k = k_pair.second;
            // chi0_2D[freq][kvec_c].create(range_all,range_all);
            pi_freq_2D[freq][kvec_c].create(range_all, range_all);

            ComplexMatrix tmp_pi_2D(range_all, range_all);
            tmp_pi_2D.zero_out();
            for (const auto &I_pair : pi_freq_k)
            {
                const auto I = I_pair.first;
                const auto pi_freq_k_I = I_pair.second;
                const size_t mu_num = atom_mu[I];
                for (const auto &J_pair : pi_freq_k_I)
                {
                    const auto J = J_pair.first;
                    const auto mat = J_pair.second;
                    const size_t nu_num = atom_mu[J];
                    // if (freq == first_freq)
                    // {
                    //     complex<double> trace_pi;
                    //     trace_pi = trace(mat);
                    //     // cout << "final pi   IJ:  " << I << "   " << J << "    kvec: " << kvec_c << "   trace: " << trace_pi << endl;
                    //     // print_complex_real_matrix(" final pi ", mat);
                    // }

                    for (size_t mu = 0; mu != mu_num; ++mu)
                    {
                        for (size_t nu = 0; nu != nu_num; ++nu)
                        {
                            // pi_freq_2D.at(freq).at(kvec_c)(part_range[I]+mu,part_range[J]+nu)+=mat(mu,nu);
                            tmp_pi_2D(part_range[I] + mu, part_range[J] + nu) += mat(mu, nu);
                        }
                    }
                    //     if(I<=J)
                    //     {
                    //         for(size_t mu=0;mu!=mu_num;++mu)
                    //         {
                    //             for(size_t nu=0;nu!=nu_num;++nu)
                    //             {
                    //                 chi0_2D.at(freq).at(kvec_c)(part_range[I]+mu,part_range[J]+nu)+=chi0_k.at(freq).at(kvec_c).at(I).at(J)(mu,nu);
                    //             }
                    //         }
                    //     }
                    //     else
                    //     {
                    //         for(size_t mu=0;mu!=mu_num;++mu)
                    //         {
                    //             for(size_t nu=0;nu!=nu_num;++nu)
                    //             {
                    //                 chi0_2D.at(freq).at(kvec_c)(part_range[I]+mu,part_range[J]+nu)+=conj(chi0_k.at(freq).at(kvec_c).at(J).at(I)(nu,mu));
                    //             }
                    //         }
                    //     }
                }
            }
            if (parall_type == "atom_pair")
            {
                para_mpi.reduce_ComplexMatrix(tmp_pi_2D, pi_freq_2D.at(freq).at(kvec_c));
            }
            else
            {
                pi_freq_2D.at(freq).at(kvec_c) = std::move(tmp_pi_2D);
            }
            // para_mpi.allreduce_ComplexMatrix(tmp_pi_2D,pi_freq_2D.at(freq).at(kvec_c));
            // cout<<"  freq:   "<<freq;
            // rt_cm_max(pi_freq_2D.at(freq).at(kvec_c));
            //  if(freq==first_freq)
            //  {
            //      cout<<"  pi_2d mat , ik_vec:  "<<kvec_c<<endl;
            //      print_complex_matrix("pi_2d",pi_freq_2D[freq][kvec_c]);
            //  }
        }
    }
    // for(auto &qvec:irk_weight)
    // {
    //     auto kvec_c=qvec.first;
    //     Vq_2D[kvec_c].create(range_all,range_all);
    //     for(auto &I_p:pi_k[first_freq][kvec_c])
    //     {
    //         const auto I=I_p.first;
    //         const size_t mu_num=atom_mu[I];
    //         for(auto &J_p:I_p.second)
    //         {
    //             const auto J=J_p.first;
    //             const size_t nu_num=atom_mu[J];
    //             if(I<=J)
    //             {
    //                 for(size_t mu=0;mu!=mu_num;++mu)
    //                 {
    //                     for(size_t nu=0;nu!=nu_num;++nu)
    //                     {
    //                         Vq_2D[kvec_c](part_range[I]+mu,part_range[J]+nu)+=(*Vq[I][J][kvec_c])(mu,nu);
    //                     }
    //                 }
    //             }
    //             else
    //             {
    //                 for(size_t mu=0;mu!=mu_num;++mu)
    //                 {
    //                     for(size_t nu=0;nu!=nu_num;++nu)
    //                     {
    //                         Vq_2D[kvec_c](part_range[I]+mu,part_range[J]+nu)+=conj((*Vq[J][I][kvec_c])(nu,mu));
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }

    // ofstream fs;
    // fs.open("librpa_pi.txt");
    // for(auto &irk:irk_weight)
    // {
    //     auto kvec_c=irk.first;
    //     fs<<kvec_c<<endl;
    //     print_complex_matrix_file("Vq_2D",Vq_2D[kvec_c],fs);
    //     print_complex_matrix_file("chi0_2D",chi0_2D[first_freq][kvec_c],fs);
    //     print_complex_matrix_file("pi_2D",pi_freq_2D[first_freq][kvec_c],fs);

    // }
    // fs.close();
    if (para_mpi.get_myid() == 0)
    {
        complex<double> tot_RPA_energy(0.0, 0.0);
        map<Vector3_Order<double>, complex<double>> cRPA_k;
        for (const auto &freq_pair : pi_freq_2D)
        {
            const auto freq = freq_pair.first;
            const auto freq_k = freq_pair.second;
            const double weight = freq_grid[freq];
            for (const auto &k_pair : freq_k)
            {
                const auto kvec_c = k_pair.first;
                const auto mat = k_pair.second;
                complex<double> rpa_for_omega_k(0.0, 0.0);
                ComplexMatrix identity(range_all, range_all);
                ComplexMatrix identity_minus_pi(range_all, range_all);

                identity.set_as_identity_matrix();

                identity_minus_pi = identity - pi_freq_2D[freq][kvec_c];
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
                trace_pi = trace(pi_freq_2D.at(freq).at(kvec_c));
                cout << "PI trace vector:" << endl;
                // if(freq==first_freq)
                // {
                //     for(int ir=0;ir!=pi_freq_2D.at(freq).at(kvec_c).nr;ir++)
                //         cout<<pi_freq_2D.at(freq).at(kvec_c)(ir,ir).real()<<"    ";
                // }
                cout << endl;
                rpa_for_omega_k = ln_det + trace_pi;
                cout << " freq:" << freq << "      rpa_for_omega_k: " << rpa_for_omega_k << "      lnt_det: " << ln_det << "    trace_pi " << trace_pi << endl;
                cRPA_k[kvec_c] += rpa_for_omega_k * weight * irk_weight[kvec_c] * double(meanfield.get_n_spins()) / TWO_PI;
                tot_RPA_energy += rpa_for_omega_k * weight * irk_weight[kvec_c] * double(meanfield.get_n_spins()) / TWO_PI;
            }
        }

        for (auto &ikvec : cRPA_k)
        {
            cout << ikvec.first << ikvec.second << endl;
        }
        cout << "gx_num_" << freq_grid.size() << "  tot_RPA_energy:  " << setprecision(8) << tot_RPA_energy << endl;
    }
}

vector<double> Cal_Periodic_Chi0::read_cosine_trans_grid(const string &file_path)
{
    ifstream infile;
    infile.open(file_path);

    vector<double> tran;
    string s;
    int n_freq = 0;
    // stringstream ss;
    // double ss_d;
    while (getline(infile, s))
    {
        if (s[0] == 'F')
        {
            n_freq++;
        }
        else
        {
            stringstream ss(s);
            double ss_d;
            ss >> ss_d;
            tran.push_back(ss_d);
        }
    }
    infile.close();

    double gap;
    double Emin, Emax;
    gap = meanfield.get_E_min_max(Emin, Emax);
    cout << "Cosine_tran_grid" << endl;
    for (int i = 0; i != tran.size(); i++)
    {
        tran[i] /= gap;
        // cout<<tran[i]<<endl;
    }
    return tran;
}

map<double, double> Cal_Periodic_Chi0::read_file_grid(const string &file_path, const char type)
{
    map<double, double> grid;
    ifstream infile;
    infile.open(file_path);

    vector<double> tran;
    string s0;
    int n_freq = 0;
    // stringstream ss;
    // double ss_d;
    string pattern = "          ";
    while (getline(infile, s0))
    {
        std::string::size_type pos;
        std::vector<std::string> result;
        s0 += pattern; //扩展字符串以方便操作

        for (int i = 0; i != s0.size(); i++)
        {
            pos = s0.find(pattern, i);
            if (pos < s0.size())
            {
                std::string s = s0.substr(i, pos - i);
                result.push_back(s);
                i = pos + pattern.size() - 1;
            }
        }
        grid.insert(pair<double, double>(std::stod(result[0]), std::stod(result[1])));
    }
    infile.close();

    double gap;
    double Emin, Emax;
    gap = meanfield.get_E_min_max(Emin, Emax);
    map<double, double> minimax_grid;
    minimax_grid.clear();
    switch (type)
    {
    case 'F':
    {
        for (auto &i_pair : grid)
            minimax_grid.insert({i_pair.first * gap, i_pair.second * gap * PI});

        cout << " MINIMAX_GRID_Freq " << endl;
        for (const auto &m : minimax_grid)
            cout << m.first << "      " << m.second << endl;
        break;
    }

    case 'T':
    {
        for (auto i_pair : grid)
            minimax_grid.insert({i_pair.first / (gap), i_pair.second / (gap)});

        cout << " MINIMAX_GRID_Tau " << endl;
        for (const auto &m : minimax_grid)
            cout << m.first << "      " << m.second << endl;
        break;
    }
    }

    return minimax_grid;
}

void Cal_Periodic_Chi0::print_matrix(const char *desc, const matrix &mat)
{
    int nr = mat.nr;
    int nc = mat.nc;
    printf("\n %s\n", desc);
    for (int i = 0; i < nr; i++)
    {
        for (int j = 0; j < nc; j++)
            printf(" %9.5f", mat.c[i * nc + j]);
        printf("\n");
    }
}
void Cal_Periodic_Chi0::print_complex_matrix(const char *desc, const ComplexMatrix &mat)
{
    int nr = mat.nr;
    int nc = mat.nc;
    printf("\n %s\n", desc);
    for (int i = 0; i < nr; i++)
    {
        for (int j = 0; j < nc; j++)
            printf("%10.6f,%8.6f ", mat.c[i * nc + j].real(), mat.c[i * nc + j].imag());
        printf("\n");
    }
}
void Cal_Periodic_Chi0::print_complex_real_matrix(const char *desc, const ComplexMatrix &mat)
{
    int nr = mat.nr;
    int nc = mat.nc;
    printf("\n %s\n", desc);
    for (int i = 0; i < nr; i++)
    {
        for (int j = 0; j < nc; j++)
            printf(" %10.6f", mat.c[i * nc + j].real());
        printf("\n");
    }
}

void Cal_Periodic_Chi0::print_complex_matrix_file(const char *desc, const ComplexMatrix &mat, ofstream &fs)
{
    int nr = mat.nr;
    int nc = mat.nc;
    fs << desc << endl;
    for (int i = 0; i < nr; i++)
    {
        for (int j = 0; j < nc; j++)
            fs << setw(15) << showpoint << fixed << setprecision(8) << mat.c[i * nc + j].real() << setw(14) << showpoint << fixed << setprecision(8) << mat.c[i * nc + j].imag();
        fs << "\n";
    }
}

void Cal_Periodic_Chi0::dgeev(matrix &mat, vector<complex<double>> &egValue, matrix &vr, int info)
{
    assert(mat.nr == mat.nc);
    vector<double> eg_r(mat.nr);
    vector<double> eg_i(mat.nr);
    matrix vl(mat.nr, mat.nc);
    vr.create(mat.nr, mat.nc);
    int lwork = 16 * mat.nr;
    vector<double> work(std::max(1, lwork));
    dgeev_("N", "V", &mat.nr, mat.c, &mat.nr, &eg_r[0], &eg_i[0], vl.c, &mat.nr, vr.c, &mat.nr, &work[0], &lwork, &info);
    for (int i = 0; i != mat.nr; i++)
        egValue.push_back(complex<double>(eg_r[i], eg_i[i]));
    vr = transpose(vr);
}

Cal_Periodic_Chi0 cal_chi0;
