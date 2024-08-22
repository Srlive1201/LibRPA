#include "librpa_main.h"

#include <algorithm>
#if defined(__MACH__)
#include <malloc/malloc.h>  // for malloc_zone_pressure_relief and malloc_default_zone
#else
#include <malloc.h>
#endif

#include "chi0.h"
#include "epsilon.h"
#include "meanfield.h"
#include "parallel_mpi.h"
#include "params.h"
#include "pbc.h"
#include "profiler.h"
#include "envs_io.h"
#include "envs_mpi.h"
#include "utils_io.h"

void librpa_main()
{
    using LIBRPA::envs::mpi_comm_global_h;
    using LIBRPA::envs::blacs_ctxt_global_h;
    using LIBRPA::utils::lib_printf;
    using LIBRPA::envs::ofs_myid;
    //printf("AFTER init MPI  myid: %d\n",mpi_comm_global_h.myid);
    // mpi_comm_global_h.barrier();

    if (mpi_comm_global_h.is_root())
    {
        LIBRPA::utils::lib_printf("LibRPA control parameters:\n");
        Params::print();
    }
    for(auto &p:atom_nw)
        ofs_myid<<"atom_nw I nw:"<<p.first<<" "<<p.second<<endl;
    for(auto &p:atom_mu)
        ofs_myid<<"atom_mu I mu:"<<p.first<<" "<<p.second<<endl;
    init_N_all_mu();
    ofs_myid<<"after init_N_all_mu"<<endl;
    for(auto &p:atom_nw)
        ofs_myid<<"atom_nw I nw:"<<p.first<<" "<<p.second<<endl;
    for(auto &p:atom_mu)
        ofs_myid<<"atom_mu I mu:"<<p.first<<" "<<p.second<<endl;
    LIBRPA::atomic_basis_wfc.set(atom_nw);
    LIBRPA::atomic_basis_abf.set(atom_mu);
   
    // for(auto &Ip:Cs)
    //     for(auto &Jp:Ip.second)
    //         for(auto &Rp:Jp.second)
    //         {
    //             auto R=Rp.first;
    //             LIBRPA::utils::lib_printf("Cs  myid : %d   I: %d, J: %d, R:(%d, %d, %d)\n",mpi_comm_global_h.myid, Ip.first,Jp.first,R.x,R.y,R.z);
    //         }

    
    if(Vq_block_loc.size()>0)
    {
        meanfield.allredue_wfc_isk();
        allreduce_atp_aux();
        allreduce_2D_coulomb_to_atompair(Vq_block_loc,Vq,Params::gf_R_threshold);
    }
    for (auto &Ip : Vq)
        for (auto &Jp : Ip.second)
            for (auto &qp : Jp.second)
            {
                auto q = qp.first;
                // LIBRPA::utils::lib_printf("allreduce Vq myid : %d   I: %d, J: %d, q:(%f, %f, %f)\n",mpi_comm_global_h.myid, Ip.first,Jp.first,q.x,q.y,q.z);
                // print_complex_matrix("vq",*qp.second);
            }
    // for(auto vq_p:Vq)
    // {
    //     auto I=vq_p.first;
    //     for(auto Ip:vq_p.second)
    //     {
    //         auto J=Ip.first;
    //         for(auto Jp:Ip.second)
    //         {
    //             auto kvec=Jp.first;
    //             cout<<"Vq  I J kvec:"<<I<<" "<<J<<" "<<kvec<<endl;
    //             print_complex_matrix("Vq",*Jp.second);
    //         }
    //     }
    // }
    //mpi_comm_global_h.barrier();
    blacs_ctxt_global_h.init();
    blacs_ctxt_global_h.set_square_grid();

    // Set output
    // The output now should be set using set_librpa_params
    // Params::output_file = "LibRPA_output.txt";

    Profiler::start("total", "Total");

    Profiler::start("driver_io_init", "Driver IO Initialization");
    //parse_inputfile_to_params(input_filename);
    // create output directory, only by the root process

    if (mpi_comm_global_h.is_root())
        system(("mkdir -p " + Params::output_dir).c_str());
    //mpi_comm_global_h.barrier();
    // Params::nfreq = stoi(argv[1]);
    // Params::gf_R_threshold = stod(argv[2]);
    Params::check_consistency();
    if (mpi_comm_global_h.is_root())
        Params::print();
    //mpi_comm_global_h.barrier();
    if (mpi_comm_global_h.is_root())
    {
        cout << "Information of mean-field starting-point" << endl;
        cout << "| number of spins: " << meanfield.get_n_spins() << endl
             << "| number of k-points: " << meanfield.get_n_kpoints() << endl
             << "| number of bands: " << meanfield.get_n_bands() << endl
             << "| number of NAOs: " << meanfield.get_n_aos() << endl
             << "| efermi (Ha): " << meanfield.get_efermi() << endl;

        double emin, emax;
        meanfield.get_E_min_max(emin, emax);
        cout << "| Minimal transition energy (Ha): " << emin << endl
             << "| Maximal transition energy (Ha): " << emax << endl << endl;
    }
    Vector3_Order<int> period {kv_nmp[0], kv_nmp[1], kv_nmp[2]};
    auto Rlist = construct_R_grid(period);

    if (mpi_comm_global_h.is_root())
    {
        cout << "Lattice vectors (Bohr)" << endl;
        latvec.print();
        cout << "Reciprocal lattice vectors (2PI Bohr^-1)" << endl;
        G.print();
        LIBRPA::utils::lib_printf("kgrids: %2d %2d %2d\n", kv_nmp[0], kv_nmp[1], kv_nmp[2]);
        cout << "k-points read (Cartisian in 2Pi Bohr^-1 | fractional):" << endl;
        for (int ik = 0; ik != meanfield.get_n_kpoints(); ik++)
        {
            LIBRPA::utils::lib_printf("ik %4d: %10.7f %10.7f %10.7f | %10.7f %10.7f %10.7f\n",
                   ik+1, klist[ik].x, klist[ik].y, klist[ik].z,
                   kfrac_list[ik].x, kfrac_list[ik].y, kfrac_list[ik].z);
        }
        cout << "R-points to compute:" << endl;
        for (int iR = 0; iR != Rlist.size(); iR++)
        {
            LIBRPA::utils::lib_printf("iR %4d: %3d %3d %3d\n", iR+1, Rlist[iR].x, Rlist[iR].y, Rlist[iR].z);
        }
        cout << endl;
    }

    const int Rt_num = Rlist.size() * Params::nfreq;
    tot_atpair = generate_atom_pair_from_nat(natom, false);
    tot_atpair_ordered = generate_atom_pair_from_nat(natom, true);
  
    if (mpi_comm_global_h.is_root())
    {
        cout << "| Number of atoms: " << natom << endl;
        cout << "| Total atom pairs (unordered): " << tot_atpair.size() << endl;
        cout << "| R-tau                       : " << Rt_num << endl;
        cout << "| Total atom pairs (ordered)  : " << tot_atpair_ordered.size() << endl;
    }
    set_parallel_routing(Params::parallel_routing, tot_atpair.size(), Rt_num, LIBRPA::parallel_routing);

    // barrier to wait for information print on master process
    //mpi_comm_global_h.barrier();

    //para_mpi.chi_parallel_type=Parallel_MPI::parallel_type::ATOM_PAIR;
    // vector<atpair_t> local_atpair;
    if(LIBRPA::parallel_routing == LIBRPA::ParallelRouting::ATOM_PAIR)
    {
        // vector<int> atoms_list(natom);
        // for(int iat=0;iat!=natom;iat++)
        //     atoms_list[iat]=iat;
        // std::array<int,3> period_arr{kv_nmp[0], kv_nmp[1], kv_nmp[2]};
        // std::pair<std::vector<int>, std::vector<std::vector<std::pair<int, std::array<int, 3>>>>> list_loc_atp
        // = RI::Distribute_Equally::distribute_atoms(mpi_comm_global_h.comm, atoms_list, period_arr, 2, false);
        // // for(auto &atp:list_loc_atp.first)
        // //     LIBRPA::utils::lib_printf("| myid: %d   atp.first: %d\n",mpi_comm_global_h.myid,atp);
        // // for(auto &atps:list_loc_atp.second)
        // //     LIBRPA::utils::lib_printf("| myid: %d   atp.second: %d\n",mpi_comm_global_h.myid,atps);
        // std::ofstream ofs("out."+std::to_string(RI::MPI_Wrapper::mpi_get_rank(mpi_comm_global)));
        // for(auto &af:list_loc_atp.first)
		//     ofs<<af<<"  ";
		// ofs<<endl;
        // for( auto &a1 : list_loc_atp.second)
        //     for(auto &a2:a1)
        //         ofs<<a2.first<<"  ("<<a2.second[0]<<", "<<a2.second[1]<<", "<<a2.second[2]<<" )"<<endl;
        auto trangular_loc_atpair= dispatch_upper_trangular_tasks(natom,blacs_ctxt_global_h.myid,blacs_ctxt_global_h.nprows,blacs_ctxt_global_h.npcols,blacs_ctxt_global_h.myprow,blacs_ctxt_global_h.mypcol);
        //local_atpair = dispatch_vector(tot_atpair, mpi_comm_global_h.myid, mpi_comm_global_h.nprocs, true);
        for(auto &iap:trangular_loc_atpair)
            local_atpair.push_back(iap);
        LIBRPA::utils::lib_printf("| process %d , local_atom_pair size:  %zu\n", mpi_comm_global_h.myid, local_atpair.size());

        LIBRPA::utils::lib_printf("| process %d, size of Cs from local_atpair: %lu\n", mpi_comm_global_h.myid, Cs_data.data_IJR.size());
        // for(auto &ap:local_atpair)
        //     LIBRPA::utils::lib_printf("   |process %d , local_atom_pair:  %d,  %d\n", mpi_comm_global_h.myid,ap.first,ap.second);
    }
    else
    {
        local_atpair = generate_atom_pair_from_nat(natom, false);
        if(Params::DFT_software == "ABACUS")
            allreduce_atp_coulomb(Vq);
    }
    // debug, check available Coulomb blocks on each process
    LIBRPA::envs::ofs_myid << "Read Coulomb blocks in process\n";
    for (const auto& IJqcoul: Vq)
    {
        const auto& I = IJqcoul.first;
        for (const auto& Jqcoul: IJqcoul.second)
        {
            const auto& J = Jqcoul.first;
            for (const auto& qcoul: Jqcoul.second)
            {
                LIBRPA::envs::ofs_myid << I << " " << J << " " << qcoul.first << "\n";
            }
        }
    }
    std::flush(LIBRPA::envs::ofs_myid);
    mpi_comm_global_h.barrier();
    Profiler::stop("driver_io_init");

    // malloc_trim(0);
    // para_mpi.mpi_barrier();
    // if(para_mpi.is_master())
    //     system("free -m");
    //erase_Cs_from_local_atp(Cs,local_atpair);
    // malloc_trim(0);
    // para_mpi.mpi_barrier();
    // if(para_mpi.is_master())
    //     system("free -m");
    //READ_Vq_Row("./", "coulomb_mat", Params::vq_threshold, Vq, local_atpair);
    /* if(argv[1][0]=='0') */
    /*     ap_chi0.chi0_main(argv[1],argv[2]);  */
    /* else */
        /* cal_chi0.chi0_main(argv[1],argv[2]);  */
    /* return 0; */
    // try the new version

    Chi0 chi0(meanfield, klist, Params::nfreq);
    chi0.gf_R_threshold = Params::gf_R_threshold;

    // build ABF IJ and qlist from Vq

    vector<Vector3_Order<double>> qlist;

    for ( auto q_weight: irk_weight)
    {
        qlist.push_back(q_weight.first);
    }

    mpi_comm_global_h.barrier(); // FIXME: barrier seems not work here...
    // cout << endl;
    if (mpi_comm_global_h.myid == 0)
        cout << "Initialization finished, start task job from myid\n";

    if ( Params::task != "exx" )
    {
        Profiler::start("chi0_build", "Build response function chi0");
        chi0.build(Cs_data, Rlist, period, local_atpair, qlist,
                   TFGrids::get_grid_type(Params::tfgrids_type), true);
        Profiler::stop("chi0_build");
    }
    if (Params::debug)
    { // debug, check chi0
        char fn[80];
        for (const auto &chi0q: chi0.get_chi0_q())
        {
            const int ifreq = chi0.tfg.get_freq_index(chi0q.first);
            for (const auto &q_IJchi0: chi0q.second)
            {
                const int iq = std::distance(klist.begin(), std::find(klist.begin(), klist.end(), q_IJchi0.first));
                for (const auto &I_Jchi0: q_IJchi0.second)
                {
                    const auto &I = I_Jchi0.first;
                    for (const auto &J_chi0: I_Jchi0.second)
                    {
                        const auto &J = J_chi0.first;
                        sprintf(fn, "chi0fq_ifreq_%d_iq_%d_I_%zu_J_%zu_id_%d.mtx", ifreq, iq, I, J, mpi_comm_global_h.myid);
                        print_complex_matrix_mm(J_chi0.second, Params::output_dir + "/" + fn, 1e-15);
                    }
                }
            }
        }
    }
    //malloc_trim(0);

    // para_mpi.mpi_barrier();
    // if(para_mpi.is_master())
    //     system("free -m");
    // FIXME: a more general strategy to deal with Cs
    // Cs is not required after chi0 is computed in rpa task
    if (Params::task == "rpa")
    {
        // for(auto &Cp:Cs)
        // {
        //     Cs.erase(Cp.first);
        // }
        Cs_data.clear();
        LIBRPA::utils::lib_printf("Cs have been cleaned!\n");
    }

    #ifndef __MACH__
    malloc_trim(0);
    #else
    malloc_zone_pressure_relief(malloc_default_zone(), 0);
    #endif
    // RPA total energy
    if ( Params::task == "rpa" )
    {
        mpi_comm_global_h.barrier();
        Profiler::start("EcRPA", "Compute RPA correlation Energy");
        CorrEnergy corr;
        if (Params::use_scalapack_ecrpa && LIBRPA::parallel_routing == LIBRPA::ParallelRouting::ATOM_PAIR)
        {
            if(meanfield.get_n_kpoints() == 1)
                corr = compute_RPA_correlation_blacs_2d_gamma_only(chi0, Vq);
            else
                corr = compute_RPA_correlation_blacs_2d(chi0, Vq);
        }
        else
            corr = compute_RPA_correlation(chi0, Vq);

        if (mpi_comm_global_h.is_root())
        {
            lib_printf("RPA correlation energy (Hartree)\n");
            lib_printf("| Weighted contribution from each k:\n");
            for (const auto& q_ecrpa: corr.qcontrib)
            {
                cout << "| " << q_ecrpa.first << ": " << q_ecrpa.second << endl;
            }
            lib_printf("| Total EcRPA: %18.9f\n", corr.value.real());
            if (std::abs(corr.value.imag()) > 1.e-3)
                lib_printf("Warning: considerable imaginary part of EcRPA = %f\n", corr.value.imag());
        }
        Profiler::stop("EcRPA");
    }
    else
    {

    }

    Profiler::stop("total");

    if (mpi_comm_global_h.is_root())
    {
        Profiler::display();
        lib_printf("libRPA finished\n");
    }

    return;
}
