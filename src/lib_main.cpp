#include "lib_main.h"
#include <algorithm>
#include <malloc.h>
#include "scalapack_connector.h"
#include <set>
#include <RI/distribute/Distribute_Equally.h>

int main(int argc, char **argv)
{
    using namespace LIBRPA;
    MPI_Wrapper::init(argc, argv);

    mpi_comm_world_h.init();
    cout << mpi_comm_world_h.str() << endl;
    mpi_comm_world_h.barrier();
    blacs_ctxt_world_h.init();
    blacs_ctxt_world_h.set_square_grid();

    prof.add(0, "total", "Total");
    prof.add(1, "chi0_main", "Chi0 object");
    prof.add(2, "cal_Green_func_R_tau", "space-time Green's function");
    prof.add(2, "cal_Green_func",       "space-time Green's function");
    prof.add(2, "R_tau_routing", "Loop over R-tau");
    prof.add(2, "atom_pair_routing", "Loop over atom pairs");
    prof.add(2, "LibRI_rouing", "Loop over LibRI");
    prof.add(3, "cal_chi0_element", "chi(tau,R,I,J)");
    prof.add(4, "X");
    prof.add(4, "O");
    prof.add(4, "N");
    prof.add(4, "Z");
    prof.add(4, "reshape_Cs", "reshape Cs");
    prof.add(4, "reshape_mat", "reshape mat");
    prof.add(2,"EcRPA");
    prof.start("total");

    int flag;
    parse_inputfile_to_params(input_filename, params);
    // NOTE: to comply with command line input, may be removed later
    InputFile inputf;
    auto parser = inputf.load(input_filename, false);
    parser.parse_int("nfreq", params.nfreq, stoi(argv[1]), flag);
    parser.parse_double("gf_R_threshold", params.gf_R_threshold, stod(argv[2]), flag);
    params.check_consistency();
    if (mpi_comm_world_h.is_root())
        params.print();
    mpi_comm_world_h.barrier();

    READ_AIMS_BAND("band_out", meanfield);
    if (mpi_comm_world_h.is_root())
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

    READ_AIMS_STRU(meanfield.get_n_kpoints(), "stru_out");
    Vector3_Order<int> period {kv_nmp[0], kv_nmp[1], kv_nmp[2]};
    auto Rlist = construct_R_grid(period);

    if (mpi_comm_world_h.is_root())
    {
        cout << "Lattice vectors (Bohr)" << endl;
        latvec.print();
        cout << "Reciprocal lattice vectors (2PI Bohr^-1)" << endl;
        G.print();
        printf("kgrids: %2d %2d %2d\n", kv_nmp[0], kv_nmp[1], kv_nmp[2]);
        cout << "k-points read (Cartisian in 2Pi Bohr^-1 | fractional):" << endl;
        for (int ik = 0; ik != meanfield.get_n_kpoints(); ik++)
        {
            printf("ik %4d: %10.7f %10.7f %10.7f | %10.7f %10.7f %10.7f\n",
                   ik+1, kvec_c[ik].x, kvec_c[ik].y, kvec_c[ik].z,
                   kfrac_list[ik].x, kfrac_list[ik].y, kfrac_list[ik].z);
        }
        cout << "R-points to compute:" << endl;
        for (int iR = 0; iR != Rlist.size(); iR++)
        {
            printf("iR %4d: %3d %3d %3d\n", iR+1, Rlist[iR].x, Rlist[iR].y, Rlist[iR].z);
        }
        cout << endl;
    }

    const int Rt_num = Rlist.size() * params.nfreq;

    READ_AIMS_EIGENVECTOR("./", meanfield);

    
    natom=get_natom_ncell_from_first_Cs_file();
    tot_atpair = generate_atom_pair_from_nat(natom, false);
    tot_atpair_ordered = generate_atom_pair_from_nat(natom, true);

    if (mpi_comm_world_h.is_root())
    {
        cout << "| Number of atoms: " << natom << endl;
        cout << "| Total atom pairs (unordered): " << tot_atpair.size() << endl;
        cout << "| R-tau                       : " << Rt_num << endl;
        cout << "| Total atom pairs (ordered)  : " << tot_atpair_ordered.size() << endl;
    }
    set_chi_parallel_type(params.chi_parallel_routing, tot_atpair.size(), Rt_num, params.use_libri_chi0);
    set_exx_parallel_type(params.exx_parallel_routing, tot_atpair.size(), Rt_num, params.use_libri_exx);
    check_parallel_type();

    // barrier to wait for information print on master process
    mpi_comm_world_h.barrier();

    //para_mpi.chi_parallel_type=Parallel_MPI::parallel_type::ATOM_PAIR;
    // vector<atpair_t> local_atpair;
    if(chi_parallel_type == parallel_type::ATOM_PAIR)
    {
        // vector<int> atoms_list(natom);
        // for(int iat=0;iat!=natom;iat++)
        //     atoms_list[iat]=iat;
        // std::array<int,3> period_arr{kv_nmp[0], kv_nmp[1], kv_nmp[2]};
        // std::pair<std::vector<int>, std::vector<std::vector<std::pair<int, std::array<int, 3>>>>> list_loc_atp
        // = RI::Distribute_Equally::distribute_atoms(mpi_comm_world_h.comm, atoms_list, period_arr, 2, false);
        // // for(auto &atp:list_loc_atp.first)
        // //     printf("| myid: %d   atp.first: %d\n",mpi_comm_world_h.myid,atp);
        // // for(auto &atps:list_loc_atp.second)
        // //     printf("| myid: %d   atp.second: %d\n",mpi_comm_world_h.myid,atps);
        // std::ofstream ofs("out."+std::to_string(RI::MPI_Wrapper::mpi_get_rank(MPI_COMM_WORLD)));
        // for(auto &af:list_loc_atp.first)
		//     ofs<<af<<"  ";
		// ofs<<endl;
        // for( auto &a1 : list_loc_atp.second)
        //     for(auto &a2:a1)
        //         ofs<<a2.first<<"  ("<<a2.second[0]<<", "<<a2.second[1]<<", "<<a2.second[2]<<" )"<<endl;
        auto trangular_loc_atpair= dispatch_upper_trangular_tasks(natom,blacs_ctxt_world_h.myid,blacs_ctxt_world_h.nprows,blacs_ctxt_world_h.npcols,blacs_ctxt_world_h.myprow,blacs_ctxt_world_h.mypcol);
        //local_atpair = dispatch_vector(tot_atpair, mpi_comm_world_h.myid, mpi_comm_world_h.nprocs, true);
        for(auto &iap:trangular_loc_atpair)
            local_atpair.push_back(iap);
        printf("| process %d , local_atom_pair size:  %zu\n", mpi_comm_world_h.myid, local_atpair.size());
        
        READ_AIMS_Cs("./", params.cs_threshold,local_atpair );
        printf("| process %d, size of Cs from local_atpair: %lu\n", LIBRPA::mpi_comm_world_h.myid, Cs.size());
        // for(auto &ap:local_atpair)
        //     printf("   |process %d , local_atom_pair:  %d,  %d\n", mpi_comm_world_h.myid,ap.first,ap.second);
        READ_Vq_Row("./", "coulomb_mat", params.vq_threshold, Vq, local_atpair);
    }
    else
    {
        local_atpair = dispatch_vector(tot_atpair, mpi_comm_world_h.myid, mpi_comm_world_h.nprocs, true);
        READ_AIMS_Cs("./", params.cs_threshold,local_atpair );
        printf("| process %d, size of Cs from local_atpair: %lu\n", LIBRPA::mpi_comm_world_h.myid, Cs.size());
        READ_Vq_Full("./", "coulomb_mat", params.vq_threshold, Vq); 
    }
    mpi_comm_world_h.barrier();
    // malloc_trim(0);
    // para_mpi.mpi_barrier();
    // if(para_mpi.is_master())
    //     system("free -m");
    //erase_Cs_from_local_atp(Cs,local_atpair);
    // malloc_trim(0);
    // para_mpi.mpi_barrier();
    // if(para_mpi.is_master())
    //     system("free -m");
    //READ_Vq_Row("./", "coulomb_mat", params.vq_threshold, Vq, local_atpair);
    /* if(argv[1][0]=='0') */
    /*     ap_chi0.chi0_main(argv[1],argv[2]);  */
    /* else */
        /* cal_chi0.chi0_main(argv[1],argv[2]);  */
    /* return 0; */
    // try the new version
    Chi0 chi0(meanfield, klist, params.nfreq);
    chi0.gf_R_threshold = params.gf_R_threshold;

    // build ABF IJ and qlist from Vq
    
    vector<Vector3_Order<double>> qlist;

    for ( auto q_weight: irk_weight)
    {
        qlist.push_back(q_weight.first);
    }

    mpi_comm_world_h.barrier(); // FIXME: barrier seems not work here...
    // cout << endl;
    if (mpi_comm_world_h.myid == 0)
        cout << "Initialization finished, start task job from myid\n";

    if ( params.task != "exx" )
    {
        chi0.build(Cs, Rlist, period, local_atpair, qlist,
                   TFGrids::get_grid_type(params.tfgrids_type), true);
    }

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
                        sprintf(fn, "chi0fq_ifreq_%d_iq_%d_I_%zu_J_%zu_id_%d.mtx", ifreq, iq, I, J, mpi_comm_world_h.myid);
                        // print_complex_matrix_mm(J_chi0.second, fn, 1e-15);
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
    if (LIBRPA::chi_parallel_type != LIBRPA::parallel_type::LIBRI_USED && params.task == "rpa")
    {
        // for(auto &Cp:Cs)
        // {
        //     Cs.erase(Cp.first);
        // }
        Cs.clear();
    }

    malloc_trim(0);
    // RPA total energy
    if ( params.task == "rpa" )
    {
        prof.start("EcRPA");
        CorrEnergy corr;
        if (params.use_scalapack_ecrpa && LIBRPA::chi_parallel_type == LIBRPA::parallel_type::ATOM_PAIR)
            corr = compute_RPA_correlation_blacs(chi0, Vq);
        else
            corr = compute_RPA_correlation(chi0, Vq);

        prof.stop("EcRPA");
        if (mpi_comm_world_h.is_root())
        {
            printf("RPA correlation energy (Hartree)\n");
            printf("| Weighted contribution from each k:\n");
            for (const auto& q_ecrpa: corr.qcontrib)
            {
                cout << "| " << q_ecrpa.first << ": " << q_ecrpa.second << endl;
            }
            printf("| Total EcRPA: %18.9f\n", corr.value.real());
            if (std::abs(corr.value.imag()) > 1.e-3)
                printf("Warning: considerable imaginary part of EcRPA = %f\n", corr.value.imag());
        }
    }
    else if ( params.task == "g0w0" )
    {
        READ_Vq_Full("./", "coulomb_cut_", params.vq_threshold, Vq_cut); 
        const auto VR = FT_Vq(Vq_cut, Rlist, true);
        auto exx = LIBRPA::Exx(meanfield, klist);
        exx.build_exx_orbital_energy(Cs, Rlist, period, VR);
        vector<std::complex<double>> epsmac_LF_imagfreq;
        map<double, atpair_k_cplx_mat_t> Wc_freq_q;
        if (params.use_scalapack_gw_wc)
            Wc_freq_q = compute_Wc_freq_q_blacs(chi0, Vq, Vq_cut, epsmac_LF_imagfreq);
        else
            Wc_freq_q = compute_Wc_freq_q(chi0, Vq, Vq_cut, epsmac_LF_imagfreq);
        { // debug, check Wc
            char fn[80];
            for (const auto &Wc: Wc_freq_q)
            {
                const int ifreq = chi0.tfg.get_freq_index(Wc.first);
                for (const auto &I_JqWc: Wc.second)
                {
                    const auto &I = I_JqWc.first;
                    for (const auto &J_qWc: I_JqWc.second)
                    {
                        const auto &J = J_qWc.first;
                        for (const auto &q_Wc: J_qWc.second)
                        {
                            const int iq = std::distance(klist.begin(), std::find(klist.begin(), klist.end(), q_Wc.first));
                            sprintf(fn, "Wcfq_ifreq_%d_iq_%d_I_%zu_J_%zu_id_%d.mtx", ifreq, iq, I, J, mpi_comm_world_h.myid);
                            print_complex_matrix_mm(*q_Wc.second, fn, 1e-15);
                        }
                    }
                }
            }
        }
        const auto Wc_tau_R = CT_FT_Wc_freq_q(Wc_freq_q, chi0.tfg, Rlist);
    }
    else if ( params.task == "exx" )
    {
        READ_Vq_Full("./", "coulomb_cut_", params.vq_threshold, Vq_cut); 
        const auto VR = FT_Vq(Vq_cut, Rlist, true);
        auto exx = LIBRPA::Exx(meanfield, kfrac_list);
        exx.build_exx_orbital_energy(Cs, Rlist, period, VR);
        // FIXME: need to reduce first when MPI is used
        // NOTE: may extract to a common function
        if (LIBRPA::mpi_comm_world_h.is_root())
        {
            for (int isp = 0; isp != meanfield.get_n_spins(); isp++)
            {
                printf("Spin channel %1d\n", isp+1);
                for (int ik = 0; ik != meanfield.get_n_kpoints(); ik++)
                {
                    cout << "k-point " << ik + 1 << ": " << kvec_c[ik] << endl;
                    printf("%-4s  %-10s  %-10s\n", "Band", "e_exx (Ha)", "e_exx (eV)");
                    for (int ib = 0; ib != meanfield.get_n_bands(); ib++)
                        printf("%4d  %10.5f  %10.5f\n", ib+1, exx.Eexx[isp][ik][ib], HA2EV * exx.Eexx[isp][ik][ib]);
                    printf("\n");
                }
            }
        }
    }

    prof.stop("total");
    if(mpi_comm_world_h.is_root())
        prof.display();

    MPI_Wrapper::finalize();
    return 0;
}
