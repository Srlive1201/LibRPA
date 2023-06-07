#include "lib_main.h"
#include "constants.h"
#include "interpolate.h"
#include "fitting.h"
#include "dielecmodel.h"
#include "scalapack_connector.h"

#include <algorithm>
#include <malloc.h>
#include <set>

int main(int argc, char **argv)
{
    using namespace LIBRPA;
    MPI_Wrapper::init(argc, argv);

    mpi_comm_world_h.init();
    cout << mpi_comm_world_h.str() << endl;
    mpi_comm_world_h.barrier();
    blacs_ctxt_world_h.init();
    blacs_ctxt_world_h.set_square_grid();

    Profiler::start("total", "Total");

    Profiler::start("driver_io_init", "Driver IO Initialization");
    parse_inputfile_to_params(input_filename);
    // create output directory, only by the root process
    if (mpi_comm_world_h.is_root())
        system(("mkdir -p " + Params::output_dir).c_str());
    mpi_comm_world_h.barrier();
    Params::nfreq = stoi(argv[1]);
    Params::gf_R_threshold = stod(argv[2]);
    Params::check_consistency();
    if (mpi_comm_world_h.is_root())
        Params::print();
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

    const int Rt_num = Rlist.size() * Params::nfreq;

    READ_AIMS_EIGENVECTOR("./", meanfield);

    natom = get_natom_ncell_from_first_Cs_file();
    tot_atpair = generate_atom_pair_from_nat(natom, false);
    tot_atpair_ordered = generate_atom_pair_from_nat(natom, true);

    if (mpi_comm_world_h.is_root())
    {
        cout << "| Number of atoms: " << natom << endl;
        cout << "| Total atom pairs (unordered): " << tot_atpair.size() << endl;
        cout << "| R-tau                       : " << Rt_num << endl;
        cout << "| Total atom pairs (ordered)  : " << tot_atpair_ordered.size() << endl;
    }
    set_chi_parallel_type(Params::chi_parallel_routing, tot_atpair.size(), Rt_num, Params::use_libri_chi0);
    set_exx_parallel_type(Params::exx_parallel_routing, tot_atpair.size(), Rt_num, Params::use_libri_exx);
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
        
        READ_AIMS_Cs("./", Params::cs_threshold,local_atpair );
        printf("| process %d, size of Cs from local_atpair: %lu\n", LIBRPA::mpi_comm_world_h.myid, Cs.size());
        // for(auto &ap:local_atpair)
        //     printf("   |process %d , local_atom_pair:  %d,  %d\n", mpi_comm_world_h.myid,ap.first,ap.second);
        READ_Vq_Row("./", "coulomb_mat", Params::vq_threshold, Vq, local_atpair);
    }
    else
    {
        local_atpair = generate_atom_pair_from_nat(natom, false);
        READ_AIMS_Cs("./", Params::cs_threshold,local_atpair );
        printf("| process %d, size of Cs from local_atpair: %lu\n", LIBRPA::mpi_comm_world_h.myid, Cs.size());
        READ_Vq_Full("./", "coulomb_mat", Params::vq_threshold, Vq); 
    }
    mpi_comm_world_h.barrier();
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

    mpi_comm_world_h.barrier(); // FIXME: barrier seems not work here...
    // cout << endl;
    if (mpi_comm_world_h.myid == 0)
        cout << "Initialization finished, start task job from myid\n";

    if ( Params::task != "exx" )
    {
        Profiler::start("chi0_build", "Build response function chi0");
        chi0.build(Cs, Rlist, period, local_atpair, qlist,
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
                        sprintf(fn, "chi0fq_ifreq_%d_iq_%d_I_%zu_J_%zu_id_%d.mtx", ifreq, iq, I, J, mpi_comm_world_h.myid);
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
    if (LIBRPA::chi_parallel_type != LIBRPA::parallel_type::LIBRI_USED && Params::task == "rpa")
    {
        // for(auto &Cp:Cs)
        // {
        //     Cs.erase(Cp.first);
        // }
        Cs.clear();
    }

    malloc_trim(0);
    // RPA total energy
    if ( Params::task == "rpa" )
    {
        mpi_comm_world_h.barrier();
        Profiler::start("EcRPA", "Compute RPA correlation Energy");
        CorrEnergy corr;
        if (Params::use_scalapack_ecrpa && LIBRPA::chi_parallel_type == LIBRPA::parallel_type::ATOM_PAIR)
            corr = compute_RPA_correlation_blacs(chi0, Vq);
        else
            corr = compute_RPA_correlation(chi0, Vq);

        Profiler::stop("EcRPA");
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
        Profiler::stop("EcRPA");
    }
    else if ( Params::task == "g0w0" )
    {
        Profiler::start("g0w0", "G0W0 quasi-particle calculation");

        Profiler::start("read_vq_cut", "Load truncated Coulomb");
        READ_Vq_Full("./", "coulomb_cut_", Params::vq_threshold, Vq_cut); 
        const auto VR = FT_Vq(Vq_cut, Rlist, true);
        Profiler::stop("read_vq_cut");

        std::vector<double> epsmac_LF_imagfreq_re;

        if (Params::replace_w_head)
        {
            std::vector<double> omegas_dielect;
            std::vector<double> dielect_func;
            read_dielec_func("dielecfunc_out", omegas_dielect, dielect_func);

            switch(Params::option_dielect_func)
            {
                case 0:
                    {
                        assert(omegas_dielect.size() == chi0.tfg.size());
                        // TODO: check if frequencies and close
                        epsmac_LF_imagfreq_re = dielect_func;
                        break;
                    }
                case 1:
                    {
                        epsmac_LF_imagfreq_re = UTILS::interp_cubic_spline(omegas_dielect, dielect_func,
                                                                           chi0.tfg.get_freq_nodes());
                        break;
                    }
                case 2:
                    {
                        UTILS::LevMarqFitting levmarq;
                        // use double-dispersion Havriliak-Negami model
                        // initialize the parameters as 1.0
                        std::vector<double> pars(DoubleHavriliakNegami::d_npar, 1);
                        pars[0] = pars[4] = dielect_func[0];
                        epsmac_LF_imagfreq_re = levmarq.fit_eval(pars, omegas_dielect, dielect_func,
                                                                 DoubleHavriliakNegami::func_imfreq,
                                                                 DoubleHavriliakNegami::grad_imfreq,
                                                                 chi0.tfg.get_freq_nodes());
                        break;
                    }
                default:
                    throw std::logic_error("Unsupported value for option_dielect_func");
            }

            if (Params::debug)
            {
                if (mpi_comm_world_h.is_root())
                {
                    printf("Dielection function parsed:\n");
                    for (int i = 0; i < chi0.tfg.get_freq_nodes().size(); i++)
                        printf("%d %f %f\n", i+1, chi0.tfg.get_freq_nodes()[i], epsmac_LF_imagfreq_re[i]);
                }
                mpi_comm_world_h.barrier();
            }
        }

        Profiler::start("g0w0_exx", "Build exchange self-energy");
        auto exx = LIBRPA::Exx(meanfield, kfrac_list);
        exx.build_exx_orbital_energy(Cs, Rlist, period, VR);
        Profiler::stop("g0w0_exx");

        Profiler::start("g0w0_wc", "Build screened interaction");
        vector<std::complex<double>> epsmac_LF_imagfreq(epsmac_LF_imagfreq_re.cbegin(), epsmac_LF_imagfreq_re.cend());
        map<double, atom_mapping<std::map<Vector3_Order<double>, matrix_m<complex<double>>>>::pair_t_old> Wc_freq_q;
        if (Params::use_scalapack_gw_wc)
            Wc_freq_q = compute_Wc_freq_q_blacs(chi0, Vq, Vq_cut, epsmac_LF_imagfreq);
        else
            Wc_freq_q = compute_Wc_freq_q(chi0, Vq, Vq_cut, epsmac_LF_imagfreq);
        Profiler::stop("g0w0_wc");

        if (Params::debug)
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
                            print_matrix_mm_file(q_Wc.second, Params::output_dir + "/" + fn, 1e-15);
                        }
                    }
                }
            }
        }
        LIBRPA::G0W0 s_g0w0(meanfield, kfrac_list, chi0.tfg);
        Profiler::start("g0w0_sigc_IJ", "Build correlation self-energy");
        s_g0w0.build_spacetime_LibRI(Cs, Wc_freq_q, Rlist, period);
        Profiler::stop("g0w0_sigc_IJ");
        Profiler::start("g0w0_sigc_rotate_KS", "Rotate self-energy, IJ -> ij -> KS");
        s_g0w0.build_sigc_matrix_KS();
        Profiler::stop("g0w0_sigc_rotate_KS");

        Profiler::start("g0w0_export_sigc_KS", "Export self-energy in KS basis");
        mpi_comm_world_h.barrier();
        if (mpi_comm_world_h.is_root())
        {
            char fn[100];
            for (const auto &ispin_sigc: s_g0w0.sigc_is_ik_f_KS)
            {
                const auto &ispin = ispin_sigc.first;
                for (const auto &ik_sigc: ispin_sigc.second)
                {
                    const auto &ik = ik_sigc.first;
                    for (const auto &freq_sigc: ik_sigc.second)
                    {
                        const auto ifreq = s_g0w0.tfg.get_freq_index(freq_sigc.first);
                        sprintf(fn, "Sigc_fk_mn_ispin_%d_ik_%d_ifreq_%d.mtx", ispin, ik, ifreq);
                        print_matrix_mm_file(freq_sigc.second, Params::output_dir + "/" + fn, 1e-10);
                    }
                }
            }
            // for aims analytic continuation reader
            write_self_energy_omega("self_energy_omega.dat", s_g0w0);
        }
        Profiler::stop("g0w0_export_sigc_KS");
        Profiler::stop("g0w0");
    }
    else if ( Params::task == "exx" )
    {
        READ_Vq_Full("./", "coulomb_cut_", Params::vq_threshold, Vq_cut); 
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

    Profiler::stop("total");
    if(mpi_comm_world_h.is_root())
    {
        Profiler::display();
        printf("libRPA finished");
    }

    MPI_Wrapper::finalize();
    return 0;
}
