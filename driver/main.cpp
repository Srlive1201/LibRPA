#include <algorithm>
#include <limits>

#include "analycont.h"
#include "chi0.h"
#include "constants.h"
#include "coulmat.h"
#include "dielecmodel.h"
#include "envs_io.h"
#include "envs_mpi.h"
#include "epsilon.h"
#include "exx.h"
#include "fitting.h"
#include "librpa.h"
#include "gw.h"
#include "inputfile.h"
#include "interpolate.h"
#include "meanfield.h"
#include "params.h"
#include "pbc.h"
#include "profiler.h"
#include "qpe_solver.h"
#include "read_data.h"
#include "stl_io_helper.h"
#include "task.h"
#include "utils_io.h"
#include "utils_mem.h"
#include "write_aims.h"

#include "task_rpa.h"

static void initialize(int argc, char **argv)
{
    using namespace LIBRPA::envs;
    using LIBRPA::utils::lib_printf;
    // MPI Initialization
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if (MPI_THREAD_MULTIPLE != provided)
    {
        lib_printf("Warning: MPI_Init_thread provide %d != required %d", provided, MPI_THREAD_MULTIPLE);
    }

    initialize_librpa_environment(MPI_COMM_WORLD, 0, 0, "");

    // Global profiler begins right after MPI is initialized
    Profiler::start("total", "Total");
}

static void finalize()
{
    using namespace LIBRPA::envs;
    using LIBRPA::utils::lib_printf;

    Profiler::stop("total");
    if (mpi_comm_global_h.is_root())
    {
        Profiler::display();
        lib_printf("libRPA finished\n");
    }

    finalize_librpa_environment();

    MPI_Finalize();
}

int main(int argc, char **argv)
{
    using LIBRPA::envs::mpi_comm_global_h;
    using LIBRPA::envs::blacs_ctxt_global_h;
    using LIBRPA::ParallelRouting;
    using LIBRPA::parallel_routing;
    using LIBRPA::task_t;
    using LIBRPA::envs::ofs_myid;
    using LIBRPA::utils::lib_printf;

    initialize(argc, argv);
    cout << mpi_comm_global_h.str() << endl;
    mpi_comm_global_h.barrier();

    /*
     * HACK: A close-to-square process grid is imposed.
     *       This might subject to performance issue when
     *       the number of processes is not dividable.
     */
    blacs_ctxt_global_h.set_square_grid();

    /*
     * Load computational parameters from input file
     */
    Profiler::start("driver_read_params", "Driver Read Input Parameters");
    parse_inputfile_to_params(input_filename);
    // Backward compatibility for parsing number of frequencies and threshold from CLI
    if (argc > 2)
    {
        Params::nfreq = stoi(argv[1]);
        Params::gf_R_threshold = stod(argv[2]);
    }
    Params::check_consistency();

    // check and set computing task
    task_t task;
    std::string task_lower = Params::task;
    transform(task_lower.begin(), task_lower.end(), task_lower.begin(), ::tolower);
    if (task_lower == "rpa")
        task = task_t::RPA;
    else if (task_lower == "g0w0")
        task = task_t::G0W0;
    else if (task_lower == "exx")
        task = task_t::EXX;
    else if (task_lower == "wc_rf")
        task = task_t::Wc_Rf;
    else if (task_lower == "print_minimax")
        task = task_t::print_minimax;
    else
        throw std::logic_error("Unknown task (" + Params::task + "). Please check your input");


    if (mpi_comm_global_h.is_root())
    {
        system(("mkdir -p " + Params::output_dir).c_str());
        Params::print();
    }
    mpi_comm_global_h.barrier();
    Profiler::stop("driver_read_params");

    Profiler::start("driver_band_out", "Driver Read Meanfield band");
    read_band("band_out", meanfield);
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
    Profiler::stop("driver_band_out");

    if (task == task_t::print_minimax)
    {
        double emin, emax;
        meanfield.get_E_min_max(emin, emax);
        TFGrids tfg(Params::nfreq);
        tfg.generate_minimax(emin, emax);
        if (mpi_comm_global_h.is_root())
            tfg.show();
        finalize();
        return 0;
    }

    Profiler::start("driver_read_common_input_data", "Driver Read Task-Common Input Data");
    read_stru(meanfield.get_n_kpoints(), "stru_out");
    Vector3_Order<int> period {kv_nmp[0], kv_nmp[1], kv_nmp[2]};
    auto Rlist = construct_R_grid(period);
    const int Rt_num = Rlist.size() * Params::nfreq;
    if (mpi_comm_global_h.is_root())
    {
        cout << "Lattice vectors (Bohr)" << endl;
        latvec.print(16);
        cout << "Reciprocal lattice vectors (2PI Bohr^-1)" << endl;
        G.print(16);
        lib_printf("kgrids: %3d %3d %3d\n", kv_nmp[0], kv_nmp[1], kv_nmp[2]);
        cout << "k-points read (Cartisian in 2Pi Bohr^-1 | fractional):" << endl;
        for (int ik = 0; ik != meanfield.get_n_kpoints(); ik++)
        {
            lib_printf("ik %4d: %12.7f %12.7f %12.7f | %12.7f %12.7f %12.7f\n",
                   ik+1, kvec_c[ik].x, kvec_c[ik].y, kvec_c[ik].z,
                   kfrac_list[ik].x, kfrac_list[ik].y, kfrac_list[ik].z);
        }
        cout << "R-points to compute:" << endl;
        for (int iR = 0; iR != Rlist.size(); iR++)
        {
            lib_printf("%4d: %3d %3d %3d\n", iR+1, Rlist[iR].x, Rlist[iR].y, Rlist[iR].z);
        }
        cout << endl;
    }
    mpi_comm_global_h.barrier();

    read_eigenvector("./", meanfield);
    get_natom_ncell_from_first_Cs_file(natom, ncell, "./", Params::binary_input);
    tot_atpair = generate_atom_pair_from_nat(natom, false);
    tot_atpair_ordered = generate_atom_pair_from_nat(natom, true);
    if (mpi_comm_global_h.is_root())
    {
        cout << "| Number of atoms: " << natom << endl;
        cout << "| Total atom pairs (unordered): " << tot_atpair.size() << endl;
        cout << "| Total atom pairs (ordered)  : " << tot_atpair_ordered.size() << endl;
        cout << "| Number of (R,t) points      : " << Rt_num << endl;
    }

    set_parallel_routing(Params::parallel_routing, tot_atpair.size(), Rt_num, parallel_routing);

    // barrier to wait for information print on master process
    mpi_comm_global_h.barrier();

    //para_mpi.chi_parallel_type=Parallel_MPI::parallel_type::ATOM_PAIR;
    // vector<atpair_t> local_atpair;
    if(parallel_routing == ParallelRouting::ATOM_PAIR)
    {
        // vector<int> atoms_list(natom);
        // for(int iat=0;iat!=natom;iat++)
        //     atoms_list[iat]=iat;
        // std::array<int,3> period_arr{kv_nmp[0], kv_nmp[1], kv_nmp[2]};
        // std::pair<std::vector<int>, std::vector<std::vector<std::pair<int, std::array<int, 3>>>>> list_loc_atp
        // = RI::Distribute_Equally::distribute_atoms(mpi_comm_global_h.comm, atoms_list, period_arr, 2, false);
        // // for(auto &atp:list_loc_atp.first)
        // //     printf("| myid: %d   atp.first: %d\n",mpi_comm_global_h.myid,atp);
        // // for(auto &atps:list_loc_atp.second)
        // //     printf("| myid: %d   atp.second: %d\n",mpi_comm_global_h.myid,atps);
        // std::ofstream ofs("out."+std::to_string(RI::MPI_Wrapper::mpi_get_rank(mpi_comm_global)));
        // for(auto &af:list_loc_atp.first)
		//     ofs<<af<<"  ";
		// ofs<<endl;
        // for( auto &a1 : list_loc_atp.second)
        //     for(auto &a2:a1)
        //         ofs<<a2.first<<"  ("<<a2.second[0]<<", "<<a2.second[1]<<", "<<a2.second[2]<<" )"<<endl;
        if (mpi_comm_global_h.is_root())
            lib_printf("Triangular dispatching of atom pairs\n");
        auto trangular_loc_atpair= dispatch_upper_trangular_tasks(natom,blacs_ctxt_global_h.myid,blacs_ctxt_global_h.nprows,blacs_ctxt_global_h.npcols,blacs_ctxt_global_h.myprow,blacs_ctxt_global_h.mypcol);
        //local_atpair = dispatch_vector(tot_atpair, mpi_comm_global_h.myid, mpi_comm_global_h.nprocs, true);
        for(auto &iap:trangular_loc_atpair)
            local_atpair.push_back(iap);
        read_Cs("./", Params::cs_threshold,local_atpair, Params::binary_input);
        // for(auto &ap:local_atpair)
        //     printf("   |process %d , local_atom_pair:  %d,  %d\n", mpi_comm_global_h.myid,ap.first,ap.second);
        read_Vq_row("./", "coulomb_mat", Params::vq_threshold, local_atpair, false);
        test_libcomm_for_system(Vq);
    }
    else if(parallel_routing == ParallelRouting::LIBRI)
    {
        if (mpi_comm_global_h.is_root()) lib_printf("Evenly distributed Cs and V for LibRI\n");
        read_Cs_evenly_distribute("./", Params::cs_threshold, mpi_comm_global_h.myid, mpi_comm_global_h.nprocs, Params::binary_input);
        // Vq distributed using the same strategy
        // There should be no duplicate for V
        auto trangular_loc_atpair= dispatch_upper_trangular_tasks(natom,blacs_ctxt_global_h.myid,blacs_ctxt_global_h.nprows,blacs_ctxt_global_h.npcols,blacs_ctxt_global_h.myprow,blacs_ctxt_global_h.mypcol);
        for(auto &iap:trangular_loc_atpair)
            local_atpair.push_back(iap);
        read_Vq_row("./", "coulomb_mat", Params::vq_threshold, local_atpair, false);
        test_libcomm_for_system(Vq);
    }
    else
    {
        if (mpi_comm_global_h.is_root()) lib_printf("Complete copy of Cs and V on each process\n");
        local_atpair = generate_atom_pair_from_nat(natom, false);
        read_Cs("./", Params::cs_threshold,local_atpair);
        read_Vq_full("./", "coulomb_mat", false);
    }

    for (int i = 0; i < mpi_comm_global_h.nprocs; i++)
    {
        if (i == mpi_comm_global_h.myid)
        {
            lib_printf("| process %d: Cs size %d from local atpair size %zu\n",
                    mpi_comm_global_h.myid, get_num_keys(Cs), local_atpair.size());
        }
        mpi_comm_global_h.barrier();
    }
    ofs_myid << "Cs size: " << get_num_keys(Cs) << ", with keys:\n";
    print_keys(ofs_myid, Cs);
    std::flush(ofs_myid);

    // debug, check available Coulomb blocks on each process
    // ofs_myid << "Read Coulomb blocks in process\n";
    // for (const auto& IJqcoul: Vq)
    // {
    //     const auto& I = IJqcoul.first;
    //     for (const auto& Jqcoul: IJqcoul.second)
    //     {
    //         const auto& J = Jqcoul.first;
    //         for (const auto& qcoul: Jqcoul.second)
    //         {
    //             ofs_myid << I << " " << J << " " << qcoul.first << "\n";
    //         }
    //     }
    // }
    // std::flush(ofs_myid);
    // mpi_comm_global_h.barrier();
    Profiler::stop("driver_read_common_input_data");

    // malloc_trim(0);
    // para_mpi.mpi_barrier();
    // if(para_mpi.is_master())
    //     system("free -m");
    //erase_Cs_from_local_atp(Cs,local_atpair);
    // malloc_trim(0);
    // para_mpi.mpi_barrier();
    // if(para_mpi.is_master())
    //     system("free -m");
    Chi0 chi0(meanfield, klist, Params::nfreq);
    chi0.gf_R_threshold = Params::gf_R_threshold;

    // build ABF IJ and qlist from Vq

    vector<Vector3_Order<double>> qlist;

    for ( auto q_weight: irk_weight)
    {
        qlist.push_back(q_weight.first);
    }

    mpi_comm_global_h.barrier(); // FIXME: barrier seems not work here...
    if (mpi_comm_global_h.myid == 0)
    {
        cout << "Initialization finished, start task job from myid\n";
    }

    if (task == task_t::RPA)
    {
        task_rpa();
        finalize();
        return 0;
    }

    if ( task != task_t::EXX )
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
    if (task == task_t::RPA)
    {
        Cs.clear();
    }

    LIBRPA::utils::release_free_mem();

    // RPA total energy
    if (task == task_t::RPA)
    {
        mpi_comm_global_h.barrier();
        Profiler::start("EcRPA", "Compute RPA correlation Energy");
        CorrEnergy corr;
        if (Params::use_scalapack_ecrpa && (LIBRPA::parallel_routing == LIBRPA::ParallelRouting::ATOM_PAIR || LIBRPA::parallel_routing == LIBRPA::ParallelRouting::LIBRI))
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
    else if (task == task_t::Wc_Rf)
    {
        Profiler::start("Wc_Rf", "Build Screened Coulomb: R and freq. space");

        Profiler::start("read_vq_cut", "Load truncated Coulomb");
        read_Vq_full("./", "coulomb_cut_", true);
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
                        epsmac_LF_imagfreq_re = LIBRPA::utils::interp_cubic_spline(omegas_dielect, dielect_func,
                                                                           chi0.tfg.get_freq_nodes());
                        break;
                    }
                case 2:
                    {
                        LIBRPA::utils::LevMarqFitting levmarq;
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
                if (mpi_comm_global_h.is_root())
                {
                    lib_printf("Dielection function parsed:\n");
                    for (int i = 0; i < chi0.tfg.get_freq_nodes().size(); i++)
                        lib_printf("%d %f %f\n", i+1, chi0.tfg.get_freq_nodes()[i], epsmac_LF_imagfreq_re[i]);
                }
                mpi_comm_global_h.barrier();
            }
        }

        Profiler::start("compute_Wc_freq_q", "Build upper half of Wc(q,w)");
        vector<std::complex<double>> epsmac_LF_imagfreq(epsmac_LF_imagfreq_re.cbegin(), epsmac_LF_imagfreq_re.cend());
        map<double, atom_mapping<std::map<Vector3_Order<double>, matrix_m<complex<double>>>>::pair_t_old> Wc_freq_q;
        if (Params::use_scalapack_gw_wc)
            Wc_freq_q = compute_Wc_freq_q_blacs(chi0, Vq, Vq_cut, epsmac_LF_imagfreq);
        else
            Wc_freq_q = compute_Wc_freq_q(chi0, Vq, Vq_cut, epsmac_LF_imagfreq);
        Profiler::stop("compute_Wc_freq_q");

        Profiler::start("construct_Wc_lower_half", "Construct Lower Half of Wc(q,w)");
        // NOTE: only upper half of Wc is built now
        //       here we recover the other half before transform to R space using the Hermitian property
        for (int ifreq = 0; ifreq < chi0.tfg.get_n_grids(); ifreq++)
        {
            const auto freq = chi0.tfg.get_freq_nodes()[ifreq];
            auto &Wc = Wc_freq_q.at(freq);
            vector<atom_t> iatoms_row;
            for (const auto &Wc_IJ: Wc)
                iatoms_row.push_back(Wc_IJ.first);
            for (auto iatom_row: iatoms_row)
            {
                vector<atom_t> iatoms_col;
                for (const auto &Wc_J: Wc.at(iatom_row))
                {
                    iatoms_col.push_back(Wc_J.first);
                }
                for (auto iatom_col: iatoms_col)
                {
                    if (iatom_row == iatom_col) continue;
                    for (const auto &Wc_q: Wc.at(iatom_row).at(iatom_col))
                    {
                        Wc[iatom_col][iatom_row][Wc_q.first] = conj(Wc_q.second);
                    }
                }
            }
        }
        Profiler::stop("construct_Wc_lower_half");

        Profiler::start("FT_Wc_freq_q", "Fourier Transform Wc(q,w) -> Wc(R,w)");
        const auto Wc_freq_MN_R = FT_Wc_freq_q(Wc_freq_q, chi0.tfg, Rlist);
        Profiler::stop("FT_Wc_freq_q");

        Profiler::start("write_Wc_freq_R", "Export Wc(R,w) to file");
        for (const auto &freq_MuNuRWc: Wc_freq_MN_R)
        {
            char fn[80];
            auto freq = freq_MuNuRWc.first;
            auto ifreq = chi0.tfg.get_freq_index(freq);
            for (const auto & Mu_NuRWc: freq_MuNuRWc.second)
            {
                auto Mu = Mu_NuRWc.first;
                // const int n_mu = atom_mu[Mu];
                for (const auto & Nu_RWc: Mu_NuRWc.second)
                {
                    auto Nu = Nu_RWc.first;
                    // const int n_nu = atom_mu[Nu];
                    for (const auto & R_Wc: Nu_RWc.second)
                    {
                        auto R = R_Wc.first;
                        auto Wc = R_Wc.second;
                        auto iteR = std::find(Rlist.cbegin(), Rlist.cend(), R);
                        auto iR = std::distance(Rlist.cbegin(), iteR);
                        sprintf(fn, "Wc_Mu_%zu_Nu_%zu_iR_%zu_ifreq_%d.mtx", Mu, Nu, iR, ifreq);
                        print_matrix_mm_file(Wc, Params::output_dir + "/" + fn, 1e-10);
                    }
                }
            }
        }
        Profiler::stop("write_Wc_freq_R");

        Profiler::stop("Wc_Rf");
    }
    else if ( task == task_t::G0W0 )
    {
        Profiler::start("g0w0", "G0W0 quasi-particle calculation");

        Profiler::start("read_vq_cut", "Load truncated Coulomb");
        read_Vq_full("./", "coulomb_cut_", true);
        const auto VR = FT_Vq(Vq_cut, Rlist, true);
        Profiler::stop("read_vq_cut");

        Profiler::start("read_vxc", "Load DFT xc potential");
        std::vector<matrix> vxc;
        int flag_read_vxc = read_vxc("./vxc_out", vxc);
        if (mpi_comm_global_h.myid == 0)
        {
            if (flag_read_vxc == 0)
            {
                cout << "* Success: Read DFT xc potential, will solve quasi-particle equation\n";
            }
            else
            {
                cout << "*   Error: Read DFT xc potential, switch off solving quasi-particle equation\n";
            }
        }
        Profiler::stop("read_vxc");

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
                        epsmac_LF_imagfreq_re = LIBRPA::utils::interp_cubic_spline(omegas_dielect, dielect_func,
                                                                           chi0.tfg.get_freq_nodes());
                        break;
                    }
                case 2:
                    {
                        LIBRPA::utils::LevMarqFitting levmarq;
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
                if (mpi_comm_global_h.is_root())
                {
                    lib_printf("Dielectric function parsed:\n");
                    for (int i = 0; i < chi0.tfg.get_freq_nodes().size(); i++)
                        lib_printf("%d %f %f\n", i+1, chi0.tfg.get_freq_nodes()[i], epsmac_LF_imagfreq_re[i]);
                }
                mpi_comm_global_h.barrier();
            }
        }

        Profiler::start("g0w0_exx", "Build exchange self-energy");
        auto exx = LIBRPA::Exx(meanfield, kfrac_list);
        exx.build_exx_orbital_energy(Cs, Rlist, period, VR);
        Profiler::stop("g0w0_exx");
        if (mpi_comm_global_h.is_root())
        {
            for (int isp = 0; isp != meanfield.get_n_spins(); isp++)
            {
                lib_printf("Spin channel %1d\n", isp+1);
                for (int ik = 0; ik != meanfield.get_n_kpoints(); ik++)
                {
                    cout << "k-point " << ik + 1 << ": " << kvec_c[ik] << endl;
                    lib_printf("%-4s  %-10s  %-10s\n", "Band", "e_exx (Ha)", "e_exx (eV)");
                    for (int ib = 0; ib != meanfield.get_n_bands(); ib++)
                        lib_printf("%4d  %10.5f  %10.5f\n", ib+1, exx.Eexx[isp][ik][ib], HA2EV * exx.Eexx[isp][ik][ib]);
                    lib_printf("\n");
                }
            }
        }
        mpi_comm_global_h.barrier();

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
                            sprintf(fn, "Wcfq_ifreq_%d_iq_%d_I_%zu_J_%zu_id_%d.mtx", ifreq, iq, I, J, mpi_comm_global_h.myid);
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

        if (Params::debug)
        { // debug, check sigc_ij
            char fn[80];
            // ofs_myid << "s_g0w0.sigc_is_f_k_IJ:\n" << s_g0w0.sigc_is_f_k_IJ << "\n";
            for (const auto &is_sigc: s_g0w0.sigc_is_f_k_IJ)
            {
                const int ispin = is_sigc.first;
                for (const auto &f_sigc: is_sigc.second)
                {
                    const auto ifreq = chi0.tfg.get_freq_index(f_sigc.first);
                    for (const auto &k_sigc: f_sigc.second)
                    {
                        const auto &k = k_sigc.first;
                        const int ik = std::distance(klist.begin(), std::find(klist.begin(), klist.end(), k));
                        for (const auto &I_sigc: k_sigc.second)
                        {
                            const auto &I = I_sigc.first;
                            for (const auto &J_sigc: I_sigc.second)
                            {
                                const auto &J = J_sigc.first;
                                sprintf(fn, "Sigcfq_ispin_%d_ifreq_%d_ik_%d_I_%zu_J_%zu_id_%d.mtx",
                                        ispin, ifreq, ik, I, J, mpi_comm_global_h.myid);
                                print_matrix_mm_file(J_sigc.second, Params::output_dir + "/" + fn, 1e-15);
                            }
                        }
                    }
                }
            }
        }

        Profiler::start("g0w0_sigc_rotate_KS", "Rotate self-energy, IJ -> ij -> KS");
        s_g0w0.build_sigc_matrix_KS();
        Profiler::stop("g0w0_sigc_rotate_KS");

        // Solve QPE only when Vxc is successfully parsed
        if (flag_read_vxc == 0)
        {
            Profiler::start("g0w0_solve_qpe", "Solve quasi-particle equation");
            if (mpi_comm_global_h.is_root())
            {
                std::cout << "Solving quasi-particle equation\n";
            }
            std::vector<cplxdb> imagfreqs;
            for (const auto &freq: chi0.tfg.get_freq_nodes())
            {
                imagfreqs.push_back(cplxdb{0.0, freq});
            }

            // TODO: parallelize analytic continuation and QPE solver among tasks
            if (mpi_comm_global_h.is_root())
            {
                map<int, map<int, map<int, double>>> e_qp_all;
                map<int, map<int, map<int, cplxdb>>> sigc_all;
                const auto efermi = meanfield.get_efermi() * 0.5;
                for (int i_spin = 0; i_spin < meanfield.get_n_spins(); i_spin++)
                {
                    for (int i_kpoint = 0; i_kpoint < meanfield.get_n_kpoints(); i_kpoint++)
                    {
                        const auto &sigc_sk = s_g0w0.sigc_is_ik_f_KS[i_spin][i_kpoint];
                        for (int i_state = 0; i_state < meanfield.get_n_bands(); i_state++)
                        {
                            // convert Ry to Ha
                            // FIXME: we should use a consistent internal unit!
                            const auto &eks_state = meanfield.get_eigenvals()[i_spin](i_kpoint, i_state) * 0.5;
                            const auto &exx_state = exx.Eexx[i_spin][i_kpoint][i_state];
                            const auto &vxc_state = vxc[i_spin](i_kpoint, i_state);
                            std::vector<cplxdb> sigc_state;
                            for (const auto &freq: chi0.tfg.get_freq_nodes())
                            {
                                sigc_state.push_back(sigc_sk.at(freq)(i_state, i_state));
                            }
                            LIBRPA::AnalyContPade pade(Params::n_params_anacon, imagfreqs, sigc_state);
                            double e_qp;
                            cplxdb sigc;
                            int flag_qpe_solver = LIBRPA::qpe_solver_pade_self_consistent(
                                pade, eks_state, efermi, vxc_state, exx_state, e_qp, sigc);
                            if (flag_qpe_solver == 0)
                            {
                                e_qp_all[i_spin][i_kpoint][i_state] = e_qp;
                                sigc_all[i_spin][i_kpoint][i_state] = sigc;
                            }
                            else
                            {
                                printf("Warning! QPE solver failed for spin %d, kpoint %d, state %d\n",
                                        i_spin+1, i_kpoint+1, i_state+1);
                                e_qp_all[i_spin][i_kpoint][i_state] = std::numeric_limits<double>::quiet_NaN();
                                sigc_all[i_spin][i_kpoint][i_state] = std::numeric_limits<cplxdb>::quiet_NaN();
                            }
                        }
                    }
                }

                // display results
                const std::string banner(90, '-');
                printf("Printing quasi-particle energy [unit: eV]\n\n");
                for (int i_spin = 0; i_spin < meanfield.get_n_spins(); i_spin++)
                {
                    for (int i_kpoint = 0; i_kpoint < meanfield.get_n_kpoints(); i_kpoint++)
                    {
                        const auto &k = kfrac_list[i_kpoint];
                        printf("spin %2d, k-point %4d: (%.5f, %.5f, %.5f) \n",
                                i_spin+1, i_kpoint+1, k.x, k.y, k.z);
                        printf("%60s\n", banner.c_str());
                        printf("%5s %16s %16s %16s %16s %16s\n", "State", "e_mf", "v_xc", "v_exx", "ReSigc", "e_qp");
                        printf("%60s\n", banner.c_str());
                        for (int i_state = 0; i_state < meanfield.get_n_bands(); i_state++)
                        {
                            const auto &eks_state = meanfield.get_eigenvals()[i_spin](i_kpoint, i_state) * RY2EV;
                            const auto &exx_state = exx.Eexx[i_spin][i_kpoint][i_state] * HA2EV;
                            const auto &vxc_state = vxc[i_spin](i_kpoint, i_state) * HA2EV;
                            const auto &resigc = sigc_all[i_spin][i_kpoint][i_state].real() * HA2EV;
                            const auto &eqp = e_qp_all[i_spin][i_kpoint][i_state] * HA2EV;
                            printf("%5d %16.5f %16.5f %16.5f %16.5f %16.5f\n",
                                   i_state+1, eks_state, vxc_state, exx_state, resigc, eqp);
                        }
                        printf("\n");
                    }
                }
            }
            Profiler::stop("g0w0_solve_qpe");
        }

        Profiler::start("g0w0_export_sigc_KS", "Export self-energy in KS basis");
        mpi_comm_global_h.barrier();
        if (mpi_comm_global_h.is_root())
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
    else if ( task == task_t::EXX )
    {
        read_Vq_full("./", "coulomb_cut_", true);
        const auto VR = FT_Vq(Vq_cut, Rlist, true);
        auto exx = LIBRPA::Exx(meanfield, kfrac_list);
        exx.build_exx_orbital_energy(Cs, Rlist, period, VR);
        // FIXME: need to reduce first when MPI is used
        // NOTE: may extract to a common function
        if (mpi_comm_global_h.is_root())
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

    finalize();
    return 0;
}
