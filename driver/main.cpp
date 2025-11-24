#include <algorithm>
#include <omp.h>

#include "../src/core/atomic_basis.h"
#include "../src/io/global_io.h"
#include "../src/mpi/global_mpi.h"
#include "../src/mpi/envs_blacs.h"
#include "../src/core/timefreq.h"
#include "../src/api/librpa.h"
#include "../src/core/geometry.h"
#include "inputfile.h"
#include "../src/core/meanfield.h"
#include "../src/core/params.h"
#include "../src/core/pbc.h"
#include "../src/utils/profiler.h"
#include "read_data.h"
#include "../src/io/stl_io_helper.h"
#include "driver_params.h"
#include "task.h"
#include "../src/utils/utils_cmake.h"
#include "../src/utils/utils_mem.h"

#include "task_rpa.h"
#include "task_exx.h"
#include "task_exx_band.h"
#include "task_gw.h"
#include "task_gw_band.h"
#include "task_screened_coulomb.h"
#include "task_test.h"

static void initialize(int argc, char **argv)
{
    using namespace librpa_int::envs;
    using librpa_int::global::lib_printf;
    // MPI Initialization
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if (MPI_THREAD_MULTIPLE != provided)
    {
        lib_printf("Warning: MPI_Init_thread provide %d != required %d", provided, MPI_THREAD_MULTIPLE);
    }

    initialize_librpa_environment(0, "");

    // Global profiler begins right after MPI is initialized
    Profiler::start("total", "Total");
}

static void finalize(bool success)
{
    using namespace librpa_int::envs;
    using namespace librpa_int::global;
    using librpa_int::global::lib_printf;

    Profiler::stop("total");
    if (mpi_comm_global_h.is_root())
    {
        if (success)
        {
            Profiler::display();
            lib_printf("libRPA finished successfully\n");
        }
        else
        {
            lib_printf("libRPA failed\n");
        }
    }

    finalize_librpa_environment();

    MPI_Finalize();
}

int main(int argc, char **argv)
{
    using librpa_int::global::mpi_comm_global_h;
    using librpa_int::global::mpi_comm_global;
    using librpa_int::envs::blacs_ctxt_global_h;
    using librpa_int::ParallelRouting;
    using librpa_int::parallel_routing;
    using librpa_int::task_t;
    using librpa_int::global::ofs_myid;
    using librpa_int::global::lib_printf;
    using librpa_int::global::lib_printf_coll;
    using librpa_int::global::lib_printf_root;
    using librpa_int::printf_comm_coll;

    initialize(argc, argv);
    lib_printf_root("Total number of tasks: %5d\n", librpa_int::global::size_global);
    lib_printf_root("Total number of nodes: %5d\n", librpa_int::global::size_inter);
    lib_printf_root("Maximumal number of threads: %3d\n", omp_get_max_threads());
    mpi_comm_global_h.barrier();
    printf_comm_coll(librpa_int::global::mpi_comm_intra_h, "Global ID of master process of node %5d : %5d\n",
                     librpa_int::global::mpi_comm_inter_h.myid, librpa_int::global::mpi_comm_global_h.myid);
    mpi_comm_global_h.barrier();
    lib_printf_coll("%s\n", mpi_comm_global_h.str().c_str());
    mpi_comm_global_h.barrier();

    // Print cmake infomation
    if (mpi_comm_global_h.is_root())
    {
        lib_printf("\n");
        librpa_int::print_cmake_info();
        lib_printf("\n");
    }
    mpi_comm_global_h.barrier();

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
    std::string task_lower = driver_params.task;
    transform(task_lower.begin(), task_lower.end(), task_lower.begin(), ::tolower);
    if (task_lower == "rpa")
        task = task_t::RPA;
    else if (task_lower == "g0w0")
        task = task_t::G0W0;
    else if (task_lower == "g0w0_band")
        task = task_t::G0W0_band;
    else if (task_lower == "exx")
        task = task_t::EXX;
    else if (task_lower == "exx_band")
        task = task_t::EXX_band;
    else if (task_lower == "wc_rf")
        task = task_t::Wc_Rf;
    else if (task_lower == "print_minimax")
        task = task_t::print_minimax;
    else if (task_lower == "test")
        task = task_t::test;
    else
        throw std::logic_error("Unknown task (" + driver_params.task + "). Please check your input");


    if (mpi_comm_global_h.is_root())
    {
        system(("mkdir -p " + Params::output_dir).c_str());
        lib_printf("===== Begin driver parameters  =====\n");
        driver_params.print();
        lib_printf("===== End driver parameters    =====\n");
        lib_printf("===== Begin control parameters =====\n");
        Params::print();
        lib_printf("===== End control parameters   =====\n");
    }
    mpi_comm_global_h.barrier();
    Profiler::stop("driver_read_params");

    Profiler::start("driver_band_out", "Driver Read Meanfield band");
    read_scf_occ_eigenvalues(driver_params.input_dir + "band_out", meanfield);
    librpa_int::envs::initialize_array_desc_wfc_global(meanfield.get_n_aos());
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

    // early exit for print_minimax task
    if (task == task_t::print_minimax)
    {
        double emin, emax;
        meanfield.get_E_min_max(emin, emax);
        TFGrids tfg(Params::nfreq);
        tfg.generate_minimax(emin, emax);
        if (mpi_comm_global_h.is_root())
            tfg.show();
        finalize(true);
        return 0;
    }

    Profiler::start("driver_read_common_input_data", "Driver Read Task-Common Input Data");
    read_stru(meanfield.get_n_kpoints(), driver_params.input_dir + "stru_out");
    Vector3_Order<int> period {kv_nmp[0], kv_nmp[1], kv_nmp[2]};
    auto Rlist = construct_R_grid(period);
    const int Rt_num = Rlist.size() * Params::nfreq;
    if (mpi_comm_global_h.is_root())
    {
        cout << "Lattice vectors (Bohr)" << endl;
        latvec.print(16);
        cout << "Reciprocal lattice vectors (2PI Bohr^-1)" << endl;
        G.print(16);
        cout << "Atom positions read (Cartisian in Bohr | fractional):" << endl;
        for (int i_at = 0; i_at != coord.size(); i_at++)
        {
            lib_printf("ia %4d: %12.7f %12.7f %12.7f | %12.7f %12.7f %12.7f\n",
                   i_at+1, coord[i_at][0], coord[i_at][1], coord[i_at][2],
                   coord_frac[i_at][0], coord_frac[i_at][1], coord_frac[i_at][2]);
        }
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

    Profiler::start("driver_read_eigenvector");
    int ret_eigenvec = read_eigenvector(driver_params.input_dir, meanfield);
    if (ret_eigenvec == 0)
    {
        lib_printf("Successfully read eigenvector files\n");
    }
    else
    {
        if (ret_eigenvec > 0)
        {
            lib_printf("Error in reading eigenvector files, return code: %d\n", ret_eigenvec);
        }
        else
        {
            lib_printf(
                "Error!!! No eigenvector files is found at directory, check if you "
                "have input files KS_eigenvector\n");
        }
        finalize(false);
        return EXIT_FAILURE;
    }
    Profiler::stop("driver_read_eigenvector");
    get_natom_ncell_from_first_Cs_file(natom, ncell, driver_params.input_dir);
    tot_atpair = generate_atom_pair_from_nat(natom, false);
    tot_atpair_ordered = generate_atom_pair_from_nat(natom, true);
    if (mpi_comm_global_h.is_root())
    {
        cout << "| Number of atoms: " << natom << endl;
        cout << "| Total atom pairs (unordered): " << tot_atpair.size() << endl;
        cout << "| Total atom pairs (ordered)  : " << tot_atpair_ordered.size() << endl;
        cout << "| Number of (R,t) points      : " << Rt_num << endl;
    }
    std::flush(cout);

    set_parallel_routing(Params::parallel_routing, tot_atpair.size(), Rt_num, parallel_routing);

    // barrier to wait for information print on master process
    mpi_comm_global_h.barrier();

    Profiler::start("driver_read_Cs_Vq");
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
        auto trangular_loc_atpair = librpa_int::dispatch_upper_triangular_tasks(
            natom, blacs_ctxt_global_h.myid, blacs_ctxt_global_h.nprows, blacs_ctxt_global_h.npcols,
            blacs_ctxt_global_h.myprow, blacs_ctxt_global_h.mypcol);
        //local_atpair = dispatch_vector(tot_atpair, mpi_comm_global_h.myid, mpi_comm_global_h.nprocs, true);
        for(auto &iap:trangular_loc_atpair)
            local_atpair.push_back(iap);

        Profiler::start("driver_read_Cs");
        read_Cs(driver_params.input_dir, Params::cs_threshold,local_atpair);
        Profiler::cease("driver_read_Cs");

        // for(auto &ap:local_atpair)
        //     printf("   |process %d , local_atom_pair:  %d,  %d\n", mpi_comm_global_h.myid,ap.first,ap.second);
        Profiler::start("driver_read_Vq");
        read_Vq_row(driver_params.input_dir, "coulomb_mat", Params::vq_threshold, local_atpair, false);
        Profiler::cease("driver_read_Vq");
        // test_libcomm_for_system(Vq);
    }
    else if(parallel_routing == ParallelRouting::LIBRI)
    {
        lib_printf_root("Evenly distributed Cs and V for LibRI\n");
        Profiler::start("driver_read_Cs");
        read_Cs_evenly_distribute(driver_params.input_dir, Params::cs_threshold, mpi_comm_global_h.myid, mpi_comm_global_h.nprocs);
        mpi_comm_global_h.barrier();
        Profiler::cease("driver_read_Cs");
        lib_printf_coll("| Process %5d: Cs with %14zu non-zero keys from local atpair size %7zu. "
                        "Data memory: %10.2f MB. Wall/CPU time [min]: %12.4f %12.4f\n",
                        mpi_comm_global_h.myid, Cs_data.n_keys(), local_atpair.size(),
                        Cs_data.n_data_bytes() * 8.0e-6,
                        Profiler::get_wall_time_last("driver_read_Cs") / 60.0,
                        Profiler::get_cpu_time_last("driver_read_Cs") / 60.0);
        // Vq distributed using the same strategy
        // There should be no duplicate for V

        Profiler::start("driver_read_Vq");
        auto trangular_loc_atpair = librpa_int::dispatch_upper_triangular_tasks(
            natom, blacs_ctxt_global_h.myid, blacs_ctxt_global_h.nprows, blacs_ctxt_global_h.npcols,
            blacs_ctxt_global_h.myprow, blacs_ctxt_global_h.mypcol);
        for(auto &iap:trangular_loc_atpair)
            local_atpair.push_back(iap);
        read_Vq_row(driver_params.input_dir, "coulomb_mat", Params::vq_threshold, local_atpair, false);
        mpi_comm_global_h.barrier();
        Profiler::cease("driver_read_Vq");
        lib_printf_coll("| Process %5d: coulomb_mat read. Wall/CPU time [min]: %12.4f %12.4f\n",
                        mpi_comm_global_h.myid,
                        Profiler::get_wall_time_last("driver_read_Vq") / 60.0,
                        Profiler::get_cpu_time_last("driver_read_Vq") / 60.0);
        // test_libcomm_for_system(Vq);
    }
    else
    {
        lib_printf_root("Complete copy of Cs and V on each process\n");
        local_atpair = generate_atom_pair_from_nat(natom, false);
        Profiler::start("driver_read_Cs");
        read_Cs(driver_params.input_dir, Params::cs_threshold, local_atpair);
        Profiler::cease("driver_read_Cs");

        Profiler::start("driver_read_Vq");
        read_Vq_full(driver_params.input_dir, "coulomb_mat", false);
        Profiler::cease("driver_read_Vq");
    }
    Profiler::stop("driver_read_Cs_Vq");
    // After succesfully loading Cs, number of auxiliary basis functions is now fixed
    librpa_int::envs::initialize_array_desc_abf_global(librpa_int::atomic_basis_abf.nb_total);

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

    // Check memory usage on each node
    mpi_comm_global_h.barrier();
    if (librpa_int::global::mpi_comm_intra_h.myid == 0)
    {
        if (mpi_comm_global_h.myid == 0)
        {
            const auto cputime = Profiler::get_cpu_time_last("driver_read_common_input_data") / 60.0;
            const auto walltime = Profiler::get_wall_time_last("driver_read_common_input_data") / 60.0;
            lib_printf("Initialization finished, Wall/CPU time [min]: %12.4f %12.4f\n", walltime, cputime);
            lib_printf("Task work begins: %s\n", task_lower.c_str());

            double freemem;
            auto flag = librpa_int::utils::get_node_free_mem(freemem);
            if (flag == 0)
            {
                lib_printf("Free memory on node %5d [GB]: %8.3f\n",
                           librpa_int::global::mpi_comm_inter_h.myid, freemem);
            }
        }
    }



    if (task == task_t::RPA)
    {
        task_rpa();
    }
    else if (task == task_t::Wc_Rf)
    {
        task_screened_coulomb_real_freq();
    }
    else if (task == task_t::G0W0)
    {
        task_g0w0();
    }
    else if (task == task_t::G0W0_band)
    {
        task_g0w0_band();
    }
    else if (task == task_t::EXX)
    {
        task_exx();
    }
    else if (task == task_t::EXX_band)
    {
        task_exx_band();
    }
    else if (task == task_t::test)
    {
        task_test();
    }

    finalize(true);
    return EXIT_SUCCESS;
}
