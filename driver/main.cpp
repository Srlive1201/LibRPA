#include "librpa.hpp"

#include "driver.h"
#include "read_data.h"
#include "inputfile.h"
#include "task.h"

#include <mpi.h>
#include <omp.h>

#include <filesystem>

// Internal headers, used here only for printing formation and some consistency check
// May move to public API later
#include "../src/utils/profiler.h"
#include "../src/utils/utils_mem.h"
#include "../src/io/fs.h"
// #include "task_qsgw.h"
// #include "task_qsgwA.h"
// #include "task_qsgw_band.h"
// #include "task_hf_band.h"
// #include "task_scRPA.h"
// #include "task_scRPA_band.h"

static void initialize(int argc, char **argv)
{
    using namespace librpa_int::global;

    // MPI Initialization
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
    if (MPI_THREAD_MULTIPLE != provided)
    {
        librpa_int::global::lib_printf("Warning: MPI_Init_thread provide %d != required %d", provided, MPI_THREAD_MULTIPLE);
    }

    librpa::init_global(LIBRPA_SWITCH_OFF);

    // Global profiler begins right after MPI is initialized
    profiler.start("driver_total", "Total for driver");

    lib_printf_root("Total number of tasks    : %5d\n", size_global);
    lib_printf_root("Total number of nodes    : %5d\n", size_inter);
    lib_printf_root("Maximal number of threads: %3d\n", omp_get_max_threads());
    mpi_comm_global_h.barrier();
    lib_printf_root("MPI tasks information:\n");
    mpi_comm_global_h.barrier();
    lib_printf_coll("| %s\n", mpi_comm_global_h.str().c_str());
    mpi_comm_global_h.barrier();
    printf_comm_root(mpi_comm_intra_h, "| Global ID of master process of node %5d : %5d\n",
                     mpi_comm_inter_h.myid, mpi_comm_global_h.myid);
    mpi_comm_global_h.barrier();

    // Print cmake infomation
    if (mpi_comm_global_h.is_root())
    {
        lib_printf("\n");
        lib_printf("%s", librpa::get_build_info());
        lib_printf("\n");
    }
    mpi_comm_global_h.barrier();

    // Create the handler before parsing any data or computation
    driver::h.init(MPI_COMM_WORLD);
}

static void finalize(bool success)
{
    using namespace librpa_int::global;

    // Free the memory space
    driver::h.free();
    profiler.stop("driver_total");

    bool is_root = mpi_comm_global_h.myid == 0;
    mpi_comm_global_h.barrier();
    librpa::finalize_global();

    if (is_root)
    {
        librpa::print_profile();
        if (success)
        {
            printf("libRPA finished successfully\n");
        }
        else
        {
            printf("libRPA failed\n");
        }
    }

    MPI_Finalize();
}

int main(int argc, char **argv)
{
    using namespace driver;
    using namespace librpa_int::global;
    using librpa_int::get_node_free_mem;

    initialize(argc, argv);

    // Parse input file with runtime options
    profiler.start("driver_read_params", "Driver Read Input Parameters");
    parse_inputfile_to_params(input_filename);
    if (mpi_comm_global_h.is_root())
    {
        lib_printf("===== Begin driver parameters  =====\n");
        lib_printf(driver_params.format().c_str());
        lib_printf("===== End driver parameters    =====\n\n");
        lib_printf("===== Begin control parameters =====\n");
        lib_printf(format_runtime_options(opts).c_str());
        lib_printf("===== End control parameters   =====\n\n");
        std::filesystem::create_directories(driver::opts.output_dir);
    }
    mpi_comm_global_h.barrier();
    profiler.stop("driver_read_params");

    const string path_stru = driver_params.input_dir + "stru_out";
    const string path_bz_sampling = driver_params.input_dir + "bz_sampling_out";
    const string path_basis = driver_params.input_dir + "basis_out";
    const string path_scf_eigen = driver_params.input_dir + "band_out";

    profiler.start("driver_read_common_input_data", "Driver Read Task-Common Input Data");
    profiler.start("driver_band_out", "DFT SCF eigenvalues/occupations");
    read_scf_occ_eigenvalues(path_scf_eigen);
    profiler.stop("driver_band_out");

    task_t task = get_task(driver_params.task);

    std::map<Vector3_Order<double>, ComplexMatrix> sinvS;

    if (task != task_t::print_minimax)
    {
        profiler.start("driver_struct", "Structure");
        read_stru(path_stru);
        profiler.stop("driver_struct");
        lib_printf_root("\n");

        profiler.start("driver_bz", "BZ sampling");
        if (librpa_int::path_exists(path_bz_sampling.c_str()))
            read_bz_sampling(path_bz_sampling);
        else
            read_bz_sampling_from_stru(path_stru);
        lib_printf_root("\n");
        profiler.stop("driver_bz");

        profiler.start("driver_basis", "Basis (wave-function and auxiliary)");
        if (librpa_int::path_exists(path_basis.c_str()))
            read_basis(path_basis);
        else
            read_basis_from_Cs(driver_params.input_dir);
        lib_printf_root("\n");
        profiler.stop("driver_basis");

        profiler.start("driver_read_eigenvector", "SCF eigenvectors");
        int ret_eigenvec = read_eigenvector(driver_params.input_dir);
        mpi_comm_global_h.barrier();
        if (ret_eigenvec == 0)
        {
            lib_printf_root("Successfully read eigenvector files\n");
        }
        else
        {
            if (ret_eigenvec > 0)
            {
                lib_printf_root("Error in reading eigenvector files (retcode %d)\n", ret_eigenvec);
            }
            else
            {
                lib_printf_root(
                    "Error!!! No eigenvector files is found at directory, check if you "
                    "have input files KS_eigenvector\n");
            }
            finalize(false);
            return EXIT_FAILURE;
        }
        profiler.stop("driver_read_eigenvector");

        profiler.start("driver_read_ri");
        read_ri(driver_params.input_dir, driver::opts.parallel_routing);
        lib_printf_root("Actual parallel routing used: %s\n", get_routing_string(driver::opts.parallel_routing).c_str());
        profiler.stop("driver_read_ri");

        const bool use_shrink_abfs = false;
        // TODO: need to include the shrinked basis information in another AtomicBasis object
        if (use_shrink_abfs)
        {
            if (mpi_comm_global_h.is_root())
            {
                std::cout << "iatom & large Nabfs: " << std::endl;
                for (auto &Imu : atom_mu)
                {
                    auto I = Imu.first;
                    auto mu = Imu.second;
                    std::cout << I << "," << mu << std::endl;
                }
            }
            // backup large atom_mu
            // atom_mu_l = atom_mu;  // TODO: replace with the actual shrinked ABFs
            read_Cs_evenly_distribute(driver_params.input_dir, driver_params.cs_threshold,
                                      mpi_comm_global_h.myid, mpi_comm_global_h.nprocs,
                                      "Cs_shrinked_data");

            profiler.start("read_shrink_sinvS_fold", "Load shrink transformation");
            // change atom_mu: number of {Mu,mu} in the later calculations
            read_shrink_sinvS(driver_params.input_dir, "shrink_sinvS_", sinvS);
            atom_mu_s = atom_mu;
            if (mpi_comm_global_h.is_root())
            {
                std::cout << "iatom & small Nabfs: " << std::endl;
                for (auto &Imu : atom_mu)
                {
                    auto I = Imu.first;
                    auto mu = Imu.second;
                    std::cout << I << "," << mu << std::endl;
                }
            }
            profiler.stop("read_shrink_sinvS");
        }
        // Vq distributed using the same strategy
        // There should be no duplicate for V
    }

    mpi_comm_global_h.barrier();
    if (librpa_int::global::mpi_comm_intra_h.myid == 0)
    {
        if (mpi_comm_global_h.myid == 0)
        {
            const auto cputime = profiler.get_cpu_time_last("driver_read_common_input_data") / 60.0;
            const auto walltime = profiler.get_wall_time_last("driver_read_common_input_data") / 60.0;
            lib_printf("Initialization finished, Wall/CPU time [min]: %12.4f %12.4f\n", walltime, cputime);
        }
        double freemem;
        auto flag = get_node_free_mem(freemem);
        if (flag == 0)
        {
            lib_printf("Free memory on node %5d [GB]: %8.3f\n", mpi_comm_inter_h.myid, freemem);
        }
    }
    lib_printf_root("Common data parsed, task %s will begin\n", driver_params.task.c_str());
    lib_printf_root("%s: %s\n", driver_params.task.c_str(), get_task_string(task).c_str());
    mpi_comm_global_h.barrier();
    profiler.stop("driver_read_common_input_data");

    run_task(task);
    mpi_comm_global_h.barrier();

    finalize(true);
    return EXIT_SUCCESS;
}
