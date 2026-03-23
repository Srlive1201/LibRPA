#include "global.h"

#include <iostream>
#ifdef LIBRPA_USE_LIBRI
#include <RI/physics/Exx.h>
#endif

// Public API headers
#include "librpa_global.h"
#include "librpa_enums.h"

// Internal headers
// #include "../core/exx.h"
// #include "../utils/constants.h"
#include "../mpi/global_mpi.h"
#include "../io/global_io.h"
#include "../utils/profiler.h"
#include "../utils/utils_cmake.h"
#include "../version.h"

const char* librpa_get_build_info(void)
{
    return librpa_int::cmake_info_storage().c_str();
}

int librpa_get_major_version(void) { return LIBRPA_MAJOR_VERSION; }

int librpa_get_minor_version(void) { return LIBRPA_MINOR_VERSION; }

int librpa_get_patch_version(void) { return LIBRPA_PATCH_VERSION; }

static void librpa_init_global_common(MPI_Comm comm, LibrpaSwitch switch_redirect_stdout,
                                      const char *redirect_path, LibrpaSwitch switch_process_output)
{
    using namespace librpa_int::global;

    init_global_mpi(comm);
    const bool redirect_stdout = switch_redirect_stdout == LIBRPA_SWITCH_ON;
    const bool process_output = switch_process_output == LIBRPA_SWITCH_ON;
    librpa_int::global::init_global_io(redirect_stdout, redirect_path, process_output);

    mpi_comm_global_h.barrier();
    lib_printf_root("Initialized LibRPA global environment\n");
    mpi_comm_global_h.barrier();
}

void librpa_init_global(LibrpaSwitch switch_redirect_stdout,
                        const char *redirect_path,
                        LibrpaSwitch switch_process_output)
{
    ::librpa_init_global_common(MPI_COMM_WORLD, switch_redirect_stdout, redirect_path, switch_process_output);
}

// Only for Fortran binding
void librpa_init_global_fortran(MPI_Fint *f_comm, LibrpaSwitch switch_redirect_stdout,
                                const char *redirect_path, LibrpaSwitch switch_process_output)
{
    MPI_Comm c_comm = MPI_Comm_f2c(*f_comm);
    ::librpa_init_global_common(c_comm, switch_redirect_stdout, redirect_path, switch_process_output);
}

void librpa_finalize_global(void)
{
    using namespace librpa_int::global;
    mpi_comm_global_h.barrier();
    lib_printf_root("Finalizing LibRPA global environment\n");
    mpi_comm_global_h.barrier();

    // print per-process profiling
    if (ofs_myid.is_open())
    {
        ofs_myid << "Print per-process profiling data" << std::endl;
        ofs_myid << profiler.get_profile_string();
    }

    finalize_global_io();
    finalize_global_mpi();
}

// static void test_bccHe_libri_exx()
// {
//     std::cout << "Hello " << __FUNCTION__ << std::endl;
// #ifdef LIBRPA_USE_LIBRI
//     RI::Exx<int, int, 3, double> exx_libri;
//     // int flag;
//     // MPI_Initialized(&flag);
//     // std::cout << " " << flag;
//     std::cout << std::endl;
//     std::map<int,std::array<double,3>> atoms_pos;
//     const int n_atoms = 2;
//     for (int i = 0; i < n_atoms; i++)
//         atoms_pos.insert(std::pair<int, std::array<double, 3>>{i, {0, 0, 0}});
//     std::array<std::array<double, 3>, 3> latvec_array;
//     latvec_array[0] = {5.7, 0.0, 0.0};
//     latvec_array[1] = {0.0, 5.7, 0.0};
//     latvec_array[2] = {0.0, 0.0, 5.7};
//     exx_libri.set_parallel(MPI_COMM_WORLD, atoms_pos, latvec_array, {2, 2, 2});
//     exx_libri.set_Cs({}, 0.0);
//     exx_libri.set_Vs({}, 0.0);
//     exx_libri.set_Ds({}, 0.0);
// #endif
// }

// static void test_bccHe_exx()
// {
//     using namespace librpa_int;
//
//     std::cout << "Hello " << __FUNCTION__ << std::endl;
//     librpa_int::MeanField mf(1, 8, 8, 8);
//     librpa_int::AtomicBasis wfc(std::vector<size_t>{4,4});
//     librpa_int::AtomicBasis aux(std::vector<size_t>{13,13});
//     librpa_int::PeriodicBoundaryData pbc;
//     pbc.set_latvec({TWO_PI, 0.0, 0.0, 0.0, TWO_PI, 0.0, 0.0, 0.0, TWO_PI});
//     pbc.set_kgrids_kvec(2, 2, 2,
//                         {
//                             0.0, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.5, 0.0, 0.0, 0.5, 0.5,
//                             0.5, 0.0, 0.0, 0.5, 0.0, 0.5, 0.5, 0.5, 0.0, 0.5, 0.5, 0.5,
//                         });
//     librpa_int::Exx exx(mf, wfc, pbc, global::mpi_comm_global_h, false);
//     exx.build(LibrpaParallelRouting::LIBRI, aux, {}, {});
// }

void librpa_test(void)
{
    std::cout << "Hello " << __FUNCTION__ << std::endl;
    // test_bccHe_libri_exx();
    // test_bccHe_exx();
}

void librpa_print_profile()
{
    librpa_int::global::profiler.display();
}
