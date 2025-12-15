#include "global.h"

// Public API headers
#include "librpa_global.h"
#include "librpa_enums.h"

// Internal headers
#include "../utils/utils_cmake.h"
#include "../version.h"

const char* librpa_get_build_info(void)
{
    return librpa_int::cmake_info_storage().c_str();
}

int librpa_get_major_version(void) { return LIBRPA_MAJOR_VERSION; }

int librpa_get_minor_version(void) { return LIBRPA_MINOR_VERSION; }

int librpa_get_patch_version(void) { return LIBRPA_PATCH_VERSION; }

void librpa_init_global(LibrpaSwitch switch_redirect_stdout,
                        const char *redirect_path,
                        LibrpaSwitch switch_process_output)
{
    using namespace librpa_int::global;

    init_global_mpi();
    const bool redirect_stdout = switch_redirect_stdout == LIBRPA_SWITCH_ON;
    const bool process_output = switch_process_output == LIBRPA_SWITCH_ON;
    librpa_int::global::init_global_io(redirect_stdout, redirect_path, process_output);

    mpi_comm_global_h.barrier();
    lib_printf_root("Initialized LibRPA global environment\n");
    mpi_comm_global_h.barrier();
}

void librpa_finalize_global(void)
{
    using namespace librpa_int::global;
    mpi_comm_global_h.barrier();
    lib_printf_root("Finalizing LibRPA global environment\n");
    mpi_comm_global_h.barrier();

    finalize_global_io();
    finalize_global_mpi();
}

void librpa_test(void)
{
}
