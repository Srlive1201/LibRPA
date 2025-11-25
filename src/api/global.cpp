#include "librpa_global.h"
#include "librpa_enums.h"

#include "../mpi/global_mpi.h"
#include "../io/global_io.h"
#include "../version.h"

inline int librpa_get_major_version(void) { return LIBRPA_MAJOR_VERSION; }

inline int librpa_get_minor_version(void) { return LIBRPA_MINOR_VERSION; }

inline int librpa_get_micro_version(void) { return LIBRPA_MICRO_VERSION; }

void librpa_init_global(LibrpaSwitch switch_redirect_stdout, const char *redirect_path)
{
    librpa_int::global::init_global_mpi();
    const bool redirect_stdout = switch_redirect_stdout == LIBRPA_SWITCH_ON;
    librpa_int::global::init_global_io(redirect_stdout, redirect_path);
}

void librpa_finalize_global(void)
{
    librpa_int::global::finalize_global_io();
    librpa_int::global::finalize_global_mpi();
}


namespace librpa
{

void init_global_env(LibrpaSwitch switch_redirect_stdout, const char *redirect_path)
{
    ::librpa_init_global(switch_redirect_stdout, redirect_path);
}

void finalize_global_env(void)
{
    ::librpa_finalize_global();
}

}
