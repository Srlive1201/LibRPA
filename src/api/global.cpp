#include "../../include/librpa_global.h"

#include "../mpi/global_mpi.h"

void librpa_init_global(LibrpaSwitch switch_redirect_stdout, const char *redirect_path)
{
    librpa_int::global::init_global_mpi();
}

void librpa_finalize_global()
{
    librpa_int::global::finalize_global_mpi();
}


namespace librpa
{

void init_global_env(LibrpaSwitch switch_redirect_stdout, const char *redirect_path)
{
    ::librpa_init_global(switch_redirect_stdout, redirect_path);
}

void finalize_global_env()
{
    ::librpa_finalize_global();
}

}
