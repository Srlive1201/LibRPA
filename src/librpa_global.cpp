#include "../include/librpa_global.h"

void librpa_init_global_env(LibrpaSwitch switch_redirect_stdout, const char *redirect_path)
{
}

void librpa_finalize_global_env()
{
}


namespace librpa
{

void init_global_env(LibrpaSwitch switch_redirect_stdout, const char *redirect_path)
{
    ::librpa_init_global_env(switch_redirect_stdout, redirect_path);
}

void finalize_global_env()
{
    ::librpa_finalize_global_env();
}

}
