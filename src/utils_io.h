#pragma once

#include "envs_io.h"

namespace LIBRPA
{

namespace utils
{

//! printf that handles the stdout redirect
template <typename... Args>
void lib_printf(const char* format, Args&&... args)
{
    envs::redirect_stdout ?
        fprintf(envs::pfile_redirect, format, std::forward<Args>(args)...) :
        printf(format, std::forward<Args>(args)...);
}

} /* end of name space utils */

} /* end of name space LIBRPA */
