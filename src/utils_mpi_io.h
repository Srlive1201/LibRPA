#pragma once

#include "utils_io.h"

#include "envs_io.h"
#include "envs_mpi.h"

namespace LIBRPA
{

namespace utils
{

//! simlar to lib_printf, but only proc 0 will dump
template <typename... Args>
void lib_printf_root(const char* format, Args&&... args)
{
    if (envs::myid_global == 0)
    {
        lib_printf(format, std::forward<Args>(args)...);
    }
}

//! simlar to lib_printf, but all processes will print in the order of myid
template <typename... Args>
void lib_printf_coll(const char* format, Args&&... args)
{
    for (int i = 0; i < envs::size_global; i++)
    {
        if (envs::myid_global == i)
        {
            lib_printf(format, std::forward<Args>(args)...);
        }
        MPI_Barrier(envs::mpi_comm_global);
    }
}

} /* end of name space utils */

} /* end of name space LIBRPA */
