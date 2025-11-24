#pragma once

#include "../utils/envs_io.h"
#include "../utils/utils_io.h"
#include "envs_mpi.h"

namespace librpa_int
{

namespace utils
{

//! simlar to lib_printf, but only proc 0 of global communicator will dump
template <typename... Args>
void lib_printf_root(const char* format, Args&&... args)
{
    if (envs::myid_global == 0)
    {
        lib_printf(format, std::forward<Args>(args)...);
    }
}

//! Similar to lib_printf_root, but one can specify any communicator
template <typename... Args>
void lib_printf_comm_root(const librpa_int::MpiCommHandler &comm_h, const char* format, Args&&... args)
{
    comm_h.check_initialized();
    if (0 == comm_h.myid)
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

//! Similar to lib_printf_cool, but one can specify any communicator
template <typename... Args>
void lib_printf_comm_coll(const librpa_int::MpiCommHandler &comm_h, const char* format, Args&&... args)
{
    comm_h.check_initialized();
    for (int i = 0; i < comm_h.nprocs; i++)
    {
        if (i == comm_h.myid)
        {
            lib_printf(format, std::forward<Args>(args)...);
        }
        comm_h.barrier();
    }
}

} /* end of name space utils */

} /* end of name space LIBRPA */
