#pragma once

#include "utils_io.h"

#include "parallel_mpi.h"


namespace LIBRPA
{

namespace utils
{


// Parallel print in the order of process ID in communicator
template <typename... Args>
void lib_printf_para(const LIBRPA::MPI_COMM_handler &comm_h, const char* format, Args&&... args)
{
    for (int i = 0; i < comm_h.nprocs; i++)
    {
        if (i == comm_h.myid)
        {
            lib_printf(format, std::forward<Args>(args)...);
        }
        comm_h.barrier();
    }
}


// Parallel print, only the master process
template <typename... Args>
void lib_printf_master(const LIBRPA::MPI_COMM_handler &comm_h, const char* format, Args&&... args)
{
    if (0 == comm_h.myid)
    {
        lib_printf(format, std::forward<Args>(args)...);
    }
}


} /* end of namespace utils */

} /* end of namespace LIBRPA */
