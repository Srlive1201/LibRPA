#pragma once
#include <fstream>

#include "../mpi/global_mpi.h"

namespace librpa_int
{

namespace global
{

//! File output stream handler of each process
extern std::ofstream ofs_myid;

//! File stream used by fprintf when stdout is redirected
extern FILE *pfile_redirect;

//! Initialize the IO environment of LibRPA
/*!
 * @param  [in]  redirect_stdout    Rlag to control to redirect stdout to file, useful when separating LibRPA output from the main program
 * @pram   [in]  redirect_path      Path of file to redirect stdout, used only when redirect_stdout is true.
 */
void init_global_io(bool redirect_stdout = false, const char *filename = "LibRPA_output.txt");

//! Check the IO environment of LibRPA is correctly initialized
bool is_io_initialized();

//! Finalize the IO environment of LibRPA
void finalize_global_io();

//! printf that handles the stdout redirect
template <typename... Args>
void lib_printf(const char* format, Args&&... args)
{
    if (pfile_redirect != nullptr)
    {
        std::fprintf(pfile_redirect, format, std::forward<Args>(args)...);
        std::fflush(pfile_redirect);
    }
    else
    {
        std::printf(format, std::forward<Args>(args)...);
    }
}

// raw string, no formatting at all
inline void lib_printf(const char* s)
{
    if (pfile_redirect)
    {
        std::fprintf(pfile_redirect, "%s", s);
        std::fflush(pfile_redirect);
    }
    else
    {
        std::printf("%s", s);
    }
}

//! simlar to global::printf, but only proc 0 of global communicator will dump
template <typename... Args>
void lib_printf_root(const char* format, Args&&... args)
{
    if (myid_global == 0)
    {
        lib_printf(format, std::forward<Args>(args)...);
    }
}

//! simlar to global::printf, but all processes will print in the order of myid
template <typename... Args>
void lib_printf_coll(const char* format, Args&&... args)
{
    for (int i = 0; i < size_global; i++)
    {
        if (myid_global == i)
        {
            lib_printf(format, std::forward<Args>(args)...);
        }
        MPI_Barrier(mpi_comm_global);
    }
}

} /* end of namespace global */

//! Similar to lib_printf_root, but one can specify any communicator
template <typename... Args>
void printf_comm_root(const librpa_int::MpiCommHandler &comm_h, const char* format, Args&&... args)
{
    comm_h.check_initialized();
    if (0 == comm_h.myid)
    {
        global::lib_printf(format, std::forward<Args>(args)...);
    }
}

//! Similar to lib_printf_cool, but one can specify any communicator
template <typename... Args>
void printf_comm_coll(const librpa_int::MpiCommHandler &comm_h, const char* format, Args&&... args)
{
    comm_h.check_initialized();
    for (int i = 0; i < comm_h.nprocs; i++)
    {
        if (i == comm_h.myid)
        {
            global::lib_printf(format, std::forward<Args>(args)...);
        }
        comm_h.barrier();
    }
}

} /* end of namespace librpa_int */
