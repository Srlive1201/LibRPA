#include "global_io.h"

#include <cstddef>
#include <cstdio>
#include <iostream>
#include <string>

namespace librpa_int
{

namespace global
{

std::ofstream ofs_myid;
FILE *pfile_redirect = nullptr;

/* ==============================================
 * Private variables starts
 */
static std::ofstream ofs;
static std::streambuf *cout_buf_old = nullptr;
static bool librpa_io_initialized = false;

/*
 * Private variables ends
 ============================================== */

void init_global_io(bool redirect_stdout, const char *filename)
{
    using namespace librpa_int::global;
    if (!is_mpi_initialized())
    {
        throw std::runtime_error(
            std::string(__FILE__) + ":" + std::to_string(__LINE__) + ":" + std::string(__FUNCTION__) + ": "
            "initialize_mpi must be called before " + std::string(__FUNCTION__));
    }

    std::string fn = "librpa_para_nprocs_" + std::to_string(size_global) +  "_myid_" + std::to_string(myid_global) + ".out";
    ofs_myid.open(fn);

    // Redirect output
    if (redirect_stdout)
    {
        // Rank 0 touches the output file
        MPI_Barrier(mpi_comm_global);
        if (myid_global == 0)
        {
            std::ofstream ofs(filename);
            ofs.close();
        }
        MPI_Barrier(mpi_comm_global);

        ofs.open(filename, std::ios::out | std::ios::app);
        cout_buf_old = std::cout.rdbuf();
        std::cout.rdbuf(ofs.rdbuf());

        // redirecting printf
        pfile_redirect = fopen(filename, "a");
    }

    librpa_io_initialized = true;
}

inline bool is_io_initialized() { return librpa_io_initialized; }

void finalize_global_io()
{
    ofs_myid.close();

    // Recover original stdout behavior
    if (pfile_redirect != nullptr)
    {
        std::cout.rdbuf(cout_buf_old);
        fclose(pfile_redirect);
        pfile_redirect = nullptr;
    }
    librpa_io_initialized = false;
}

} /* end of namespace global */

} /* end of namespace librpa_int */
