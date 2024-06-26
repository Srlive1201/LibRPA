#include "envs_io.h"

#include <string>
#include <iostream>
#include <cstdio>

#include "envs_mpi.h"


namespace LIBRPA
{

namespace envs
{

std::ofstream ofs_myid;
bool redirect_stdout = false;
FILE *pfile_redirect;

/* ==============================================
 * Private variables starts
 */
static bool librpa_io_initialized = false;
static std::streambuf *cout_buf_old;
/*
 * Private variables ends
 ============================================== */

void initialize_io(bool redirect_stdout_in, const char *filename)
{
    if (!is_mpi_initialized())
    {
        throw std::runtime_error(
            std::string(__FILE__) + ":" + std::to_string(__LINE__) + ":" + std::string(__FUNCTION__) + ": "
            "initialize_mpi must be called before " + std::string(__FUNCTION__));
    }

    std::string fn = "librpa_para_nprocs_" + std::to_string(size_global) +  "_myid_" + std::to_string(myid_global) + ".out";
    ofs_myid.open(fn);

    redirect_stdout = redirect_stdout_in;

    if (redirect_stdout)
    {
        // Redirect cout stream
        std::ofstream ofs(filename);
        // Save the original cout buffer
        cout_buf_old = std::cout.rdbuf();
        std::cout.rdbuf(ofs.rdbuf());

        // for redirecting printf by using fprintf(pfile_redirect)
        pfile_redirect = fopen(filename, "w");
    }

    librpa_io_initialized = true;
}

bool is_io_initialized()
{
    return librpa_io_initialized;
}

void finalize_io()
{
    ofs_myid.close();

    // Recover original stdout behavior
    if (redirect_stdout)
    {
        std::cout.rdbuf(cout_buf_old);
        fclose(pfile_redirect);
    }

    redirect_stdout = false;
    librpa_io_initialized = false;
}

} /* end of namespace envs */

} /* end of namespace LIBRPA */
