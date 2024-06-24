#include "envs_mpi.h"

#include <stdexcept>

namespace LIBRPA
{

namespace envs
{

int myid_global = 0;

int size_global = 1;

MPI_Comm mpi_comm_global;

std::string procname = "localhost";

static bool librpa_mpi_initialized = false;

void initialize_mpi(MPI_Comm mpi_comm_global_in)
{
    int flag;
    MPI_Initialized(&flag);
    if (flag == 0)
    {
        throw std::runtime_error(
            std::string(__FILE__) + ":" + std::to_string(__LINE__) + ":" + std::string(__FUNCTION__) + ": "
            "MPI_Init or MPI_Init_thread must be called before " + std::string(__FUNCTION__));
    }

    mpi_comm_global = mpi_comm_global_in;

    MPI_Comm_rank(mpi_comm_global_in, &myid_global);
    MPI_Comm_size(mpi_comm_global_in, &size_global);

    char name[MPI_MAX_PROCESSOR_NAME];
    int length;
    MPI_Get_processor_name(name, &length);
    procname = name;

    librpa_mpi_initialized = true;
}

bool is_mpi_initialized()
{
    return librpa_mpi_initialized;
}

void finalize_mpi()
{
    librpa_mpi_initialized = false;
}

} /* end of namespace envs */

} /* end of namespace LIBRPA */
