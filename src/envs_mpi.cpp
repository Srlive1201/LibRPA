#include "envs_mpi.h"

#include <stdexcept>

namespace LIBRPA
{

namespace envs
{

int myid_global = 0;

int size_global = 1;

MPI_Comm mpi_comm_global;

LIBRPA::MPI_COMM_handler mpi_comm_global_h;

LIBRPA::BLACS_CTXT_handler blacs_ctxt_global_h;

std::string procname = "not-init";

static bool librpa_mpi_initialized = false;

void initialize_mpi(const MPI_Comm &mpi_comm_global_in)
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

    MPI_Comm_rank(mpi_comm_global, &myid_global);
    MPI_Comm_size(mpi_comm_global, &size_global);

    char name[MPI_MAX_PROCESSOR_NAME];
    int length;
    MPI_Get_processor_name(name, &length);
    procname = name;

    // Initialize handlers of global MPI communicator and BLACS context
    mpi_comm_global_h.reset_comm(mpi_comm_global);
    mpi_comm_global_h.init();

    blacs_ctxt_global_h.reset_comm(mpi_comm_global);
    blacs_ctxt_global_h.init();

    /*
     * HACK: A close-to-square process grid is imposed.
     *       This might subject to performance issue when
     *       the number of processes is not dividable.
     */
    blacs_ctxt_global_h.set_square_grid();

    librpa_mpi_initialized = true;
}

bool is_mpi_initialized()
{
    return librpa_mpi_initialized;
}

void finalize_mpi()
{
    procname = "not-init";
    librpa_mpi_initialized = false;
}

} /* end of namespace envs */

} /* end of namespace LIBRPA */
