#include "envs_mpi.h"

#include <stdexcept>

namespace LIBRPA
{

namespace envs
{

MPI_Comm mpi_comm_global;
LIBRPA::MPI_COMM_handler mpi_comm_global_h;
int myid_global = 0;
int size_global = 1;

MPI_Comm mpi_comm_intra;
LIBRPA::MPI_COMM_handler mpi_comm_intra_h;
int myid_intra = 0;
int size_intra = 1;

MPI_Comm mpi_comm_inter;
LIBRPA::MPI_COMM_handler mpi_comm_inter_h;
int myid_inter = 0;
int size_inter = 1;

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
    // Initialize handlers of global MPI communicator
    mpi_comm_global_h.reset_comm(mpi_comm_global);
    mpi_comm_global_h.init();
    myid_global = mpi_comm_global_h.myid;
    size_global = mpi_comm_global_h.nprocs;
    procname = mpi_comm_global_h.procname;

    // Communicators for shared-memory operations
    // Set up inter-node and intra-node communicators under the global communicator
    flag = MPI_Comm_split_type(mpi_comm_global, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &mpi_comm_intra);
    if (flag != 0)
    {
        throw std::runtime_error(
            std::string(__FILE__) + ":" + std::to_string(__LINE__) + ":" + std::string(__FUNCTION__) + ": "
            "intra-node communicator creation failed");
    }
    mpi_comm_intra_h.reset_comm(mpi_comm_intra);
    mpi_comm_intra_h.init();
    myid_intra = mpi_comm_intra_h.myid;
    size_intra = mpi_comm_intra_h.nprocs;

    flag = MPI_Comm_split(mpi_comm_global, myid_intra, myid_global, &mpi_comm_inter);
    if (flag != 0)
    {
        throw std::runtime_error(
            std::string(__FILE__) + ":" + std::to_string(__LINE__) + ":" + std::string(__FUNCTION__) + ": "
            "inter-node communicator creation failed");
    }
    mpi_comm_inter_h.reset_comm(mpi_comm_inter);
    mpi_comm_inter_h.init();
    myid_inter = mpi_comm_inter_h.myid;
    size_inter = mpi_comm_inter_h.nprocs;

    librpa_mpi_initialized = true;
}

bool is_mpi_initialized()
{
    return librpa_mpi_initialized;
}

void finalize_mpi()
{
    procname = "not-init";
    mpi_comm_global_h.reset_comm();

    mpi_comm_inter_h.reset_comm();
    mpi_comm_intra_h.reset_comm();
    MPI_Comm_free(&mpi_comm_intra);
    MPI_Comm_free(&mpi_comm_inter);

    librpa_mpi_initialized = false;
}

} /* end of namespace envs */

} /* end of namespace LIBRPA */
