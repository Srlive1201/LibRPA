#include "global_mpi.h"

#include <stdexcept>

namespace librpa_int
{

namespace global
{

MPI_Comm mpi_comm_global;
MpiCommHandler mpi_comm_global_h;
int myid_global = 0;
int size_global = 1;

std::string procname = "not-init";

MPI_Comm mpi_comm_intra;
MpiCommHandler mpi_comm_intra_h;
int myid_intra = 0;
int size_intra = 1;

MPI_Comm mpi_comm_inter;
MpiCommHandler mpi_comm_inter_h;
int myid_inter = 0;
int size_inter = 1;

static bool mpi_initialized = false;

void init_global_mpi()
{
    int flag;
    MPI_Initialized(&flag);
    if (flag == 0)
    {
        throw std::runtime_error(
            std::string(__FILE__) + ":" + std::to_string(__LINE__) + ":" + std::string(__FUNCTION__) + ": "
            "MPI_Init or MPI_Init_thread must be called before " + std::string(__FUNCTION__));
    }

    mpi_comm_global = MPI_COMM_WORLD;
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

    mpi_initialized = true;
}

bool is_mpi_initialized()
{
    return mpi_initialized;
}

void finalize_global_mpi()
{
    procname = "not-init";
    mpi_comm_global_h.reset_comm();
    mpi_comm_inter_h.reset_comm();
    mpi_comm_intra_h.reset_comm();
    MPI_Comm_free(&mpi_comm_intra);
    MPI_Comm_free(&mpi_comm_inter);

    mpi_initialized = false;
}

}

}
