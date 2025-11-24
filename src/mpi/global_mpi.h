#pragma once
#include <string>

#include "../mpi/base_mpi.h"

namespace librpa_int
{

namespace global
{

//! Global communicator
extern MPI_Comm mpi_comm_global;
//! Handler of the global communicator
extern MpiCommHandler mpi_comm_global_h;
//! Rank of process in the global communciator
extern int myid_global;
//! Number of processes in the global communciator
extern int size_global;

//! Name of process
extern std::string procname;

//! Intra-node communicator
extern MPI_Comm mpi_comm_intra;
//! Handler of the intra-node communicator
extern MpiCommHandler mpi_comm_intra_h;
//! Rank of process in the intra-node communciator
extern int myid_intra;
//! Number of processes in the intra-node communciator
extern int size_intra;

//! Inter-node communicator
extern MPI_Comm mpi_comm_inter;
//! Handler of the inter-node communicator
extern MpiCommHandler mpi_comm_inter_h;
//! Rank of process in the intra-node communciator
extern int myid_inter;
//! Number of processes in the intra-node communciator
extern int size_inter;

void init_global_mpi();

bool is_mpi_initialized();

void finalize_global_mpi();

}

}
