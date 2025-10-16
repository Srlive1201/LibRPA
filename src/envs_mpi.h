/*!
 * @file      envs_mpi.h
 * @brief     Environment setup regarding MPI.
 * @author    Min-Ye Zhang
 * @date      2024-07-02
 */
#pragma once
#include <string>
#include <mpi.h>

#include "base_mpi.h"

//! @namespace LIBRPA           Global namespace for LibRPA functions and classes
//! @namespace LIBRPA::envs     Facilities regarding computation envrionment
//! @namespace LIBRPA::utils    Utilities functions
//! @namespace LIBRPA::app      Functions for applications, e.g. RPA correlation, exact-exchange matrix

namespace LIBRPA
{

namespace envs
{

//! Global communicator
extern MPI_Comm mpi_comm_global;
//! Handler of the global communicator
extern LIBRPA::MPI_COMM_handler mpi_comm_global_h;
//! Rank of process in the global communciator
extern int myid_global;
//! Number of processes in the global communciator
extern int size_global;

//! Name of process
extern std::string procname;

//! Intra-node communicator
extern MPI_Comm mpi_comm_intra;
//! Handler of the intra-node communicator
extern LIBRPA::MPI_COMM_handler mpi_comm_intra_h;
//! Rank of process in the intra-node communciator
extern int myid_intra;
//! Number of processes in the intra-node communciator
extern int size_intra;

//! Inter-node communicator
extern MPI_Comm mpi_comm_inter;
//! Handler of the inter-node communicator
extern LIBRPA::MPI_COMM_handler mpi_comm_inter_h;
//! Rank of process in the intra-node communciator
extern int myid_inter;
//! Number of processes in the intra-node communciator
extern int size_inter;

//! Initialize the MPI environment of LibRPA
/*!
 * @param  [in]  mpi_comm_global_in    Global MPI communicator
 */
void initialize_mpi(const MPI_Comm &mpi_comm_global_in);

//! Check whether the MPI environment of LibRPA is correctly initialized
bool is_mpi_initialized();

//! Finalize the MPI environment of LibRPA
void finalize_mpi();

} /* end of namespace envs */

} /* end of namespace LIBRPA */
