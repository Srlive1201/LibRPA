#pragma once
#include <string>
#include <mpi.h>

namespace LIBRPA
{

namespace envs
{

//! Global communicator
extern MPI_Comm mpi_comm_global;

//! Rank of process in the global communciator
extern int myid_global;

//! Number of processes in the global communciator
extern int size_global;

//! Name of process
extern std::string procname;

//! Initialize the MPI environment of LibRPA
/*!
 * @param  [in]  mpi_comm_global_in    Global MPI communicator
 */
void initialize_mpi(MPI_Comm mpi_comm_global_in);

//! Check the MPI environment of LibRPA is correctly initialized
bool is_mpi_initialized();

//! Finalize the MPI environment of LibRPA
void finalize_mpi();

} /* end of namespace envs */

} /* end of namespace LIBRPA */
