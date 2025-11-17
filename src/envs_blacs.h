#pragma once

#include "base_blacs.h"

namespace LIBRPA
{

namespace envs
{

//! Handler of the global BLACS context
extern LIBRPA::BLACS_CTXT_handler blacs_ctxt_global_h;

//! Array descriptor for matrix within wave function basis that shared by all MPI tasks
extern LIBRPA::Array_Desc array_desc_wfc_global;
//! Array descriptor for matrix within auxiliary basis that shared by all MPI tasks
extern LIBRPA::Array_Desc array_desc_abf_global;

//! Initialize the MPI environment of LibRPA
/*!
 * @param  [in]  mpi_comm_global_in    Global MPI communicator
 */
void initialize_blacs(const MPI_Comm &mpi_comm_global_in);

//! Check if the blacs environment is initialized
bool is_blacs_initialized();

//! Initialize global array descriptor for matrix of wave function basis size
/*!
 * The maximal block size will be used to ensure that the row/column indices
 * are continous on each process.
 *
 * @param  [in]  n_wfc    size of wave function basis
 *
 * @warn The BLACS environment must be initialized first (by initialize_blacs)
 */
void initialize_array_desc_wfc_global(const size_t &n_wfc);

//! Initialize global array descriptor for matrix of wave function basis size
/*!
 * The maximal block size will be used to ensure that the row/column indices
 * are continous on each process.
 *
 * @param  [in]  n_abf    size of auxiliary basis
 *
 * @warn The BLACS environment must be initialized first (by initialize_blacs)
 */
void initialize_array_desc_abf_global(const size_t &n_abf);

void finalize_blacs();

} /* end of namespace envs */

} /* end of namespace LIBRPA */
