#pragma once
#include <array>
#include <complex>
#include <memory>
#include <valarray>
#include <mpi.h>
#include "params.h"

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \brief set dimension parameters of the system
 */
void set_dimension(int nspins, int nkpts, int nstates, int nbasis, int natoms);

void set_wg_ekb_efermi(int nspins, int nkpts, int nstates, double* wg, double* ekb, double efermi);

/*!
 * \brief set AO basis for wave function expansion
 */
void set_ao_basis_wfc(int is, int ik, double* wfc_real, double* wfc_imag);

void set_latvec_and_G(double* lat_mat, double* G_mat);
/*!
 * \brief set kmesh grids
 */
void set_kgrids_kvec_tot(int nk1, int nk2, int nk3, double* kvecs);

void set_ibz2bz_index_and_weight(const int nk_irk, const int* ibz2bz_index, const double* wk_irk);

/*!
 * \brief set auxiliary AO basis
 */
void set_ao_basis_aux(int I, int J, int nbasis_i, int nbasis_j, int naux_mu, int* R, double* Cs_in);

void set_aux_coulomb_k_atom_pair(int I, int J, int naux_mu, int naux_nu, int ik, double* Vq_real_in, double* Vq_imag_in);

void set_librpa_params();
void run_librpa_main(MPI_Comm comm_in);

#ifdef __cplusplus
}
#endif
