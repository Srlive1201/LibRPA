#pragma once
#include <array>
#include <complex>
#include <memory>
#include <valarray>

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \brief set dimension parameters of the system
 */
void set_dimension(int nspins, int nkpts, int nstates, int nbasis);

/*!
 * \brief set kmesh grids
 */
void set_kgrids(int nk1, int nk2, int nk3);

/*!
 * \brief set individual kpoints
 */
void set_kpoints();

/*!
 * \brief set AO basis for wave function expansion
 */
void set_ao_basis_wfc();

/*!
 * \brief set auxiliary AO basis
 */
void set_ao_basis_aux();

#ifdef __cplusplus
}
#endif
