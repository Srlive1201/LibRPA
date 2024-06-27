#pragma once
#include <array>
#include <complex>
#include <memory>
#include <valarray>
#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * \brief C struct to handle input parameters
 *
 * The struct should have essentially the same members as the Params C++ struct in params.h.
 * However, not all members are implemented here, because some parameters are used to control
 * how to read the file, i.e. for the driver. They are not relevant in API calls.
 *
 * A better solution is to move these parameters to other place.
 */
struct LibRPAParams
{
    // char array types
    char task[20];
    char output_file[100];
    char output_dir[100];
    // char parallel_routing[20];
    // char tfgrids_type[20];

    // integer types
    int nfreq;

    // integer types, but correspondent to bool in Params
    // int debug;
    // int use_scalapack_ecrpa;

    // double types
    // double gf_R_threshold;
    // double cs_threshold;
    // double vq_threshold;
    // double sqrt_coulomb_threshold;
    // double libri_chi0_threshold_C;
    // double libri_chi0_threshold_G;
    // double libri_exx_threshold_CSM;
    // double libri_exx_threshold_C;
    // double libri_exx_threshold_D;
    // double libri_exx_threshold_V;
};

/*!
 * \brief set dimension parameters of the system
 */
void set_dimension(int nspins, int nkpts, int nstates, int nbasis, int natoms);

/*!
 * \brief initialize the environment of LibRPA calculation
 */
void initialize_librpa_environment(
        MPI_Comm comm_in, int is_fortran_comm,
        int redirect_stdout, const char *output_filename);

/*!
 * \brief finalie the environment of LibRPA calculation
 */
void finalize_librpa_environment();

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
void set_ao_basis_aux(int I, int J, int nbasis_i, int nbasis_j, int naux_mu, int* R, double* Cs_in, int insert_index_only);

void set_aux_coulomb_k_atom_pair(int I, int J, int naux_mu, int naux_nu, int ik, double* Vq_real_in, double* Vq_imag_in);
void set_aux_coulomb_k_2D_block(int ik, int max_naux, int mu_begin, int mu_end, int nu_begin, int nu_end, double* Vq_real_in, double* Vq_imag_in );
void set_librpa_params(LibRPAParams *params);
void get_default_librpa_params(LibRPAParams *params);
void run_librpa_main();

#ifdef __cplusplus
}
#endif
