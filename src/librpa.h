/*!
 * @file      librpa.h
 * @brief     C interface of LibRPA
 * @author    LibRPA developers
 * @date      2024-06-28
 */
#pragma once
#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

/*!
 * @brief C struct to handle input parameters for LibRPA calculation
 */
/*
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
    char parallel_routing[20];
    char tfgrids_type[20];

    char DFT_software[20];
    // integer types
    int nfreq;

    // integer types, but correspondent to bool in Params
    int debug;
    int use_scalapack_ecrpa;

    // double types
    double gf_R_threshold;
    double cs_threshold;
    double vq_threshold;
    double sqrt_coulomb_threshold;
    double libri_chi0_threshold_C;
    double libri_chi0_threshold_G;
    double libri_exx_threshold_CSM;
    double libri_exx_threshold_C;
    double libri_exx_threshold_D;
    double libri_exx_threshold_V;
    double libri_gw_threshold_C;
    double libri_gw_threshold_G;
    double libri_gw_threshold_W;
};

/*!
 * @brief Initialize the environment of LibRPA calculation
 *
 * @param[in] comm_global_in     Global MPI communicator
 * @param[in] is_fortran_comm    Flag to identify whether the input communicator is Fortran.
 * @param[in] redirect_stdout    Flag to control whether output LibRPA will be printed to a file
 * @param[in] output_filename    Name of file for redirected output. Only used when redirect_stdout is true
 *
 * @todo is_fortran_comm could be removed when Fortran binding is implemented.
 */
void initialize_librpa_environment(
        MPI_Comm comm_global_in, int is_fortran_comm,
        int redirect_stdout, const char *output_filename);

/*!
 * @brief Finalize the environment of LibRPA calculation
 */
void finalize_librpa_environment();

/*!
 * @brief Set dimension parameters of the system
 *
 * @param[in] nspins    Number of spin channels
 * @param[in] nkpts     Number of k-points, on which the Kohn-Sham eigenvectors are computed
 * @param[in] nstates   Number of states
 * @param[in] nbasis    Total number of AO basis
 * @param[in] natoms    Total number of atoms
 */
void set_dimension(int nspins, int nkpts, int nstates, int nbasis, int natoms);

/*!
 * @brief Set eigenvalues, occupation number and fermi level.
 *
 * @param[in] nspins    Number of spin channels
 * @param[in] nkpts     Number of k-points, on which the Kohn-Sham eigenvectors are computed
 * @param[in] nstates   Number of states
 * @param[in] wg        Unnormalized occupation number, [nspins][nkpts][nstates]
 * @param[in] ekb       Eigenvalues in Hartree unit, [nspins][nkpts][nstates]
 * @param[in] efermi    Fermi level in Hartree unit
 */
void set_wg_ekb_efermi(int nspins, int nkpts, int nstates, double* wg, double* ekb, double efermi);

/*!
 * @brief Set wave function expansion of AO basis for a particular spin channel and k-point
 * @param[in] ispin      index of spin channel
 * @param[in] ik         index of k-point
 * @param[in] wfc_real   real part of wave function, [nstates][nbasis]
 * @param[in] wfc_imag   imaginary part of wave function, [nstates][nbasis]
 */
void set_ao_basis_wfc(int ispin, int ik, double* wfc_real, double* wfc_imag);

/*!
 * @brief Set the real-space and reciprocal-space lattice vectors
 *
 * @param[in] lat_mat    pointer to array of real-space lattice vectors, in Bohr unit
 * @param[in] G_mat      pointer to array of reciprocal-space lattice vectors, in Bohr^{-1} unit
 */
void set_latvec_and_G(double lat_mat[9], double G_mat[9]);

/*!
 * @brief Set k-mesh grids
 *
 * @param[in] nk1      Number of k-grids along the 1st reciprocal lattice vector
 * @param[in] nk2      Number of k-grids along the 2nd reciprocal lattice vector
 * @param[in] nk3      Number of k-grids along the 3rd reciprocal lattice vector
 * @param[in] kvecs    Coordinates of k-mesh vectors, in inverse Bohr unit, [nk1*nk2*nk3][3]
 */
void set_kgrids_kvec_tot(int nk1, int nk2, int nk3, double* kvecs);

/*!
 * @brief Set the mapping of irreducible k-points to the whole k-points set and their weights.
 *
 * @param[in] nk_irk          Number of irreducible k-points
 * @param[in] ibz2bz_index    Mapping of irreducible k-points to the full set, [nk_irk]
 * @param[in] wk_irk          Weights of irreducible k-points, [nk_irk]
 */
void set_ibz2bz_index_and_weight(const int nk_irk, const int* ibz2bz_index, const double* wk_irk);

/*!
 * @brief Insert the atom index and set the local RI triple coefficients
 *
 * @param[in] I                    Index of atom, where one basis and the auxiliary basis reside.
 * @param[in] J                    Index of atom, where the other basis resides.
 * @param[in] nbasis_i             Number of basis functions centered at atom I
 * @param[in] nbasis_j             Number of basis functions centered at atom J
 * @param[in] naux_mu              Number of auxliary basis functions centered at atom I
 * @param[in] R                    Real-space lattice vector where atom J resides, [3]
 * @param[in] Cs_in                Local RI triple coefficients, [nbasis_i][nbasis_j][naux_mu]
 * @param[in] insert_index_only    Only insert the indices of atoms but skip setting Cs_in
 */
void set_ao_basis_aux(int I, int J, int nbasis_i, int nbasis_j, int naux_mu, int* R, double* Cs_in, int insert_index_only);

/*!
 * @brief Set the atom-pair block of bare Coulomb matrix in auxiliary basis
 *
 * @param[in] ik            Index of k-point of parsed Coulomb matrix
 * @param[in] I             Index of atom for basis functions of the row indices
 * @param[in] J             Index of atom for basis functions of the column indices
 * @param[in] naux_mu       Number of auxliary basis functions centered at atom I
 * @param[in] naux_nu       Number of auxliary basis functions centered at atom J
 * @param[in] Vq_real_in    Real part of the Coulomb matrix block
 * @param[in] Vq_imag_in    Imaginary part of the Coulomb matrix block
 */
void set_aux_bare_coulomb_k_atom_pair(int ik, int I, int J, int naux_mu, int naux_nu, double* Vq_real_in, double* Vq_imag_in);

/*!
 * @brief Set the atom-pair block of truncated (cut) Coulomb matrix in auxiliary basis
 *
 * @param[in] ik            Index of k-point of parsed Coulomb matrix
 * @param[in] I             Index of atom for basis functions of the row indices
 * @param[in] J             Index of atom for basis functions of the column indices
 * @param[in] naux_mu       Number of auxliary basis functions centered at atom I
 * @param[in] naux_nu       Number of auxliary basis functions centered at atom J
 * @param[in] Vq_real_in    Real part of the truncated Coulomb matrix block
 * @param[in] Vq_imag_in    Imaginary part of the truncated Coulomb matrix block
 */
void set_aux_cut_coulomb_k_atom_pair(int ik, int I, int J, int naux_mu, int naux_nu, double* Vq_real_in, double* Vq_imag_in);

/*!
 * @brief Set the bare Coulomb matrix in auxiliary basis under BLACS 2D block format
 *
 * @param[in] ik            Index of k-point of parsed Coulomb matrix
 * @param[in] max_naux      Total number of auxiliary basis functions
 * @param[in] mu_begin      Starting row index
 * @param[in] mu_end        End row index
 * @param[in] nu_begin      Starting column index
 * @param[in] nu_end        End column index
 * @param[in] Vq_real_in    Real part of the Coulomb matrix block
 * @param[in] Vq_imag_in    Imaginary part of the Coulomb matrix block
 */
void set_aux_bare_coulomb_k_2D_block(int ik, int max_naux, int mu_begin, int mu_end, int nu_begin, int nu_end, double* Vq_real_in, double* Vq_imag_in);

/*!
 * @brief Set the truncated (cut) Coulomb matrix in auxiliary basis under BLACS 2D block format
 *
 * @param[in] ik            Index of k-point of parsed Coulomb matrix
 * @param[in] max_naux      Total number of auxiliary basis functions
 * @param[in] mu_begin      Starting row index
 * @param[in] mu_end        End row index
 * @param[in] nu_begin      Starting column index
 * @param[in] nu_end        End column index
 * @param[in] Vq_real_in    Real part of the Coulomb matrix block
 * @param[in] Vq_imag_in    Imaginary part of the Coulomb matrix block
 */
void set_aux_cut_coulomb_k_2D_block(int ik, int max_naux, int mu_begin, int mu_end, int nu_begin, int nu_end, double* Vq_real_in, double* Vq_imag_in);

/*!
 * @brief Set the control parameters of LibRPA
 *
 * @param[in] params    Struct containing the control parameters
 */
void set_librpa_params(LibRPAParams *params);

/*!
 * @brief Obtain the default controlling parameters of LibRPA
 *
 * @param[in] params    Struct containing the control parameters
 */
void get_default_librpa_params(LibRPAParams *params);

/*!
 * @brief Obtain the frequency grid points
 *
 * @param[in] params    Struct containing the control parameters
 */
void get_frequency_grids(int ngrid, double *freqeuncy_grids);

void run_librpa_main();

/*!
 * @brief compute and return the RPA correlation energy
 *
 * @param[out] rpa_corr                RPA correlation energy, in Hartree unit. Complex number represented as double array [2].
 * @param[out] rpa_corr_irk_contrib    Weighted contribution to RPA correlation energy from each irreducible k-point.
 *                                     Complex array represented as double array [2*n_irk_points]
 *
 * @warning
 * Current implementation will free the RI coefficients array. So this can only be called once.
 */
void get_rpa_correlation_energy(double *rpa_corr, double *rpa_corr_irk_contrib);

/*!
 * @brief Compute the exact exchange (EXX) energy for states at specified k-points
 *
 * @param[in]  i_state_low       The lowest index of state (included) to compute
 * @param[in]  i_state_high      The highest index of state (excluded) to compute
 * @param[in]  n_kpoints_task    The number of k-points to return in the called process.
 *                               When equal to 0, an empty vector will be returned.
 *                               When less than 0, all k-points will be computed.
 *                               Otherwise, the states at k-points whose indices
 *                               are stored in `i_kpoints_task` will be computed.
 * @param[in]  i_kpoints_task    The indices of k-points to compute EXX energy.
 * @param[out] exx               Exchange energy. It should have length of at least
 *                               `n_spins` * `n_kpoints_task` * (`i_state_high` - `i_state_low`).
 *                               When `n_kpoints_task` < 0, the length should be at least
 *                               `n_spins` * `n_kpoints` * (`i_state_high` - `i_state_low`).
 */
void compute_exx_orbital_energy(int i_state_low, int i_state_high,
                                int n_kpoints_task, const int *i_kpoints_task,
                                double *exx);

// /*
//  * @brief compute the screened Coulomb matrix in auxiliary basis representation
//  */
// void compute_screened_coulomb_abf(double _Complex *screened_coul);

// void compute_gw_quasiparticle_energy_kgrid(int n_qp_state_low, int n_qp_state_high, const double *vxc, double _Complex *screened_coul);

#ifdef __cplusplus
}
#endif
