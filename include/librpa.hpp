#pragma once
/**
 * @file librpa.hpp
 * @brief C++ wrapper API for LibRPA.
 *
 * This file provides a convenient C++ interface to LibRPA functionality.
 * It wraps the C API with STL containers for easier use.
 *
 * Usage:
 * @code
 * #include "librpa.hpp"
 * @endcode
 *
 * @note This header includes librpa_options.h and librpa_handler.h which provide the underlying C types.
 */

#include "librpa_enums.h"
#include "librpa_options.h"
#include "librpa_handler.h"

#include <vector>
#include <complex>
#include <utility>

/**
 * @brief Main namespace for LibRPA C++ API.
 */
namespace librpa
{

/* Enums */
/**
 * @brief Parallel routing strategy (C++ alias).
 * @see LibrpaParallelRouting
 */
typedef LibrpaParallelRouting ParallelRouting;

typedef LibrpaTimeFreqGrid TimeFreqGrid;  ///< Time/frequency grid type (C++ alias)

typedef LibrpaSwitch Switch;  ///< Boolean switch type (C++ alias)

typedef LibrpaKind Kind;  ///< DFT code kind type (C++ alias, reserved)

typedef LibrpaVerbose Verbose;  ///< Verbosity level type (C++ alias)


/* Options */
/**
 * @brief C++ wrapper for runtime options.
 *
 * Provides a convenient C++ interface to LibrpaOptions with constructor
 * initialization and setter methods.
 *
 * @note Straightforward inheritance. DO NOT add extra member variables here,
 *       which will break the inherited data layout. New control options
 *       should be put under the LibrpaOptions C structure.
 */
class Options : public ::LibrpaOptions
{
public:
    /** @brief Construct with default options. */
    Options() { ::librpa_init_options(this); }

    /**
     * @brief Set output directory.
     *
     * @param[in] output_dir   Path for output data files.
     */
    void set_output_dir(const char *output_dir) { ::librpa_set_output_dir(this, output_dir); }
};


/* Global environment functions */

/**
 * @brief Get build information string.
 * @return C-string containing build information.
 * @see librpa_get_build_info
 */
const char* get_build_info(void);

/**
 * @brief Get major version number.
 * @return Major version (X in X.Y.Z).
 * @see librpa_get_major_version
 */
int get_major_version(void);

/**
 * @brief Get minor version number.
 * @return Minor version (Y in X.Y.Z).
 * @see librpa_get_minor_version
 */
int get_minor_version(void);

/**
 * @brief Get patch version number.
 * @return Patch version (Z in X.Y.Z).
 * @see librpa_get_patch_version
 */
int get_patch_version(void);

/**
 * @brief Initialize the global LibRPA environment.
 *
 * Must be called after MPI_Init() and before any other LibRPA functions.
 *
 * @param switch_redirect_stdout If true, redirect stdout to a file.
 * @param redirect_path          Path for redirected output (default: "stdout").
 * @param switch_process_output  If true, enable per-process output (default: true).
 * @see librpa_init_global
 */
void init_global(Switch switch_redirect_stdout = LIBRPA_SWITCH_OFF, const char *redirect_path = "stdout",
                 Switch switch_process_output = LIBRPA_SWITCH_ON);

/**
 * @brief Finalize the global LibRPA environment.
 *
 * Must be called after all LibRPA operations are complete.
 * @see librpa_finalize_global
 */
void finalize_global(void);

/**
 * @brief Run internal self-tests.
 * @see librpa_test
 */
void test(void);

/**
 * @brief Print profiling information.
 * @see librpa_print_profile
 */
void print_profile(void);

/* C++ handler definition, input and compute methods declarations */

/**
 * @brief Main handler class for LibRPA computations.
 *
 * This class wraps the C handler with RAII semantics for convenient
 * resource management. It provides methods to set input data and
 * perform RPA/EXX/G0W0 calculations.
 *
 * Typical usage:
 * @code{.cpp}
 * librpa::init_global();
 *
 * librpa::Handler h(MPI_COMM_WORLD);
 * h.set_scf_dimension(nspins, nkpts, nstates, nbasis);
 * // ... set more input data ...
 *
 * double ec_rpa = h.get_rpa_correlation_energy(opts, ibzk_contrib);
 * h.free();
 *
 * librpa::finalize_global();
 * @endcode
 */
class Handler
{
private:
    LibrpaHandler *h_ = nullptr;
public:
    /** @brief Default constructor (creates null handler). */
    Handler(): h_(nullptr) {};

    /**
     * @brief Construct and initialize handler with MPI communicator.
     * @param comm MPI communicator for parallel computation.
     */
    Handler(int comm);

    /** @brief Get underlying C handler (for internal use). */
    const LibrpaHandler *get_c_handler() const { return h_; }

    /**
     * @brief Initialize handler with given MPI communicator.
     * @param comm MPI communicator.
     */
    void init(int comm);

    /**
     * @brief Free handler and release resources.
     */
    void free();

    /**
     * @brief Defatul destructor - automatically frees handler if not already freed.
     *
     * @note Reply on free() to explicitly release resources.
     */
    ~Handler();

    /* Input (set) functions */

    /** @brief Set SCF wavefunction dimension (numbers of spins, k-points, states, basis). */
    void set_scf_dimension(int nspins, int nkpts, int nstates, int nbasis);

    /** @brief Set occupation numbers, eigenvalues, and Fermi level. */
    void set_wg_ekb_efermi(int nspins, int nkpts, int nstates, const double *wg, const double *ekb,
                           double efermi);

    /** @brief Set wavefunction coefficients (separated real/imag arrays). */
    void set_wfc(int ispin, int ik, int nstates_local, int nbasis_local, const double *wfc_real,
                 const double *wfc_imag);

    /** @brief Set wavefunction coefficients (packed complex array). */
    void set_wfc_packed(int ispin, int ik, int nstates_local, int nbasis_local,
                        const std::complex<double> *wfc);

    /** @brief Set atomic orbital basis size for wavefunctions. */
    void set_ao_basis_wfc(const std::vector<size_t> &nbs_wfc);

    /** @brief Set auxiliary atomic orbital basis size. */
    void set_ao_basis_aux(const std::vector<size_t> &nbs_aux);

    /** @brief Set lattice vectors and reciprocal lattice vectors. */
    void set_latvec_and_G(const double lat_mat[9], const double G_mat[9]);

    /** @brief Set atom types and Cartesian coordinates. */
    void set_atoms(const std::vector<int> &types, const std::vector<double> &pos_cart);

    /** @brief Set k-point grid vectors. */
    void set_kgrids_kvec(int nk1, int nk2, int nk3, const double *kvecs);

    /** @brief Set mapping from full k-point grid to irreducibleBZ. */
    void set_ibz_mapping(const std::vector<int> &map_ibzk);

    /** @brief Set local RI coefficients. */
    void set_lri_coeff(LibrpaParallelRouting routing, int I, int J, int nbasis_i, int nbasis_j,
                       int naux_mu, const int R[3], const double *Cs_in);

    /** @brief Set bare Coulomb matrix elements (atom-pair format). */
    void set_aux_bare_coulomb_k_atom_pair(int ik, int I, int J, int naux_mu, int naux_nu,
                                          const double *Vq_real_in, const double *Vq_imag_in,
                                          double vq_threshold);

    /** @brief Set truncated Coulomb matrix elements (atom-pair format). */
    void set_aux_cut_coulomb_k_atom_pair(int ik, int I, int J, int naux_mu, int naux_nu,
                                         const double *Vq_real_in, const double *Vq_imag_in,
                                         double vq_threshold);

    /** @brief Set bare Coulomb matrix elements (2D block format). */
    void set_aux_bare_coulomb_k_2d_block(int ik, int mu_begin, int mu_end, int nu_begin, int nu_end,
                                         const double *Vq_real_in, const double *Vq_imag_in);

    /** @brief Set truncated Coulomb matrix elements (2D block format). */
    void set_aux_cut_coulomb_k_2d_block(int ik, int mu_begin, int mu_end, int nu_begin, int nu_end,
                                        const double *Vq_real_in, const double *Vq_imag_in);

    /** @brief Set dielectric function on imaginary frequency axis. */
    void set_dielect_func_imagfreq(const std::vector<double> &omegas_imag,
                                   const std::vector<double> &dielect_func);

    /** @brief Set k-points for band structure calculations. */
    void set_band_kvec(int n_kpts_band, const double *kfrac_band);

    /** @brief Set occupation numbers and eigenvalues for band k-points. */
    void set_band_occ_eigval(int n_spins, int n_kpts_band, int n_states, const double *occ,
                             const double *eig);

    /** @brief Set wavefunction for band k-point (separated real/imag). */
    void set_wfc_band(int ispin, int ik_band, int nstates_local, int nbasis_local,
                      const double *wfc_real, const double *wfc_imag);

    /** @brief Set wavefunction for band k-point (packed complex). */
    void set_wfc_band_packed(int ispin, int ik_band, int nstates_local, int nbasis_local,
                             const std::complex<double> *wfc);

    /** @brief Reset band structure data. */
    void reset_band_data();

    /* Compute (build/get) functions */

    /** @brief Construct and return frequency grids. */
    void get_imaginary_frequency_grids(const Options &opts,
                                       std::vector<double> &omegas, std::vector<double> &weights);

    /** @brief Compute RPA correlation energy.
     *
     * @param[in] opts Runtime options.
     * @param[out] rpa_corr_ibzk_contrib Complex RPA correlation contribution per k-point.
     * @return Total RPA correlation energy.
     */
    double get_rpa_correlation_energy(const Options &opts,
                                      std::vector<std::complex<double>> &rpa_corr_ibzk_contrib);

    /** @brief Build exact-exchange matrix in real space.
     *
     * @param[in] opts  Runtime options.
     * */
    void build_exx(const Options &opts);

    /** @brief Get exact-exchange potential for k-grid states.
     *
     * @param[in] opts         Runtime options.
     * @param[in] n_spins      Number of spin channels.
     * @param[in] iks_this     List of k-point indices computed on this process.
     * @param[in] i_state_low  First state index (inclusive).
     * @param[in] i_state_high Last state index (exclusive).
     * @return Exact-exchange potentials for selected states.
     */
    std::vector<double>
    get_exx_pot_kgrid(const Options &opts, const int n_spins, const std::vector<int> &iks_this,
                      int i_state_low, int i_state_high);

    /** @brief Get exact-exchange potential for band k-points.
     *
     * @param[in] opts           Runtime options.
     * @param[in] n_spins        Number of spin channels.
     * @param[in] iks_band_this  List of band k-point indices on this process.
     * @param[in] i_state_low    First state index (inclusive).
     * @param[in] i_state_high   Last state index (exclusive).
     *
     * @return Exact-exchange potentials for band states.
     */
    std::vector<double>
    get_exx_pot_band_k(const Options &opts, const int n_spins, const std::vector<int> &iks_band_this,
                       int i_state_low, int i_state_high);

    /** @brief Build G0W0 self-energy matrix in real space.
     *
     * @param[in] opts  Runtime options.
     * */
    void build_g0w0_sigma(const Options &opts);

    /** @brief Get G0W0 correlation self-energy for k-grid states.
     *
     * @param[in] opts         Runtime options.
     * @param[in] n_spins      Number of spin channels.
     * @param[in] iks_this     List of k-point indices on this process.
     * @param[in] i_state_low  First state index (inclusive).
     * @param[in] i_state_high Last state index (exclusive).
     * @param[in] vxc          XC potential for selected states.
     * @param[in] vexx         Exact-exchange potential for selected states.
     *
     * @return Correlation self-energy for selected states.
     */
    std::vector<std::complex<double>>
    get_g0w0_sigc_kgrid(const Options &opts, const int n_spins, const std::vector<int> &iks_this,
                        int i_state_low, int i_state_high, const std::vector<double> &vxc, const std::vector<double> &vexx);

    /** @brief Get G0W0 correlation self-energy for band k-points.
     *
     * @param[in] opts            Runtime options.
     * @param[in] n_spins         Number of spin channels.
     * @param[in] iks_band_this   List of band k-point indices on this process.
     * @param[in] i_state_low     First state index (inclusive).
     * @param[in] i_state_high    Last state index (exclusive).
     * @param[in] vxc_band        XC potential for band states.
     * @param[in] vexx_band       Exact-exchange potential for band states.
     *
     * @return Correlation self-energy for band states.
     */
    std::vector<std::complex<double>>
    get_g0w0_sigc_band_k(const Options &opts, const int n_spins, const std::vector<int> &iks_band_this,
                         int i_state_low, int i_state_high, const std::vector<double> &vxc_band, const std::vector<double> &vexx_band);

    /* Utility functions */
};

}
