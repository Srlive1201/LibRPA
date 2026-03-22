#pragma once

#include "librpa_enums.h"
#include "librpa_options.h"
#include "librpa_handler.h"

#include <vector>
#include <complex>
#include <utility>

namespace librpa
{

/* Enums */
typedef LibrpaParallelRouting ParallelRouting;

typedef LibrpaTimeFreqGrid TimeFreqGrid;

typedef LibrpaSwitch Switch;

typedef LibrpaKind Kind;

typedef LibrpaVerbose Verbose;


/* Options */
/* Straighforward inheritance
   IMPORTANT: DO NOT add extra member variables here, which will break the inheritated data layout
   New control options should be put under the LibrpaOptions C structure. */
class Options : public ::LibrpaOptions
{
public:
    Options() { ::librpa_init_options(this); }
    void set_output_dir(const char *output_dir) { ::librpa_set_output_dir(this, output_dir); }
};


/* Global environment functions */
const char* get_build_info(void);

int get_major_version(void);

int get_minor_version(void);

int get_patch_version(void);

void init_global(Switch switch_redirect_stdout = LIBRPA_SWITCH_OFF, const char *redirect_path = "stdout",
                 Switch switch_process_output = LIBRPA_SWITCH_ON);

void finalize_global(void);

void test(void);

void print_profile(void);

/* C++ handler definition, input and compute methods declarations */

// Declaration wrapper macro, without and with options
#define LIBRPA_CPP_H_METHOD_DECL_WRAP(ret, func, ...) ret func(__VA_ARGS__)
#define LIBRPA_CPP_H_METHOD_DECL_WRAP_WOPTS(ret, func, ...) ret func(const librpa::Options &opts, __VA_ARGS__)

class Handler
{
private:
    LibrpaHandler *h_ = nullptr;
public:
    Handler(): h_(nullptr) {};
    Handler(int comm);
    const LibrpaHandler *get_c_handler() const { return h_; }
    void create(int comm);
    void free();
    ~Handler();

    /* Input (set) functions */

    LIBRPA_CPP_H_METHOD_DECL_WRAP(void, set_scf_dimension,
                                  int nspins, int nkpts, int nstates, int nbasis);
    LIBRPA_CPP_H_METHOD_DECL_WRAP(void, set_wg_ekb_efermi,
                                  int nspins, int nkpts, int nstates, const double* wg, const double* ekb, double efermi);
    LIBRPA_CPP_H_METHOD_DECL_WRAP(void, set_wfc,
                                  int ispin, int ik, int nstates_local, int nbasis_local,
                                  const double* wfc_real, const double* wfc_imag);
    LIBRPA_CPP_H_METHOD_DECL_WRAP(void, set_wfc_packed,
                                  int ispin, int ik, int nstates_local, int nbasis_local, const std::complex<double> *wfc);
    LIBRPA_CPP_H_METHOD_DECL_WRAP(void, set_ao_basis_wfc, const std::vector<size_t> &nbs_wfc);
    LIBRPA_CPP_H_METHOD_DECL_WRAP(void, set_ao_basis_aux, const std::vector<size_t> &nbs_aux);
    LIBRPA_CPP_H_METHOD_DECL_WRAP(void, set_latvec_and_G,
                                  const double lat_mat[9], const double G_mat[9]);
    LIBRPA_CPP_H_METHOD_DECL_WRAP(void, set_atoms,
                                  const std::vector<int> &types, const std::vector<double> &pos_cart);
    LIBRPA_CPP_H_METHOD_DECL_WRAP(void, set_kgrids_kvec,
                                  int nk1, int nk2, int nk3, const double* kvecs);
    LIBRPA_CPP_H_METHOD_DECL_WRAP(void, set_ibz_mapping,
                                  const std::vector<int> &map_ibzk);
    LIBRPA_CPP_H_METHOD_DECL_WRAP(void, set_lri_coeff,
                                  LibrpaParallelRouting routing,
                                  int I, int J, int nbasis_i, int nbasis_j, int naux_mu, const int R[3], const double* Cs_in);
    LIBRPA_CPP_H_METHOD_DECL_WRAP(void, set_aux_bare_coulomb_k_atom_pair,
                                  int ik, int I, int J, int naux_mu, int naux_nu,
                                  const double* Vq_real_in, const double* Vq_imag_in, double vq_threshold);
    LIBRPA_CPP_H_METHOD_DECL_WRAP(void, set_aux_cut_coulomb_k_atom_pair,
                                  int ik, int I, int J, int naux_mu, int naux_nu,
                                  const double* Vq_real_in, const double* Vq_imag_in, double vq_threshold);
    LIBRPA_CPP_H_METHOD_DECL_WRAP(void, set_aux_bare_coulomb_k_2d_block,
                                  int ik, int mu_begin, int mu_end, int nu_begin, int nu_end,
                                  const double* Vq_real_in, const double* Vq_imag_in);
    LIBRPA_CPP_H_METHOD_DECL_WRAP(void, set_aux_cut_coulomb_k_2d_block,
                                  int ik, int mu_begin, int mu_end, int nu_begin, int nu_end,
                                  const double* Vq_real_in, const double* Vq_imag_in);
    LIBRPA_CPP_H_METHOD_DECL_WRAP(void, set_dielect_func_imagfreq,
                                  const std::vector<double> &omegas_imag, const std::vector<double> &dielect_func);

    LIBRPA_CPP_H_METHOD_DECL_WRAP(void, set_band_kvec,
                                  int n_kpts_band, const double* kfrac_band);
    LIBRPA_CPP_H_METHOD_DECL_WRAP(void, set_band_occ_eigval,
                                  int n_spins, int n_kpts_band, int n_states,
                                  const double* occ, const double* eig);
    LIBRPA_CPP_H_METHOD_DECL_WRAP(void, set_wfc_band,
                                  int ispin, int ik_band, int nstates_local, int nbasis_local,
                                  const double* wfc_real, const double* wfc_imag);
    LIBRPA_CPP_H_METHOD_DECL_WRAP(void, set_wfc_band_packed,
                                  int ispin, int ik_band, int nstates_local, int nbasis_local, const std::complex<double> *wfc);
    LIBRPA_CPP_H_METHOD_DECL_WRAP(void, reset_band_data);

    /* Compute (build/get) functions */
    //! Construct and return frequency grids
    std::pair<std::vector<double>, std::vector<double>>
    get_imaginary_frequency_grids(const Options &opts);

    //! Compute RPA correlation energy
    LIBRPA_CPP_H_METHOD_DECL_WRAP_WOPTS(double, get_rpa_correlation_energy,
                                        std::vector<std::complex<double>> &rpa_corr_ibzk_contrib);
    //! Compute real-space exact-exchange matrix
    void build_exx(const Options &opts);

    std::vector<double>
    get_exx_pot_kgrid(const Options &opts, const int n_spins, const std::vector<int> &iks_this,
                      int i_state_low, int i_state_high);

    std::vector<double>
    get_exx_pot_band_k(const Options &opts, const int n_spins, const std::vector<int> &iks_band_this,
                       int i_state_low, int i_state_high);

    //! Compute real-space G0W0 self-energy matrix
    void build_g0w0_sigma(const Options &opts);

    std::vector<std::complex<double>>
    get_g0w0_sigc_kgrid(const Options &opts, const int n_spins, const std::vector<int> &iks_this,
                        int i_state_low, int i_state_high, const std::vector<double> &vxc, const std::vector<double> &vexx);

    std::vector<std::complex<double>>
    get_g0w0_sigc_band_k(const Options &opts, const int n_spins, const std::vector<int> &iks_band_this,
                         int i_state_low, int i_state_high, const std::vector<double> &vxc_band, const std::vector<double> &vexx_band);

    /* Utility functions */
};

#undef LIBRPA_CPP_H_METHOD_DECL_WRAP

}
