// Public API headers
#include "librpa.hpp"
#include "librpa_global.h"
#include "librpa_input.h"
#include "librpa_compute.h"

#include <cstddef>

// Internal headers
#include "../utils/error.h"

namespace librpa
{

/* Global environment functions */
const char* get_build_info(void) { return ::librpa_get_build_info(); }

int get_major_version(void) { return ::librpa_get_major_version(); }

int get_minor_version(void) { return ::librpa_get_minor_version(); }

int get_micro_version(void) { return ::librpa_get_micro_version(); }

void init_global(LibrpaSwitch switch_redirect_stdout, const char *redirect_path,
                 LibrpaSwitch switch_process_output)
{
    ::librpa_init_global(switch_redirect_stdout, redirect_path, switch_process_output);
}

void finalize_global(void) { ::librpa_finalize_global(); }


/* Constructor and destructor definition for the C++ handler object */
Handler::Handler(int comm) : h_(nullptr) { this->create(comm); }

void Handler::create(int comm)
{
    if (h_) throw LIBRPA_RUNTIME_ERROR("double creating handler");  // safe guard
    h_ = ::librpa_create_handler(comm);
}

void Handler::free()
{
    if (h_) ::librpa_destroy_handler(h_);
    h_ = nullptr;
}

Handler::~Handler() { free(); }

/* Input methods definition for the C++ handler object */

// The following macros are helper to implement the methods of the C++ handler class
// that wraps the corresponding `librpa_` prefixed C-API functions.
#define LIBRPA_C_NAME(name) librpa_##name
#define LIBRPA_UNPAREN(x) LIBRPA_UNPAREN_ x
#define LIBRPA_UNPAREN_(...) __VA_ARGS__

#define LIBRPA_CPP_H_METHOD_DEF_WRAP_VOID(cpp_name, params, args) \
    void Handler::cpp_name(LIBRPA_UNPAREN(params))                       \
    {                                                                    \
        ::LIBRPA_C_NAME(cpp_name)(this->h_, LIBRPA_UNPAREN(args));       \
    }

#define LIBRPA_CPP_H_METHOD_DEF_WRAP_VOID_WOPT(cpp_name, params, args)        \
    void Handler::cpp_name(const LibrpaOptions* opts, LIBRPA_UNPAREN(params)) \
    {                                                                         \
        ::LIBRPA_C_NAME(cpp_name)(this->h_, opts, LIBRPA_UNPAREN(args));      \
    }

#define LIBRPA_CPP_H_METHOD_DEF_WRAP_RET_NOOPT(ret, cpp_name, params, args) \
    ret Handler::cpp_name(LIBRPA_UNPAREN(params))                           \
    {                                                                       \
        return ::LIBRPA_C_NAME(cpp_name)(this->h_, LIBRPA_UNPAREN(args));   \
    }

#define LIBRPA_CPP_H_METHOD_DEF_WRAP_RET(ret, cpp_name, params, args)           \
    ret Handler::cpp_name(const LibrpaOptions* opts, LIBRPA_UNPAREN(params))    \
    {                                                                           \
        return ::LIBRPA_C_NAME(cpp_name)(this->h_, opts, LIBRPA_UNPAREN(args)); \
    }

// LIBRPA_CPP_H_METHOD_DEF_WRAP_VOID(
//     set_dimension,
//     (int nspins, int nkpts, int nstates, int nbasis),
//     (nspins, nkpts, nstates, nbasis)
// )
// // Equivalent, write it explicitly to tell the language server that librpa_input.h is needed
void Handler::set_scf_dimension(int nspins, int nkpts, int nstates, int nbasis)
{
    ::librpa_set_scf_dimension(this->h_, nspins, nkpts, nstates, nbasis);
}

LIBRPA_CPP_H_METHOD_DEF_WRAP_VOID(
    set_wg_ekb_efermi,
    (int nspins, int nkpts, int nstates, const double* wg, const double* ekb, double efermi),
    (nspins, nkpts, nstates, wg, ekb, efermi)
)

LIBRPA_CPP_H_METHOD_DEF_WRAP_VOID(
    set_wfc,
    (int ispin, int ik, const double* wfc_real, const double* wfc_imag),
    (ispin, ik, wfc_real, wfc_imag)
)

LIBRPA_CPP_H_METHOD_DEF_WRAP_VOID(
    set_ao_basis_wfc,
    (const std::vector<size_t> &nbs_wfc),
    (nbs_wfc.size(), nbs_wfc.data())
)

LIBRPA_CPP_H_METHOD_DEF_WRAP_VOID(
    set_ao_basis_aux,
    (const std::vector<size_t> &nbs_aux),
    (nbs_aux.size(), nbs_aux.data())
)

LIBRPA_CPP_H_METHOD_DEF_WRAP_VOID(
    set_latvec_and_G,
    (const double lat_mat[9], const double G_mat[9]),
    (lat_mat, G_mat)
)

LIBRPA_CPP_H_METHOD_DEF_WRAP_VOID(
    set_atoms,
    (const std::vector<int> &types, const std::vector<double> &pos_cart),
    (types.size(), types.data(), pos_cart.data())
)

LIBRPA_CPP_H_METHOD_DEF_WRAP_VOID(
    set_kgrids_kvec,
    (int nk1, int nk2, int nk3, const double* kvecs),
    (nk1, nk2, nk3, kvecs)
)

LIBRPA_CPP_H_METHOD_DEF_WRAP_VOID(
    set_ibz_mapping,
    (const std::vector<int> &map_ibzk),
    (map_ibzk.data())
)

LIBRPA_CPP_H_METHOD_DEF_WRAP_VOID(
    set_lri_coeff,
    (LibrpaParallelRouting routing, int I, int J, int nbasis_i, int nbasis_j, int naux_mu, const int R[3], const double* Cs_in),
    (routing, I, J, nbasis_i, nbasis_j, naux_mu, R, Cs_in)
)

LIBRPA_CPP_H_METHOD_DEF_WRAP_VOID(
    set_aux_bare_coulomb_k_atom_pair,
    (int ik, int I, int J, int naux_mu, int naux_nu, const double* Vq_real_in, const double* Vq_imag_in, double vq_threshold),
    (ik, I, J, naux_mu, naux_nu, Vq_real_in, Vq_imag_in, vq_threshold)
)

LIBRPA_CPP_H_METHOD_DEF_WRAP_VOID(
    set_aux_cut_coulomb_k_atom_pair,
    (int ik, int I, int J, int naux_mu, int naux_nu, const double* Vq_real_in, const double* Vq_imag_in, double vq_threshold),
    (ik, I, J, naux_mu, naux_nu, Vq_real_in, Vq_imag_in, vq_threshold)
)

LIBRPA_CPP_H_METHOD_DEF_WRAP_VOID(
    set_aux_bare_coulomb_k_2d_block,
    (int ik, int max_naux, int mu_begin, int mu_end, int nu_begin, int nu_end, const double* Vq_real_in, const double* Vq_imag_in),
    (ik, max_naux, mu_begin, mu_end, nu_begin, nu_end, Vq_real_in, Vq_imag_in)
)

LIBRPA_CPP_H_METHOD_DEF_WRAP_VOID(
    set_aux_cut_coulomb_k_2d_block,
    (int ik, int max_naux, int mu_begin, int mu_end, int nu_begin, int nu_end, const double* Vq_real_in, const double* Vq_imag_in),
    (ik, max_naux, mu_begin, mu_end, nu_begin, nu_end, Vq_real_in, Vq_imag_in)
)

}
