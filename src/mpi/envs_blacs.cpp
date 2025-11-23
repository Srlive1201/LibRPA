#include "envs_blacs.h"

namespace LIBRPA
{

namespace envs
{

static bool librpa_blacs_initialized = false;

LIBRPA::BLACS_CTXT_handler blacs_ctxt_global_h;
LIBRPA::Array_Desc array_desc_wfc_global;
LIBRPA::Array_Desc array_desc_abf_global;

void initialize_blacs(const MPI_Comm &mpi_comm_global_in)
{
    blacs_ctxt_global_h.reset_comm(mpi_comm_global_in);
    blacs_ctxt_global_h.init();

    /*
     * HACK: A close-to-square process grid is imposed.
     *       This might subject to performance issue when
     *       the number of processes is not dividable.
     */
    blacs_ctxt_global_h.set_square_grid();

    librpa_blacs_initialized = true;
}

static void _initialize_array_desc_global(Array_Desc &ad, const size_t &n)
{
    assert(is_blacs_initialized());
    ad.reset_handler(blacs_ctxt_global_h);
    ad.init_1b1p(n, n, 0, 0);
}

void initialize_array_desc_wfc_global(const size_t &n_wfc)
{
    _initialize_array_desc_global(array_desc_wfc_global, n_wfc);
}

void initialize_array_desc_abf_global(const size_t &n_abf)
{
    _initialize_array_desc_global(array_desc_abf_global, n_abf);
}

bool is_blacs_initialized()
{
    return librpa_blacs_initialized;
}

void finalize_blacs()
{
    blacs_ctxt_global_h.exit();
    librpa_blacs_initialized = false;
}

} /* end of namespace envs */

} /* end of namespace LIBRPA */
