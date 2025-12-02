#include "dataset.h"

namespace librpa_int
{

Dataset::Dataset(MPI_Comm comm)
    : comm_h(comm, true),
      blacs_ctxt_h(comm),
      desc_abf(),
      basis_wfc(),
      basis_aux(),
      atoms(),
      pbc(),
      mf(),
      tfg(),
      cs_data(),
      vq(),
      vq_cut(),
      vq_block_loc(),
      vq_cut_block_loc(),
      p_exx(nullptr),
      p_chi0(nullptr),
      p_g0w0(nullptr)
{
    blacs_ctxt_h.init();
    // TODO more flexible process grid initialization
    blacs_ctxt_h.set_square_grid();
}

void Dataset::free()
{
}

}
