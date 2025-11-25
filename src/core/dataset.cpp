#include "dataset.h"

namespace librpa_int
{

Dataset::Dataset(MPI_Comm comm)
    : comm_h(comm),
      blacs_ctxt_h(comm),
      basis_wfc(),
      basis_abf(),
      atoms(),
      pbc(),
      mf(),
      tfg(),
      cs_data(),
      vq(),
      vq_cut(),
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
