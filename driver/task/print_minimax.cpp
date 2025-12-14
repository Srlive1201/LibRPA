#include "../task.h"
#include "../driver.h"
#include "../../src/api/instance_manager.h"

void driver::task_print_minimax()
{
    using namespace librpa_int;
    double emin, emax;

    auto pds = librpa_int::api::get_dataset_instance(h);
    auto &mf = pds->mf;
    mf.get_E_min_max(emin, emax);

    auto &tfg = pds->tfg;
    tfg.reset(driver::opts.nfreq);
    tfg.generate_minimax(emin, emax);
    if (global::mpi_comm_global_h.is_root()) tfg.show();
}
