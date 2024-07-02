#include "app_rpa.h"

#include "utils_mem.h"

#include "chi0.h"
#include "epsilon.h"
#include "envs_mpi.h"
#include "parallel_mpi.h"
#include "params.h"
#include "pbc.h"
#include "profiler.h"
#include "ri.h"

#include "stl_io_helper.h"

namespace LIBRPA
{

namespace app
{

void get_rpa_correlation_energy_(std::complex<double> &rpa_corr,
                                 std::vector<std::complex<double>> &rpa_corr_irk_contrib)
{
    using LIBRPA::envs::mpi_comm_global_h;
    using LIBRPA::utils::lib_printf;

    Vector3_Order<int> period {kv_nmp[0], kv_nmp[1], kv_nmp[2]};
    auto Rlist = construct_R_grid(period);
    const int Rt_num = Rlist.size() * Params::nfreq;

    tot_atpair = generate_atom_pair_from_nat(natom, false);

    set_parallel_routing(Params::parallel_routing, tot_atpair.size(), Rt_num, LIBRPA::parallel_routing);

    Chi0 chi0(meanfield, klist, Params::nfreq);
    chi0.gf_R_threshold = Params::gf_R_threshold;
    vector<Vector3_Order<double>> qlist;
    for ( auto q_weight: irk_weight)
    {
        qlist.push_back(q_weight.first);
    }

    mpi_comm_global_h.barrier();
    if (mpi_comm_global_h.myid == 0)
    {
        cout << "Initialization finished, start task job from myid\n";
    }

    Profiler::start("chi0_build", "Build response function chi0");
    chi0.build(Cs, Rlist, period, local_atpair, qlist, TFGrids::get_grid_type(Params::tfgrids_type), true);
    Profiler::stop("chi0_build");

    if (Params::debug)
    { // debug, check chi0
        char fn[80];
        for (const auto &chi0q: chi0.get_chi0_q())
        {
            const int ifreq = chi0.tfg.get_freq_index(chi0q.first);
            for (const auto &q_IJchi0: chi0q.second)
            {
                const int iq = std::distance(klist.begin(), std::find(klist.begin(), klist.end(), q_IJchi0.first));
                for (const auto &I_Jchi0: q_IJchi0.second)
                {
                    const auto &I = I_Jchi0.first;
                    for (const auto &J_chi0: I_Jchi0.second)
                    {
                        const auto &J = J_chi0.first;
                        sprintf(fn, "chi0fq_ifreq_%d_iq_%d_I_%zu_J_%zu_id_%d.mtx", ifreq, iq, I, J, mpi_comm_global_h.myid);
                        print_complex_matrix_mm(J_chi0.second, Params::output_dir + "/" + fn, 1e-15);
                    }
                }
            }
        }
    }

    // NOTE: Cs is cleaned up.
    // This means that the behavior will be undefined if this function is called again
    Cs.clear();
    LIBRPA::utils::release_free_mem();

    mpi_comm_global_h.barrier();
    Profiler::start("EcRPA", "Compute RPA correlation Energy");
    CorrEnergy corr;
    if (Params::use_scalapack_ecrpa && (LIBRPA::parallel_routing == LIBRPA::ParallelRouting::ATOM_PAIR || LIBRPA::parallel_routing == LIBRPA::ParallelRouting::LIBRI))
    {
        if(meanfield.get_n_kpoints() == 1)
            corr = compute_RPA_correlation_blacs_2d_gamma_only(chi0, Vq);
        else
            corr = compute_RPA_correlation_blacs_2d(chi0, Vq);
    }
    else
        corr = compute_RPA_correlation(chi0, Vq);

    rpa_corr = corr.value;

    cout << qlist << "\n";
    for (const auto &irk_corr: corr.qcontrib)
    {
        auto ite = std::find(qlist.cbegin(), qlist.cend(), irk_corr.first);
        // cout << irk_corr.first << " " << irk_corr.second << "\n";
        const auto iq = std::distance(qlist.cbegin(), ite);
        // cout << iq << " " << irk_corr.second << "\n";
        rpa_corr_irk_contrib[iq] = irk_corr.second;
    }

    Profiler::stop("EcRPA");
}

} /* end of namespace app */

} /* end of namespace LIBRPA */
