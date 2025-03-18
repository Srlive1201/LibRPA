#include "app_rpa.h"

#include "chi0.h"
#include "envs_mpi.h"
#include "epsilon.h"
#include "parallel_mpi.h"
#include "params.h"
#include "pbc.h"
#include "profiler.h"
#include "ri.h"
#include "stl_io_helper.h"
#include "utils_mem.h"
#include "utils_timefreq.h"

namespace LIBRPA
{

namespace app
{

void get_rpa_correlation_energy_(std::complex<double> &rpa_corr,
                                 std::vector<std::complex<double>> &rpa_corr_irk_contrib,
                                 std::map<Vector3_Order<double>, ComplexMatrix> &sinvS,
                                 const std::string &input_dir, const bool use_shrink_abfs)
{
    using LIBRPA::envs::mpi_comm_global_h;
    using LIBRPA::utils::lib_printf;

    Vector3_Order<int> period{kv_nmp[0], kv_nmp[1], kv_nmp[2]};
    auto Rlist = construct_R_grid(period);
    const int Rt_num = Rlist.size() * Params::nfreq;

    tot_atpair = generate_atom_pair_from_nat(natom, false);

    set_parallel_routing(Params::parallel_routing, tot_atpair.size(), Rt_num,
                         LIBRPA::parallel_routing);

    // Build time-frequency objects
    auto tfg = utils::generate_timefreq_grids(Params::nfreq, Params::tfgrids_type, meanfield);

    Chi0 chi0(meanfield, klist, tfg);
    chi0.gf_R_threshold = Params::gf_R_threshold;
    vector<Vector3_Order<double>> qlist;
    for (auto q_weight : irk_weight)
    {
        qlist.push_back(q_weight.first);
    }

    mpi_comm_global_h.barrier();

    auto atom_mu_s = atom_mu;
    if (use_shrink_abfs)
    {
        // replace atom_mu by atom_mu_l to construct chi0 due to LRI error
        atom_mu = atom_mu_l;
        LIBRPA::atomic_basis_abf.set(atom_mu);
        atom_mu_part_range.resize(atom_mu.size());
        atom_mu_part_range[0] = 0;
        for (int I = 1; I != atom_mu.size(); I++)
            atom_mu_part_range[I] = atom_mu.at(I - 1) + atom_mu_part_range[I - 1];

        N_all_mu = atom_mu_part_range[natom - 1] + atom_mu[natom - 1];
    }

    Profiler::start("chi0_build", "Build response function chi0");
    chi0.build(Cs_data, Rlist, period, local_atpair, qlist);
    Profiler::stop("chi0_build");
    auto &chi0_q = chi0.get_chi0_q();
    print_complex_matrix_mm(chi0_q.at(tfg.get_freq_nodes()[0]).at(qlist[0]).at(0).at(0), "chi0",
                            0.);

    if (use_shrink_abfs)
    {
        print_complex_matrix_mm(sinvS.at(qlist[0]), "sinvS", 0.);
        auto identity = sinvS.at(qlist[0]) * transpose(sinvS.at(qlist[0]), true);
        print_complex_matrix_mm(identity, "sinvS_Ssinv", 0.);
        //  change atom_mu: number of {Mu,mu} in the later calculations
        atom_mu = atom_mu_s;
        LIBRPA::atomic_basis_abf.set(atom_mu);
        atom_mu_part_range.resize(atom_mu.size());
        atom_mu_part_range[0] = 0;
        for (int I = 1; I != atom_mu.size(); I++)
            atom_mu_part_range[I] = atom_mu.at(I - 1) + atom_mu_part_range[I - 1];

        N_all_mu = atom_mu_part_range[natom - 1] + atom_mu[natom - 1];
        Profiler::start("shrink_chi0_abfs", "Do shrink transformation");
        chi0.shrink_abfs_chi0(sinvS, qlist, atom_mu_l);
        Profiler::stop("shrink_chi0_abfs");
        sinvS.clear();
    }

    if (Params::debug)
    {  // debug, check chi0
        char fn[80];
        for (const auto &chi0q : chi0.get_chi0_q())
        {
            const int ifreq = chi0.tfg.get_freq_index(chi0q.first);
            for (const auto &q_IJchi0 : chi0q.second)
            {
                const int iq = std::distance(klist.begin(),
                                             std::find(klist.begin(), klist.end(), q_IJchi0.first));
                for (const auto &I_Jchi0 : q_IJchi0.second)
                {
                    const auto &I = I_Jchi0.first;
                    for (const auto &J_chi0 : I_Jchi0.second)
                    {
                        const auto &J = J_chi0.first;
                        sprintf(fn, "chi0fq_ifreq_%d_iq_%d_I_%zu_J_%zu_id_%d.mtx", ifreq, iq, I, J,
                                mpi_comm_global_h.myid);
                        print_complex_matrix_mm(J_chi0.second, Params::output_dir + "/" + fn,
                                                1e-15);
                    }
                }
            }
        }
    }

    // NOTE: Cs is cleaned up.
    // This means that the behavior will be undefined if this function is called again
    Cs_data.clear();
    LIBRPA::utils::release_free_mem();

    mpi_comm_global_h.barrier();
    Profiler::start("EcRPA", "Compute RPA correlation Energy");
    CorrEnergy corr;
    if (Params::use_scalapack_ecrpa &&
        (LIBRPA::parallel_routing == LIBRPA::ParallelRouting::ATOM_PAIR ||
         LIBRPA::parallel_routing == LIBRPA::ParallelRouting::LIBRI))
    {
        if (meanfield.get_n_kpoints() == 1)
        {
            print_complex_matrix_mm(*Vq.at(0).at(0).at(qlist[0]), "coulmat", 0.);
            corr = compute_RPA_correlation_blacs_2d_gamma_only(chi0, Vq);
        }
        else
        {
            corr = compute_RPA_correlation_blacs_2d(chi0, Vq);
        }
    }
    else
        corr = compute_RPA_correlation(chi0, Vq);

    rpa_corr = corr.value;

    cout << qlist << "\n";
    for (const auto &irk_corr : corr.qcontrib)
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
