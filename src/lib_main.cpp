#include "lib_main.h"
#include <algorithm>
#include <malloc.h>
#include "scalapack_connector.h"
#include <set>
int main(int argc, char **argv)
{
    para_mpi.mpi_init(argc,argv);

    prof.add(0, "total", "Total");
    prof.add(1, "chi0_main", "Chi0 object");
    prof.add(2, "cal_Green_func_R_tau", "space-time Green's function");
    prof.add(2, "cal_Green_func",       "space-time Green's function");
    prof.add(2, "R_tau_routing", "Loop over R-tau");
    prof.add(2, "atom_pair_rouing", "Loop over atom pairs");
    prof.add(2, "LibRI_rouing", "Loop over LibRI");
    prof.add(3, "cal_chi0_element", "chi(tau,R,I,J)");
    prof.add(4, "X");
    prof.add(4, "O");
    prof.add(4, "N");
    prof.add(4, "Z");
    prof.add(4, "reshape_Cs", "reshape Cs");
    prof.add(4, "reshape_mat", "reshape mat");

    prof.start("total");

    int flag;
    parse_inputfile_to_params(input_filename, params);
    // NOTE: to comply with command line input, may be removed later
    InputFile inputf;
    auto parser = inputf.load(input_filename, false);
    parser.parse_int("nfreq", params.nfreq, stoi(argv[1]), flag);
    parser.parse_double("gf_R_threshold", params.gf_R_threshold, stod(argv[2]), flag);
    params.check_consistency();
    if (para_mpi.is_master())
        params.print();
    para_mpi.mpi_barrier();
    para_mpi.set_blacs_parameters();

    READ_AIMS_BAND("band_out", meanfield);
    if (para_mpi.get_myid() == 0)
    {
        cout << "Information of mean-field starting-point" << endl;
        cout << "| number of spins: " << meanfield.get_n_spins() << endl
             << "| number of k-points: " << meanfield.get_n_kpoints() << endl
             << "| number of bands: " << meanfield.get_n_bands() << endl
             << "| number of NAOs: " << meanfield.get_n_aos() << endl
             << "| efermi (Ha): " << meanfield.get_efermi() << endl;

        double emin, emax;
        meanfield.get_E_min_max(emin, emax);
        cout << "| Minimal transition energy (Ha): " << emin << endl
             << "| Maximal transition energy (Ha): " << emax << endl << endl;
    }

    READ_AIMS_STRU(meanfield.get_n_kpoints(), "stru_out");
    Vector3_Order<int> period {kv_nmp[0], kv_nmp[1], kv_nmp[2]};
    auto Rlist = construct_R_grid(period);

    if (para_mpi.get_myid() == 0)
    {
        cout << "Lattice vectors (Bohr)" << endl;
        latvec.print();
        cout << "Reciprocal lattice vectors (2PI Bohr^-1)" << endl;
        G.print();
        printf("kgrids: %2d %2d %2d\n", kv_nmp[0], kv_nmp[1], kv_nmp[2]);
        cout << "k-points read (Cartisian in 2Pi Bohr^-1 | fractional):" << endl;
        for (int ik = 0; ik != meanfield.get_n_kpoints(); ik++)
        {
            printf("ik %4d: %10.7f %10.7f %10.7f | %10.7f %10.7f %10.7f\n",
                   ik+1, kvec_c[ik].x, kvec_c[ik].y, kvec_c[ik].z,
                   kfrac_list[ik].x, kfrac_list[ik].y, kfrac_list[ik].z);
        }
        cout << "R-points to compute:" << endl;
        for (int iR = 0; iR != Rlist.size(); iR++)
        {
            printf("iR %4d: %3d %3d %3d\n", iR+1, Rlist[iR].x, Rlist[iR].y, Rlist[iR].z);
        }
        cout << endl;
    }

    const int Rt_num = Rlist.size() * params.nfreq;

    READ_AIMS_EIGENVECTOR("./", meanfield);

    READ_AIMS_Cs("./", params.cs_threshold);

    tot_atpair = generate_atom_pair_from_nat(natom, false);
    if (para_mpi.is_master())
        cout << "| Natoms: " << natom << "   tot_atpairs:  " << tot_atpair.size() << endl;
    // barrier to wait for information print on master process
    para_mpi.mpi_barrier();
    
    para_mpi.set_chi_parallel_type(tot_atpair.size(),Rt_num,params.use_libri_chi0);
    //para_mpi.chi_parallel_type=Parallel_MPI::parallel_type::ATOM_PAIR;
    //vector<atpair_t> local_atpair;
    if(para_mpi.chi_parallel_type==Parallel_MPI::parallel_type::ATOM_PAIR)
    {
        local_atpair = dispatch_vector(tot_atpair, para_mpi.get_myid(), para_mpi.get_size(), true);
        printf("| process %d , local_atom_pair size:  %zu\n", para_mpi.get_myid(), local_atpair.size());
        // for(auto &ap:local_atpair)
        //     printf(" |process %d , local_atom_pair:  %d,  %d\n",para_mpi.get_myid(),ap.first,ap.second);
        READ_Vq_Row("./", "coulomb_mat", params.vq_threshold, Vq, local_atpair);
    }
    else
    {
        READ_Vq_Full("./", "coulomb_mat", params.vq_threshold, Vq); 
        local_atpair = get_atom_pair(Vq);
    }
    // malloc_trim(0);
    // para_mpi.mpi_barrier();
    // if(para_mpi.is_master())
    //     system("free -m");
    erase_Cs_from_local_atp(Cs,local_atpair);
    // malloc_trim(0);
    // para_mpi.mpi_barrier();
    // if(para_mpi.is_master())
    //     system("free -m");
    //READ_Vq_Row("./", "coulomb_mat", params.vq_threshold, Vq, local_atpair);
    /* if(argv[1][0]=='0') */
    /*     ap_chi0.chi0_main(argv[1],argv[2]);  */
    /* else */
        /* cal_chi0.chi0_main(argv[1],argv[2]);  */
    /* return 0; */
    // try the new version
    Chi0 chi0(meanfield, klist, params.nfreq);
    chi0.gf_R_threshold = params.gf_R_threshold;

    // build ABF IJ and qlist from Vq
    
    vector<Vector3_Order<double>> qlist;

    for ( auto q_weight: irk_weight)
    {
        qlist.push_back(q_weight.first);
    }

    chi0.build(Cs, Rlist, period, local_atpair, qlist,
               TFGrids::get_grid_type(params.tfgrids_type), true);
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
                        sprintf(fn, "chi0fq_ifreq_%d_iq_%d_I_%zu_J_%zu_id_%d.mtx", ifreq, iq, I, J, para_mpi.get_myid());
                        // print_complex_matrix_mm(J_chi0.second, fn, 1e-15);
                    }
                }
            }
        }
    }
    //malloc_trim(0);
    
    // para_mpi.mpi_barrier();
    // if(para_mpi.is_master())
    //     system("free -m");
    for(auto &Cp:Cs)
    {
        Cs.erase(Cp.first);
       
    }
    malloc_trim(0);
    // RPA total energy
    if ( params.task == "rpa" )
    {
        CorrEnergy corr;
        if (params.use_scalapack_ecrpa && para_mpi.chi_parallel_type==Parallel_MPI::parallel_type::ATOM_PAIR)
            corr = compute_RPA_correlation_blacs(chi0, Vq);
        else
            corr = compute_RPA_correlation(chi0, Vq);

        if (para_mpi.get_myid() == 0)
        {
            printf("RPA correlation energy (Hartree)\n");
            printf("| Weighted contribution from each k:\n");
            for (const auto& q_ecrpa: corr.qcontrib)
            {
                cout << "| " << q_ecrpa.first << ": " << q_ecrpa.second << endl;
            }
            printf("| Total EcRPA: %18.9f\n", corr.value.real());
            if (std::abs(corr.value.imag()) > 1.e-3)
                printf("Warning: considerable imaginary part of EcRPA = %f\n", corr.value.imag());
        }
    }
    else if ( params.task == "g0w0" )
    {
        READ_Vq_Full("./", "coulomb_cut_", params.vq_threshold, Vq_cut); 
        const auto VR = FT_Vq(Vq_cut, Rlist, true);
        auto exx = LIBRPA::Exx(meanfield, klist);
        exx.build_exx_orbital_energy(Cs, Rlist, period, VR);
        // const auto Wc_freq_q = compute_Wc_freq_q(chi0, Vq, Vq_cut);
        // const auto Wc_tau_R = CT_FT_Wc_freq_q(Wc_freq_q, chi0.tfg, Rlist);
    }
    else if ( params.task == "exx" )
    {
        READ_Vq_Full("./", "coulomb_cut_", params.vq_threshold, Vq_cut); 
        const auto VR = FT_Vq(Vq_cut, Rlist, true);
        auto exx = LIBRPA::Exx(meanfield, kfrac_list);
        exx.build_exx_orbital_energy(Cs, Rlist, period, VR);
        // FIXME: need to reduce first when MPI is used
        // NOTE: may extract to a common function
        if (para_mpi.is_master())
        {
            for (int isp = 0; isp != meanfield.get_n_spins(); isp++)
            {
                printf("Spin channel %1d\n", isp+1);
                for (int ik = 0; ik != meanfield.get_n_kpoints(); ik++)
                {
                    cout << "k-point " << ik + 1 << ": " << kvec_c[ik] << endl;
                    printf("%-4s  %-10s  %-10s\n", "Band", "e_exx (Ha)", "e_exx (eV)");
                    for (int ib = 0; ib != meanfield.get_n_bands(); ib++)
                        printf("%4d  %10.5f  %10.5f\n", ib+1, exx.Eexx[isp][ik][ib], HA2EV * exx.Eexx[isp][ik][ib]);
                    printf("\n");
                }
            }
        }
    }

    else if ( params.task == "rpa_force" )
    {
        CorrEnergy corr;
        if (params.use_scalapack_ecrpa)
            corr = compute_RPA_correlation_blacs(chi0, Vq);
        else
            corr = compute_RPA_correlation(chi0, Vq);
       cout << "energy in force " << endl; 

        if (para_mpi.get_myid() == 0)
        {
            printf("RPA correlation energy (Hartree)\n");
            printf("| Weighted contribution from each k:\n");
            for (const auto& q_ecrpa: corr.qcontrib)
            {
                cout << "| " << q_ecrpa.first << ": " << q_ecrpa.second << endl;
            }
            printf("| Total EcRPA: %18.9f\n", corr.value.real());
            if (std::abs(corr.value.imag()) > 1.e-3)
                printf("Warning: considerable imaginary part of EcRPA = %f\n", corr.value.imag());
        }

	       
    }

    prof.stop("total");
    if(para_mpi.get_myid()==0)
        prof.display();

    return 0;
}
