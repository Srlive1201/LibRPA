#include "lib_main.h"
#include "scalapack_connector.h"
int main(int argc, char **argv)
{
    para_mpi.mpi_init(argc,argv);
    cout<<"  MPI process_num: "<<para_mpi.get_myid()<<endl;

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
    auto parser = inputf.load(input_filename, false);
    parser.parse_string("task", task, "rpa", flag);
    parser.parse_double("cs_threshold", cs_threshold, 1e-6, flag);
    parser.parse_double("vq_threshold", vq_threshold, 1e-6, flag);
    parser.parse_double("sqrt_coulomb_threshold", sqrt_coulomb_threshold, 1e-4, flag);

    para_mpi.set_blacs_parameters();
    
    READ_AIMS_BAND("band_out", meanfield);
    READ_AIMS_STRU("stru_out");
    READ_AIMS_EIGENVECTOR("./", meanfield);

    READ_AIMS_Cs("./", cs_threshold);
    READ_AIMS_Vq("./", "coulomb_mat_", vq_threshold, Vq); 
    /* if(argv[1][0]=='0') */
    /*     ap_chi0.chi0_main(argv[1],argv[2]);  */
    /* else */
        /* cal_chi0.chi0_main(argv[1],argv[2]);  */
    /* return 0; */
    // try the new version
    Chi0 chi0(meanfield, klist, stoi(argv[1]));
    parser.parse_double("gf_threshold", chi0.gf_R_threshold, stod(argv[2]), flag);

    Vector3_Order<int> period {kv_nmp[0], kv_nmp[1], kv_nmp[2]};
    auto Rlist = construct_R_grid(period);
    cout << "Rlist to compute:" << endl;
    for (int iR = 0; iR != Rlist.size(); iR++)
    {
        printf("iR %5d: R = (%2d, %2d, %2d)\n", iR, Rlist[iR].x, Rlist[iR].y, Rlist[iR].z);
    }
    // build ABF IJ and qlist from Vq
    vector<atpair_t> atpairs_ABF = get_atom_pair(Vq);
    vector<Vector3_Order<double>> qlist;

    for ( auto q_weight: irk_weight)
    {
        qlist.push_back(q_weight.first);
    }
    chi0.build(Cs, Rlist, period, atpairs_ABF, qlist, TFGrids::GRID_TYPES::Minimax, true);
    // RPA total energy

    if ( task == "rpa" )
        compute_RPA_correlation(chi0, Vq);
    else if ( task == "g0w0" )
    {
        READ_AIMS_Vq("./", "coulomb_cut_", vq_threshold, Vq_cut); 
        const auto Wc_freq_q = compute_Wc_freq_q(chi0, Vq, Vq_cut);
        const auto Wc_tau_R = CT_FT_Wc_freq_q(Wc_freq_q, chi0.tfg, Rlist);
        const auto VR = FT_Vq(Vq_cut, Rlist);
    }

    prof.stop("total");
    if(para_mpi.get_myid()==0)
        prof.display();

    return 0;
}
