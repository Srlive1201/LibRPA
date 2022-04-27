#include "lib_main.h"

int main(int argc, char **argv)
{
    prof.add(0, "total", "Total");
    prof.add(1, "chi0_main", "Chi0 object");
    prof.add(2, "cal_Green_func_R_tau", "space-time Green's function");
    prof.add(2, "R_tau_routing", "Loop over R-tau");
    prof.add(2, "atom_pair_rouing", "Loop over atom pairs");
    prof.add(3, "cal_chi0_element", "chi(tau,R,I,J)");
    prof.add(4, "X");
    prof.add(4, "O");
    prof.add(4, "N");
    prof.add(4, "Z");
    prof.add(4, "reshape_Cs", "reshape Cs");
    prof.add(4, "reshape_mat", "reshape mat");

    prof.start("total");

    para_mpi.mpi_init(argc,argv);
    cout<<"  MPI process_num: "<<para_mpi.get_myid()<<endl;
    READ_AIMS_BAND("band_out");
    READ_AIMS_STRU("stru_out");
    READ_AIMS_EIGENVECTOR("./",cal_chi0.wfc_k); 
    READ_AIMS_EIGENVECTOR("./",ap_chi0.wfc_k);
    READ_AIMS_Cs("./");
    READ_AIMS_Vq("./"); 
    if(argv[1][0]=='0')
        ap_chi0.chi0_main(argv[1],argv[2]); 
    else
        cal_chi0.chi0_main(argv[1],argv[2]); 
   
    prof.stop("total");
    prof.display();

    return 0;
}
