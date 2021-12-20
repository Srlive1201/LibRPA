#include "lib_main.h"
#include "global_class.h"
int main(int argc, char **argv)
{
    para_mpi.mpi_init(argc,argv);
    cout<<"  MPI process_num: "<<para_mpi.get_myid()<<endl;
    READ_AIMS_BAND("band_out");
    READ_AIMS_STRU("stru_out");
    READ_AIMS_EIGENVECTOR("./",cal_chi0.wfc_k); 
    READ_AIMS_EIGENVECTOR("./",ap_chi0.wfc_k);
    READ_AIMS_Cs("./");
    READ_AIMS_Vq("./"); 
    if(argv[1][0]=='0')
        ap_chi0.main();
    else
        cal_chi0.chi0_main(argv[1]); 
    
    return 0;
}
