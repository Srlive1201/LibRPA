#include "lib_main.h"
int main(int argc, char **argv)
{
    READ_AIMS_BAND("band_out");
    READ_AIMS_STRU("stru_out");
    READ_AIMS_Cs("Cs_data.txt");
    //READ_AIMS_Vq("Coulomb_q_mat");
    READ_AIMS_Vq_real("real_coulomb.txt"); 
    READ_AIMS_EIGENVECTOR("KS_eigenvector",cal_chi0.wfc_k); 
    READ_AIMS_EIGENVECTOR("KS_eigenvector",ap_chi0.wfc_k);
    cout<<" argv:"<<argv[1]<<endl;
    if(argv[1][0]=='1')
        ap_chi0.main();
    if(argv[1][0]=='2')
        cal_chi0.chi0_main(); 
    
    return 0;
}
