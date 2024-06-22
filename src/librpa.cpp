/**
 * @file librpa.cpp
 * @author Rong Shi (srlive@mail.ustc.edu.cn)
 * @brief 
 * @version 0.1
 * @date 2023-08-18
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#include "librpa.h"
#include "matrix_m.h"
#include "meanfield.h"
#include <stdlib.h>
#include "pbc.h"
#include "constants.h"
#include "ri.h"
#include "librpa_main.h"
#include "params.h"
#include "parallel_mpi.h"
using LIBRPA::mpi_comm_world_h;
void test_interface(int test_int, char* test_str)
{
    cout<<"In LibRPA: "<<test_int<<"   "<<test_str<<endl;
}
void set_dimension(int nspins, int nkpts, int nstates, int nbasis,int natoms)
{
    cout<<"In LibRPA nspin: "<<nspins<<"  nkpt: "<<nkpts<<endl;
    printf("In LibRPA nspin: %d,  nkpts: %d nstate: %d, nbasis: %d, natom: %d\n",nspins,nkpts,nstates,nbasis,natoms);
    meanfield.set(nspins, nkpts, nstates, nbasis);
    natom=natoms;
    
}

void set_wg_ekb_efermi(int nspins, int nkpts, int nstates, double* wg, double* ekb, double efermi)
{
    meanfield.get_efermi() =  efermi *2.0 ;
    auto & eskb = meanfield.get_eigenvals();
    auto & swg = meanfield.get_weight();
    int length_kb=nkpts*nstates;
    for(int is= 0;is!=nspins;is++)
    {
        memcpy(eskb[is].c,ekb+length_kb*is,length_kb*sizeof(double));
        memcpy(swg[is].c,wg+length_kb*is,length_kb*sizeof(double));
        eskb[is]*=2;
        swg[is]*=(1.0/nkpts);
        
    }
    // for(int is=0;is!=nspins;is++)
    // {
    //     print_matrix(" eskb",eskb[is]);
    //     print_matrix(" swg",swg[is]);
    // }
}

void set_ao_basis_wfc(int is, int ik, double* wfc_real, double* wfc_imag)
{
    printf("is: %d, ik: %d\n",is,ik);
    // int length_ib_iw=meanfield.get_n_bands()*meanfield.get_n_aos();
    // vector<double> vec_wfc_real(wfc_real,wfc_real+length_ib_iw);
    // vector<double> vec_wfc_imag(wfc_imag,wfc_imag+length_ib_iw);
    auto & wfc = meanfield.get_eigenvectors();
    for(int i=0;i!=meanfield.get_n_bands()*meanfield.get_n_aos();i++)
    {
        printf("In ao wfc: %f, %f\n",wfc_real[i],wfc_imag[i]);
        wfc.at(is).at(ik).c[i]=complex<double>(wfc_real[i],wfc_imag[i]);
    }
    // print_complex_matrix("wfc_isk", wfc.at(is).at(ik));
}

void set_latvec_and_G(double* lat_mat, double* G_mat)
{
    latvec.e11 = lat_mat[0];
    latvec.e12 = lat_mat[1];
    latvec.e13 = lat_mat[2];

    latvec.e21 = lat_mat[3];
    latvec.e22 = lat_mat[4];
    latvec.e23 = lat_mat[5];

    latvec.e31 = lat_mat[6];
    latvec.e32 = lat_mat[7];
    latvec.e33 = lat_mat[8];
    
    latvec /= ANG2BOHR;
    lat_array[0] = {latvec.e11,latvec.e12,latvec.e13};
    lat_array[1] = {latvec.e21,latvec.e22,latvec.e23};
    lat_array[2] = {latvec.e31,latvec.e32,latvec.e33};

    G.e11 = G_mat[0];
    G.e12 = G_mat[1];
    G.e13 = G_mat[2];
   
    G.e21 = G_mat[3];
    G.e22 = G_mat[4];
    G.e23 = G_mat[5];
 
    G.e31 = G_mat[6];
    G.e32 = G_mat[7];
    G.e33 = G_mat[8];

    G /= TWO_PI;
    G *= ANG2BOHR;
    
    printf(" LibRPA_lat : %f, %f, %f\n",latvec.e11,latvec.e12,latvec.e13);
    //latvec.print();
    //G.print();
}

void set_kgrids_kvec_tot(int nk1, int nk2, int nk3, double* kvecs)
{
    kv_nmp[0] = nk1;
    kv_nmp[1] = nk2;
    kv_nmp[2] = nk3;

    for(int ik=0;ik!=meanfield.get_n_kpoints();ik++)
    {
        double kx=kvecs[ik*3];
        double ky=kvecs[ik*3+1];
        double kz=kvecs[ik*3+2];
        // printf("ik: %d, (%f, %f, %f)\n", ik,kx,ky,kz);
        Vector3_Order<double> kvec_tmp{kx,ky,kz};
        kvec_tmp*=(ANG2BOHR / TWO_PI);
        klist.push_back(kvec_tmp);
        kfrac_list.push_back(latvec * kvec_tmp);
        
    }
}

void set_ibz2bz_index_and_weight(const int nk_irk, const int* ibz2bz_index, const double* wk_irk)
{
   // printf(" nks_irk: %d\n",nk_irk);
    for(int ik_ibz=0;ik_ibz!=nk_irk;ik_ibz++)
    {
        Vector3_Order<double> kvec_ibz=klist[ibz2bz_index[ik_ibz]];
        klist_ibz.push_back(kvec_ibz);
        irk_weight.insert(pair<Vector3_Order<double>, double>(kvec_ibz, wk_irk[ik_ibz]));
        // printf("ibz2bz:  %d   kvec_ibz:( %f, %f, %f)\n",ibz2bz_index[ik_ibz],kvec_ibz.x,kvec_ibz.y,kvec_ibz.z);
        //printf("irk_weight: %f\n",irk_weight[kvec_ibz]);
        // for (int ik = 0; ik != meanfield.get_n_kpoints(); ik++)
        // {
        //     if(klist[ik]==kvec_ibz)
        //         map_irk_ks[kvec_ibz].push_back(klist[ik]);
        // }
    }
}

void set_ao_basis_aux(int I, int J, int nbasis_i, int nbasis_j, int naux_mu, int* R, double* Cs_in)
{
    atom_nw.insert(pair<atom_t, int>(I, nbasis_i));
    atom_mu.insert(pair<atom_t, int>(I, naux_mu));
    Vector3_Order<int> box(R[0],R[1],R[2]);
        // cout<< ia1<<ia2<<box<<endl;
    // printf("Cs_in size: %zu",sizeof(Cs_in) / sizeof(Cs_in[0]));
    int cs_size=nbasis_i*nbasis_j*naux_mu;
    shared_ptr<matrix> cs_ptr = make_shared<matrix>();
    cs_ptr->create(nbasis_i * nbasis_j, naux_mu);
    memcpy((*cs_ptr).c,Cs_in,sizeof(double)*cs_size);
    //(*cs_ptr).c=Cs_in;
    Cs[I][J][box] = cs_ptr;
    // printf("Cs out:\n");
    // for(int i=0;i!=cs_size;i++)
    //     printf("   %f",(*Cs[I][J][box]).c[i]);
    
}

void set_aux_coulomb_k_atom_pair(int I, int J, int naux_mu, int naux_nu, int ik, double* Vq_real_in, double* Vq_imag_in)
{
    // printf("I,J,mu,nu: %d  %d  %d  %d\n",I,J, naux_mu,naux_nu);
    // printf("atom_mu nu: %d %d\n",atom_mu[I],atom_mu[J]);
    // printf("Vq threshold : %f\n",Params::vq_threshold);
    Vector3_Order<double> qvec(klist[ik]);
    //printf("qvec: %f,%f,%f\n",qvec.x,qvec.y,qvec.z);
    shared_ptr<ComplexMatrix> vq_ptr = make_shared<ComplexMatrix>();
    vq_ptr->create(naux_mu, naux_nu);
    // printf("vq_ptr_size: %d\n",(*vq_ptr).size);
    // printf(" vq_in 0 0 :  %f\n",Vq_real_in[0]);
    for (int i_mu = 0; i_mu != naux_mu; i_mu++)
    {
        for (int i_nu = 0; i_nu != naux_nu; i_nu++)
        {
           // printf("mu,nu:  %d %d, vq_real: %f,  vq_imag: %f\n",i_mu,i_nu,Vq_real_in[i_nu+i_mu*naux_nu],Vq_imag_in[i_nu+i_mu*naux_nu]);

            (*vq_ptr)(i_mu, i_nu)=complex<double>(Vq_real_in[i_nu+i_mu*naux_nu],Vq_imag_in[i_nu+i_mu*naux_nu]);
        }
    }
    if ((*vq_ptr).real().absmax() >= Params::vq_threshold)
    {
        Vq[I][J][qvec] = vq_ptr;
    }
    //print_complex_matrix("Vq",(*Vq[I][J][qvec]));
}

void set_aux_coulomb_k_2D_block(int ik, int max_naux, int mu_begin, int mu_end, int nu_begin, int nu_end, double* Vq_real_in, double* Vq_imag_in )
{
    int brow = mu_begin - 1;
    int erow = mu_end - 1;
    int bcol = nu_begin - 1;
    int ecol = nu_end - 1;
      
    Vector3_Order<double> qvec(klist[ik]);
    printf("qvec: %f,%f,%f\n",qvec.x,qvec.y,qvec.z);
    shared_ptr<ComplexMatrix> vq_ptr = make_shared<ComplexMatrix>();
    vq_ptr->create(max_naux, max_naux);

    if (!Vq_block_loc.count(qvec))
    {
        Vq_block_loc[qvec].create(max_naux, max_naux);
    }
    int ii=0;
    for (int i_mu = brow; i_mu <= erow; i_mu++)
    {
        for (int i_nu = bcol; i_nu <= ecol; i_nu++)
        {
            Vq_block_loc[qvec](i_mu, i_nu) = complex<double>(Vq_real_in[ii], Vq_imag_in[ii]); 
            ii+=1;
        }
    }
}

void set_librpa_params()
{
    Params::nfreq=12;
    Params::print();
}

void run_librpa_main(MPI_Comm comm_in, int is_fortran_comm)
{
    printf("begin run librpa\n");
    // std::ofstream outputFile("LibRPA_cout.txt");
    // std::streambuf* originalCoutBuffer = std::cout.rdbuf();
   // std::cout.rdbuf(outputFile.rdbuf());
    librpa_main(comm_in);
    
    //std::cout.rdbuf(originalCoutBuffer);
    //outputFile.close();
}