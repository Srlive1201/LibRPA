#include "ri.h"

#include <memory.h>

#include "envs_mpi.h"
#include "utils_io.h"
#include "params.h"
#include "pbc.h"

int n_irk_points;
int natom;
int ncell;
std::vector<Vector3_Order<double>> irk_points;
map<Vector3_Order<double>, double> irk_weight;
map<atom_t, size_t> atom_nw;
map<atom_t, size_t> atom_mu;
map<atom_t, size_t> atom_nw_loc;
map<atom_t, size_t> atom_mu_loc;
vector<size_t> atom_mu_part_range;
int N_all_mu;

vector<atpair_t> tot_atpair;
vector<atpair_t> local_atpair;
vector<atpair_t> tot_atpair_ordered;

atpair_R_mat_t Cs;
atpair_k_cplx_mat_t Vq;
atpair_k_cplx_mat_t Vq_cut;
map<Vector3_Order<double>, ComplexMatrix> Vq_block_loc;
map<Vector3_Order<double>, ComplexMatrix> Vq_cut_block_loc;

int atom_iw_loc2glo(const int &atom_index, const int &iw_lcoal)
{
    int nb = 0;
    for (int ia = 0; ia != atom_index; ia++)
        nb += atom_nw[ia];
    return iw_lcoal + nb;
}

int atom_mu_loc2glo(const int &atom_index, const int &mu_lcoal)
{
    // int nb = 0;
    // for (int ia = 0; ia != atom_index; ia++)
    //     nb += atom_mu[ia];
    // return mu_lcoal + nb;
    return atom_mu_part_range[atom_index]+mu_lcoal;
}

int atom_mu_glo2loc(const int &glo_index, int &mu_index)
{
    for(int I=0;I!=atom_mu.size();I++)
    {
        if((mu_index=glo_index-atom_mu_part_range[I])<atom_mu[I])
            return I;
    }
    throw invalid_argument("invalid glo_index");
}

vector<int> get_part_range()
{
    vector<int> part_range;
    part_range.resize(atom_mu.size());
    part_range[0] = 0;

    int count_range = 0;
    for (int Mu = 0; Mu != atom_mu.size() - 1; Mu++)
    {
        count_range += atom_mu[Mu];
        part_range[Mu + 1] = count_range;
    }
    return part_range;
}

matrix reshape_Cs(size_t n1, size_t n2, size_t n3, const shared_ptr<matrix> &Csmat) //(n1*n2,n3) -> (n2,n1*n3)
{
    const auto length = sizeof(double) * n3;
    const auto n13 = n1 * n3;
    const double *m_ptr = (*Csmat).c;
    matrix m_new(n2, n1 * n3, false);
    for (size_t i1 = 0; i1 != n1; ++i1)
    {
        double *m_new_ptr = m_new.c + i1 * n3;
        for (size_t i2 = 0; i2 != n2; ++i2, m_ptr += n3, m_new_ptr += n13)
            memcpy(m_new_ptr, m_ptr, length);
    }
    return m_new;
}

matrix reshape_mat(const size_t n1, const size_t n2, const size_t n3, const matrix &mat) //(n1,n2*n3) -> (n2,n1*n3)
{
    const auto length = sizeof(double) * n3;
    const auto n13 = n1 * n3;
    const double *m_ptr = mat.c;
    matrix m_new(n2, n1 * n3, false);
    for (size_t i1 = 0; i1 != n1; ++i1)
    {
        double *m_new_ptr = m_new.c + i1 * n3;
        for (size_t i2 = 0; i2 != n2; ++i2, m_ptr += n3, m_new_ptr += n13)
            memcpy(m_new_ptr, m_ptr, length);
    }
    return m_new;
}

matrix reshape_mat_21(const size_t n1, const size_t n2, const size_t n3, const matrix &mat) //(n1,n2*n3) -> (n1*n2,n3)
{
    const auto length = sizeof(double) * n1 * n2 * n3;
    const auto n13 = n1 * n2 * n3;
    const double *m_ptr = mat.c;
    matrix m_new(n1 * n2, n3, false);
    double *m_new_ptr = m_new.c;
    memcpy(m_new_ptr, m_ptr, length);
    return m_new;
}

void init_N_all_mu()
{
    using LIBRPA::envs::mpi_comm_global;
    using LIBRPA::envs::mpi_comm_global_h;
    using LIBRPA::utils::lib_printf;
    
    lib_printf("begin init_N_all_mu   atom_mu.size: %d  myid: %d\n",atom_mu.size(), mpi_comm_global_h.myid);
    if(Params::DFT_software == "ABACUS")
    {
        mpi_comm_global_h.barrier();
        vector<size_t> loc_mu(natom);
        vector<size_t> loc_nw(natom);
        for(const auto &nw_p:atom_nw)
        {
            loc_nw[nw_p.first]=nw_p.second;
        }
        for(const auto &mu_p:atom_mu)
        {
            loc_mu[mu_p.first]=mu_p.second;
        }
        LIBRPA::utils::lib_printf(" mid init_N_all_mu\n");
        vector<size_t> glo_nw(natom);
        vector<size_t> glo_mu(natom);
        LIBRPA::utils::lib_printf(" mid init_N_all_mu myid: %d\n",mpi_comm_global_h.myid);
        MPI_Allreduce(loc_nw.data(),glo_nw.data(),natom,MPI_UNSIGNED_LONG_LONG,MPI_MAX,mpi_comm_global);
        MPI_Allreduce(loc_mu.data(),glo_mu.data(),natom,MPI_UNSIGNED_LONG_LONG,MPI_MAX,mpi_comm_global);

        atom_nw.clear();
        atom_mu.clear();
        for(int i =0; i!=natom;i++ )
        {
            LIBRPA::utils::lib_printf("I: %d ,  glo_nw: %d",i,glo_nw[i]);
            LIBRPA::utils::lib_printf("I: %d ,  glo_mu: %d",i,glo_mu[i]);
            atom_nw.insert(pair<atom_t,size_t>(i,glo_nw[i]));
            atom_mu.insert(pair<atom_t,size_t>(i,glo_mu[i]));
        }
        allreduce_atp_aux();
    }
    
    atom_mu_part_range.resize(atom_mu.size());
    atom_mu_part_range[0]=0;
    for(int I=1;I!=atom_mu.size();I++)
        atom_mu_part_range[I]=atom_mu.at(I-1)+atom_mu_part_range[I-1];
    
    N_all_mu=atom_mu_part_range[natom-1]+atom_mu[natom-1];
    // LIBRPA::utils::lib_printf("end init_N_all_mu, atom_mu.size: %d\n",atom_mu.size());
  //  MPI_Barrier(mpi_comm_global);
}

void allreduce_atp_aux()
{
    using LIBRPA::envs::mpi_comm_global;
    using LIBRPA::envs::mpi_comm_global_h;

    for(int I=0;I!=natom;I++)
        for(int J=0;J!=natom;J++)
        {
            auto nbasI=atom_nw[I];
            auto nbasJ=atom_nw[J];
            auto nauxI=atom_mu[I];
            int cs_size=nbasI*nbasJ*nauxI;

            int Rsize=Cs[I][J].size();
            int loc_count=3* Rsize;

            std::vector<int> send_R_data;
            std::vector<int> all_R_counts(mpi_comm_global_h.nprocs, 0);

            all_R_counts[mpi_comm_global_h.myid] = loc_count;

            for (const auto& Rp : Cs[I][J]) 
            {   
                auto R=Rp.first;
                // if(I==0 && J==0)
                //     LIBRPA::utils::lib_printf(" send R  Cs  myid : %d   I: %d, J: %d, R:(%d, %d, %d)\n",mpi_comm_global_h.myid, I,J,R.x,R.y,R.z);
                send_R_data.push_back(Rp.first.x);
                send_R_data.push_back(Rp.first.y);
                send_R_data.push_back(Rp.first.z);
            }
            MPI_Allgather(&loc_count, 1, MPI_INT, all_R_counts.data(), 1, MPI_INT, mpi_comm_global);
            std::vector<int> displs(mpi_comm_global_h.nprocs);
            int total_R_count = 0;
            for (int i = 0; i < mpi_comm_global_h.nprocs; ++i) {
                displs[i] = total_R_count;
                total_R_count +=  all_R_counts[i];
            }
            std::vector<int> all_R_data(total_R_count);

            MPI_Allgatherv(send_R_data.data(), loc_count, MPI_INT, all_R_data.data(), all_R_counts.data(), displs.data(), MPI_INT, mpi_comm_global);
            
            int nR=total_R_count/3;
            for(int iR=0;iR!=nR;iR++)
            {
                Vector3_Order<int> Rv(all_R_data[iR*3],all_R_data[iR*3+1],all_R_data[iR*3+2]);
                // if(I==0 && J==0)
                //     LIBRPA::utils::lib_printf(" Rv  Cs  myid : %d   I: %d, J: %d, R:(%d, %d, %d)\n",mpi_comm_global_h.myid, I,J,Rv.x,Rv.y,Rv.z);
                shared_ptr<matrix> cs_ptr = make_shared<matrix>();
                cs_ptr->create(nbasI*nbasJ, nauxI);
                matrix loc_cs(nbasI*nbasJ, nauxI);
                if(Cs[I][J].count(Rv))
                    memcpy(loc_cs.c,(*Cs[I][J].at(Rv)).c,sizeof(double)*cs_size);
                mpi_comm_global_h.allreduce_matrix(loc_cs,(*cs_ptr));
                if(!Cs[I][J].count(Rv))
                    Cs[I][J][Rv] = cs_ptr;
            }
        }
}

void allreduce_2D_coulomb_to_atompair(map<Vector3_Order<double>, ComplexMatrix> &Vq_loc, atpair_k_cplx_mat_t &coulomb_mat,double threshold )
{
    using LIBRPA::envs::mpi_comm_global_h;
    using LIBRPA::utils::lib_printf;

    lib_printf("Begin allreduce_2D_coulomb_to_atompair!\n");

    map<Vector3_Order<double>, ComplexMatrix> Vq_glo;
    for(auto &qvec:klist_ibz)
    {
        Vq_glo[qvec].create(N_all_mu,N_all_mu);
        if(!Vq_block_loc.count(qvec))
        {
            Vq_block_loc[qvec].create(N_all_mu, N_all_mu);
        }
        
        mpi_comm_global_h.allreduce_ComplexMatrix(Vq_block_loc[qvec],Vq_glo[qvec]);
        // print_complex_matrix("vq_glo", Vq_glo[qvec]);
    }
    size_t vq_save = 0;
    size_t vq_discard = 0;
    for (auto &vf_p : Vq_glo)
    {
        auto qvec = vf_p.first;
        // cout << "Qvec:" << qvec << endl;
        
        for (int I = 0; I != atom_mu.size(); I++)
            for (int J = 0; J != atom_mu.size(); J++)
            {
                if (I > J)
                    continue;
                shared_ptr<ComplexMatrix> vq_ptr = make_shared<ComplexMatrix>();
                vq_ptr->create(atom_mu[I], atom_mu[J]);
                // vq_ptr_tran->create(atom_mu[J],atom_mu[I]);
                // cout << "I J: " << I << "  " << J << "   mu,nu: " << atom_mu[I] << "  " << atom_mu[J] << endl;
                for (int i_mu = 0; i_mu != atom_mu[I]; i_mu++)
                {

                    for (int i_nu = 0; i_nu != atom_mu[J]; i_nu++)
                    {
                        //(*vq_ptr)(i_mu, i_nu) = vf_p.second(atom_mu_loc2glo(J, i_nu), atom_mu_loc2glo(I, i_mu)); ////for aims
                        (*vq_ptr)(i_mu, i_nu) = vf_p.second(atom_mu_loc2glo(I, i_mu), atom_mu_loc2glo(J, i_nu)); // for abacus
                    }
                }

                // if (I == J)
                // {
                //     (*vq_ptr).set_as_identity_matrix();
                // }

                if ((*vq_ptr).real().absmax() >= threshold)
                {
                    coulomb_mat[I][J][qvec] = vq_ptr;
                    vq_save++;
                }
                else
                {
                    vq_discard++;
                }
            }
    }

    lib_printf("End allreduce_2D_coulomb_to_atompair!\n");
}

void allreduce_atp_coulomb( atpair_k_cplx_mat_t &coulomb_mat )
{
    using LIBRPA::envs::mpi_comm_global;
    using LIBRPA::envs::mpi_comm_global_h;
    using LIBRPA::utils::lib_printf;

    printf("Begin allreduce_atp_coulomb!\n");
    // mpi_comm_world_h.barrier();
    map<Vector3_Order<double>, ComplexMatrix> Vq_glo;
    for(auto &qvec:klist_ibz)
    {
        for(int I=0;I!=atom_mu.size();I++)
        {
            for(int J=I;J!=atom_mu.size();J++)
            {
                shared_ptr<ComplexMatrix> vq_recv_ptr = make_shared<ComplexMatrix>();
                vq_recv_ptr->create(atom_mu[I], atom_mu[J]);
                ComplexMatrix tmp_send(atom_mu[I], atom_mu[J]);
                // ComplexMatrix tmp_recv(atom_mu[I], atom_mu[J]);
                if(coulomb_mat[I][J].count(qvec)!=0)
                    tmp_send=*coulomb_mat.at(I).at(J).at(qvec);
                    
                printf("I: %d, J: %d, mu:%d, nu:%d, myid: %d\n",I,J,atom_mu[I], atom_mu[J],mpi_comm_global_h.myid);
                mpi_comm_global_h.allreduce_ComplexMatrix(tmp_send,*vq_recv_ptr);
                if(coulomb_mat[I][J].count(qvec)==0)
                    coulomb_mat[I][J][qvec]=vq_recv_ptr;
                printf("end I: %d, J: %d,  myid: %d\n",I,J,mpi_comm_global_h.myid);
                
                
                // coulomb_mat[I][J][qvec]=vq_recv_ptr;
                    
            }
        }
        
    }

    printf("End allreduce_atp_coulomb!\n");
    // mpi_comm_world_h.barrier();
}