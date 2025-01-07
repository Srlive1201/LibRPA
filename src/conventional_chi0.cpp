#include "conventional_chi0.h"
#include "lapack_connector.h"
//#include "src_global/matrix.h"
#include<fstream>
#include<sstream>
#include<iomanip>
#include<cmath>
#include<cstring>
#include<omp.h>
//#include<iostream>
#include "parallel_mpi.h"
#include "envs_mpi.h"
#include "pbc.h"
#include "constants.h"
#include "params.h"
using LIBRPA::envs::mpi_comm_global_h;
using LIBRPA::ParallelRouting;
using LIBRPA::parallel_routing;
using LIBRPA::envs::ofs_myid;
using LIBRPA::utils::lib_printf;


void Conventional_Chi0::build(const Cs_LRI &Cs,
                   const vector<Vector3_Order<int>> &Rlist,
                   const Vector3_Order<int> &R_period,
                   const vector<atpair_t> &atpair_ABF,
                   const vector<Vector3_Order<double>> &qlist,
                   TFGrids::GRID_TYPES gt)
{
    Tot_band_num=mf.get_n_bands();
    cout<<"Tot_band_num: "<<Tot_band_num<<endl;
    cout<<" NSPIN= "<<mf.get_n_spins()<<endl;
    for(int is=0;is!=mf.get_n_spins();is++)
    {
        cout<<"tag_0"<<endl;
        int temp_occ=0;
        cout<<"1"<<endl;
        for(int i_band=0;i_band!=Tot_band_num;i_band++)
        {
            cout<<"2"<<endl;
            if(mf.get_weight()[0](is,i_band)>=1e-6)    //kv.nks == NSPIN ?
            //if(wg(is,i_band)>=1) 
            {    
                temp_occ++;
            }
        }
        occ_num[is]=temp_occ;
        //occ_num.insert(pair<size_t,size_t>(is,temp_occ));
        //unocc_num.insert(pair<size_t,size_t>(is,Tot_band_num-temp_occ));
        //tot_ia.insert(pair<size_t,size_t>(is,temp_occ*(Tot_band_num-temp_occ)));
        unocc_num[is]=Tot_band_num-occ_num[is];
        tot_ia[is]=occ_num[is]*unocc_num[is];
        cout<<"i_spin:"<<is<<"   occ_num: "<<occ_num[is]<<"     unocc_num:  "<<unocc_num[is]<<"    tot_ia"<<tot_ia[is]<<endl;
    }
    switch (gt)
    {
        case ( TFGrids::GRID_TYPES::Minimax):
        {
            double emin, emax;
            mf.get_E_min_max(emin, emax);
            tfg.generate_minimax(emin, emax);
            break;
        }
        case ( TFGrids::GRID_TYPES::EvenSpaced_TF ):
        {
            // WARN: only used to debug, adjust emin and interval manually
            tfg.generate_evenspaced_tf(0.005, 0.0, 0.005, 0.0);
            break;
        }
        default:
            throw invalid_argument("chi0 on requested grid is not implemented");
    }
    if (mpi_comm_global_h.is_root())
        tfg.show();
    mpi_comm_global_h.barrier();

    aims_C_tilde(Cs,qlist);

    complex<double> tot_cRPA(0.0,0.0);
    cout<<"contional chi0_freq:"<<endl;
    
    auto part_range=LIBRPA::atomic_basis_abf.get_part_range();
    
    char fn[80];
    for(int ifreq=0;ifreq!=tfg.size();ifreq++)
    {
        double freq=tfg.get_freq_nodes()[ifreq];
        double weight=tfg.get_freq_weights()[ifreq];
        map<size_t,map<size_t,map<size_t,ComplexMatrix>>> chi_omega_mat;
        for (size_t atom_pair = 0; atom_pair != atpair_ABF.size(); atom_pair++)
        {
            auto Mu = atpair_ABF[atom_pair].first;
            auto Nu = atpair_ABF[atom_pair].second;
           
            map<size_t,ComplexMatrix> chi_omega_mat_tmp(cal_ap_chi0_mat(freq,Mu,Nu,qlist));
            for (auto &pair : chi_omega_mat_tmp) {
                chi_omega_mat[pair.first][Mu][Nu] = std::move(pair.second);
            }
        }
        for(const auto &ikp:chi_omega_mat)
        {
            const auto ik=ikp.first;
            for(auto &I_p:ikp.second)
            {
                const auto I=I_p.first;
                for(auto &J_p:I_p.second)
                {
                    const auto J=J_p.first;
                    sprintf(fn, "chi0fq_ifreq_%d_iq_%d_I_%zu_J_%zu_id_%d.mtx", ifreq, ik, I, J, mpi_comm_global_h.myid);
                    print_complex_matrix_mm(J_p.second, Params::output_dir + "/" + fn, 1e-15);
                }
            }
        }
        map<size_t,map<size_t,ComplexMatrix>> pi_mat=cal_pi_mat(chi_omega_mat[0],freq);
     
        complex<double> part_cRPA=cal_cRPA_omega(pi_mat,freq,part_range);
        tot_cRPA+=part_cRPA*weight/TWO_PI;
    }
    cout<<"ap_cRPA:  "<<setprecision(8)<<tot_cRPA<<endl;

}

complex<double> Conventional_Chi0::cal_cRPA_omega(map<size_t,map<size_t,ComplexMatrix>> &pi_mat, const double &freq,const vector<atom_t> &part_range)
{
    cout<<" cal_cRPA_omega"<<endl;
    ComplexMatrix pi_2D(N_all_mu,N_all_mu);
    cout<<"N_all_mu: "<<N_all_mu<<endl;
    //complex<double> c_spin(0.5,0.0);
    for(const auto &I_pair:pi_mat)
    {
        const auto I=I_pair.first;
        for(const auto &J_pair:I_pair.second)
        {
            const auto J=J_pair.first;
            const auto &mat=J_pair.second;
            for(int i=0;i!=atom_mu[I];i++)
            {
                for(int j=0;j!=atom_mu[J];j++)
                {
                    //cout<<"   I J i j  "<<I<<J<<i<<j<<"     part_range_I: "<<part_range[I]<<"   mat.nr  nc:"<<mat.nr<<"  "<<mat.nc<<endl;
                    pi_2D(part_range[I]+i,part_range[J]+j)=mat(i,j);
                }
            }
        }
    }
    // if(abs(freq-0.000440961)<1e-6)
    // {
    //     print_complex_matrix( "pi_2D", pi_2D );    
    // }
    
    
    complex<double> rpa_for_omega(0.0,0.0);
	ComplexMatrix identity(N_all_mu,N_all_mu);
	ComplexMatrix identity_minus_pi(N_all_mu,N_all_mu);

	identity.set_as_identity_matrix();

	identity_minus_pi=identity-pi_2D;
    
	complex<double> det_for_rpa(1.0,0.0);
//	det_for_rpa=complexmatrix_det(range_all, identity_minus_pi);
	int info_LU =0;
	int *ipiv=new int[N_all_mu];
	LapackConnector::zgetrf(N_all_mu,N_all_mu,identity_minus_pi,N_all_mu,ipiv,&info_LU);
	for(int ib=0;ib!=N_all_mu;ib++)
	{
		if(ipiv[ib] !=(ib+1)) det_for_rpa=-det_for_rpa*identity_minus_pi(ib,ib);
		else det_for_rpa=det_for_rpa*identity_minus_pi(ib,ib);
	}
	delete[] ipiv;
	
	complex<double> trace_pi;
	complex<double> ln_det;
	ln_det=std::log(det_for_rpa);
	trace_pi=trace(pi_2D);
	rpa_for_omega=ln_det+trace_pi;
	cout<<"  freq: "<<freq<<"    "<<rpa_for_omega<<"     ln_det: "<<ln_det<<"      tr_pi: "<<trace_pi<<endl;
	return rpa_for_omega;
    
}

map<size_t,map<size_t,ComplexMatrix>> Conventional_Chi0::cal_pi_mat(map<size_t,map<size_t,ComplexMatrix>> &chi0_omega_mat, const double &freq)
{
    cout<<"cal_pi_mat"<<endl;
    map<size_t,map<size_t,ComplexMatrix>> pi_mat;
    
    for(const auto &I_pair:chi0_omega_mat)
    {
        const auto I=I_pair.first;
        //const auto adj_I=Abfs::get_adjs(I);
        const size_t I_aux_num=atom_mu[I];
        for(const auto &J_pair:I_pair.second)
        {
            const auto J=J_pair.first;
            const size_t J_aux_num=atom_mu[J];
            pi_mat[I][J].create(I_aux_num,J_aux_num);
            if(I!=J)
                pi_mat[J][I].create(J_aux_num,I_aux_num);
        }
    }

    for(auto &I_p:Vq)
    {
        const auto I=I_p.first;
        for(auto &Q_p:Vq)
        {
            const auto Q=Q_p.first;
            for(auto &J_p:I_p.second)
            {
                const auto J=J_p.first;
                for(auto &kvec_p:J_p.second)
                {
                    Vector3_Order<double> ik_vec=kvec_p.first;
                    auto Vq_mat=*kvec_p.second;
                    if(J<=Q)
                    {
                        //cout<<"   1";
                        pi_mat.at(I).at(Q)+=Vq_mat*chi0_omega_mat.at(J).at(Q);
                        // if(freq==first_freq)
                        // {
                        //     cout<<" IJQ: "<<I<<" "<<J<<" "<<Q<<endl;
                        //     print_complex_real_matrix("vq:",Vq_mat);
                        //     print_complex_real_matrix("chi0:",chi0_omega_mat.at(J).at(Q));
                        //     print_complex_real_matrix("pi_mat:",pi_mat.at(I).at(Q));
                        // }
                    }
                    else
                    {
                        //cout<<"   2";
                        pi_mat.at(I).at(Q)+=Vq_mat*transpose(chi0_omega_mat.at(Q).at(J),1);
                    }
                    
                    if(I!=J)
                    {
                        if(I<=Q)
                        {
                           // cout<<"   3"<<"  pi: "<<pi_k[freq][ik_vec][J][Q].nr<<pi_k[freq][ik_vec][J][Q].nc<<"  vq nr nc: "<<(transpose(Vq_mat,1)).nr<<"  "<<(transpose(Vq_mat,1)).nc<<"   chi0:"<<chi0_k[freq][ik_vec][I][Q].nr<<"  "<<chi0_k[freq][ik_vec][I][Q].nc<<endl;
                            pi_mat.at(J).at(Q)+=transpose(Vq_mat,1)*chi0_omega_mat.at(I).at(Q);
                        }
                        else
                        {
                            //cout<<"   4";
                            pi_mat.at(J).at(Q)+=transpose(Vq_mat,1)*transpose(chi0_omega_mat.at(Q).at(I),1);
                        }
                    }   
                }
            }
        }
    }
    return pi_mat;
}


map<size_t,ComplexMatrix> Conventional_Chi0::cal_ap_chi0_mat(const double &freq ,const size_t &Mu, const size_t &Nu,const vector<Vector3_Order<double>> &qlist)
{
    //cout<<"cal_chi0_mat"<<endl;
    const size_t I_aux_num=atom_mu[Mu];
    const size_t J_aux_num=atom_mu[Nu];

    map<size_t,ComplexMatrix> chi0_k_mat; 
    for(int ik=0;ik!=mf.get_n_kpoints();ik++)
        chi0_k_mat[ik].create(I_aux_num,J_aux_num);                  
    
    // assert(n_kpoints==1);
    for(size_t is=0;is!=mf.get_n_spins();is++)
    {   
        if (tot_ia[is]==0) continue;
        //cout<<" i_spin: "<<is<<endl;
        for(int ik=0;ik!=mf.get_n_kpoints();ik++)
        {
            for(int iq=0;iq!=mf.get_n_kpoints();iq++)
            {
                ComplexMatrix cx_omega_band_mat(construct_omega_band_mat(is,ik,iq,freq));
                if(freq== tfg.get_freq_nodes()[0])
                {
                    print_complex_real_matrix("cx_omega_band_mat",cx_omega_band_mat);
                }
                
                if(ik==iq)
                {
                    ComplexMatrix C_omega=(transpose(C_tilde_map.at(is).at(ik).at(iq).at(Mu),false)) * cx_omega_band_mat;
                    chi0_k_mat.at(ik)+=2*conj(C_omega)*C_tilde_map.at(is).at(iq).at(ik).at(Nu);
                    //chi0_mat[I][J]+=C_omega*reshape_complexmat(occ_num[is],unocc_num[is],atom_mu[J],C_tilde_map.at(is).at(J));//
                   // cout<<"  Cs_step: I J: "<<I<<J<<endl;
                    //chi0_mat[I][J]=C_omega*reshape_complexmat(unocc_num,occ_num,C_tilde_map.at(J).nc,C_tilde_map.at(J));  //careabout second C_tilde_map (mn,mu)  band_index !!!
                          
                }    
                // else{
                //     ComplexMatrix C_omega=(transpose(C_tilde_map.at(is).at(ik).at(iq).at(Mu),false)) * cx_omega_band_mat;
                //     chi0_k_mat.at(ik)+=conj(C_omega)*C_tilde_map.at(is).at(iq).at(ik).at(Nu);
                // }
            }
        }
    }
    return chi0_k_mat;
}


matrix Conventional_Chi0::construct_omega_band_mat(const size_t &is, const size_t ik, const size_t iq, const double &omega)
{
    //cout<<"construct_omega_band_mat"<<endl;
    matrix omega_band_mat(tot_ia[is],tot_ia[is]);
    //complex<double> i_omega=complex<double>(0.0,omega);
    for(int i_occ=0;i_occ!=occ_num[is];i_occ++)
    {
        for(int i_unocc=occ_num[is];i_unocc!=Tot_band_num;i_unocc++)
        {
            double tran_energy=0.5*(mf.get_eigenvals()[is](ik,i_occ)-mf.get_eigenvals()[is](iq,i_unocc));
            //omega_band_mat(unocc_num[is]*i_occ+i_unocc-occ_num[is],unocc_num[is]*i_occ+i_unocc-occ_num[is])= \
                            (2*(wg(is,i_occ)-wg(is,i_unocc))/(tran_energy-i_omega)).real(); 
            omega_band_mat(unocc_num[is]*i_occ+i_unocc-occ_num[is],unocc_num[is]*i_occ+i_unocc-occ_num[is])= \
                            2*(mf.get_weight()[is](ik,i_occ)-mf.get_weight()[is](iq,i_unocc))*tran_energy/(tran_energy*tran_energy+omega*omega);   //!!! consider the spin!
        }
    }
    return omega_band_mat;
}

void Conventional_Chi0::aims_C_tilde(const Cs_LRI &Cs, const vector<Vector3_Order<double>> &qlist)
{
    map<size_t,map<size_t,map<size_t,map<size_t,ComplexMatrix>>>> aims_Cc_tmp=aims_Cc(Cs,qlist);
    cout<<" Begin aims_C_tilde"<<endl;
    int iter_nks=0;
    for(int is=0;is!=mf.get_n_spins();is++)
    {
        for(int ik=0;ik!=mf.get_n_kpoints(); ++ik)
        {
            for(int iq=0;iq!=mf.get_n_kpoints(); ++iq)
            {
                for(int I=0;I!=natom;I++)
                {
                    C_tilde_map[is][ik][iq][I].create((occ_num[is])*unocc_num[is],atom_mu[I]);   
                   // cout<<"   ik iq I  "<<ik<<"    "<<iq<<"   "<<I<<"    dim: "<<C_tilde_map.at(is).at(ik).at(iq).at(I).nr<<"  "<<aims_Ct.at(is).at(ik).at(iq).at(I).nc<<endl;
                    int i_num=atom_nw[I];
                    for(int i_occ=0;i_occ!=occ_num[is];i_occ++)
                    {
                        for(int i_unocc=occ_num[is];i_unocc!=Tot_band_num;i_unocc++)
                        {
                            for(int i_orb=0;i_orb!=i_num;i_orb++)
                            {
                                const size_t iwt1=atom_iw_loc2glo(I,i_orb);
                                for(int mu=0;mu!=atom_mu[I];mu++)
                                {
                                    //cout<<"        "<<i_occ*unocc_num[is]+i_unocc-occ_num[is]<<"  "<<mu<<endl;
                                    
                                    C_tilde_map.at(is).at(ik).at(iq).at(I)(i_occ*unocc_num[is]+i_unocc-occ_num[is],mu)+=mf.get_eigenvectors()[is][ik](i_unocc,iwt1)*aims_Cc_tmp.at(is).at(ik).at(I).at(i_occ)(i_orb,mu) \
                                    + aims_Cc_tmp.at(is).at(iq).at(I).at(i_unocc)(i_orb,mu)*mf.get_eigenvectors()[is][ik](i_occ,iwt1);
                                }
                            }
                        }
                    }
                    
                    print_complex_real_matrix("aims_C_tilde",C_tilde_map.at(is).at(ik).at(iq).at(I));
                }
            }
        }
        iter_nks+=mf.get_n_kpoints();  
    }
}

map<size_t,map<size_t,map<size_t,map<size_t,ComplexMatrix>>>> Conventional_Chi0::aims_Cc(const Cs_LRI &Cs, const vector<Vector3_Order<double>> &qlist)
{
    cout<<"aims_Cc"<<endl;
    map<size_t,map<size_t,map<size_t,map<size_t,ComplexMatrix>>>> Cc_map;
    const auto &LRI_Cs = Cs.data_IJR;     
    int iter_nks=0;
    for(int is=0;is!=mf.get_n_spins();is++)
    {
        for(int ik=0; ik!=mf.get_n_kpoints(); ++ik)
        {
            for(const auto &I_pair:LRI_Cs)
            {
                const auto I=I_pair.first;
                for(int i_state=0;i_state!=mf.get_n_bands();i_state++)
                {
                    Cc_map[is][ik][I][i_state].create(atom_nw[I],atom_mu[I]);
                    cout<<"ik   I  i_state: "<<ik<<"    "<<I<<"    "<<i_state<<endl;
                    for(const auto &J_pair:I_pair.second)
                    {
                        const auto J=J_pair.first;
                        for(const auto &R_pair:J_pair.second)
                        {
                            const auto R=R_pair.first;
                            const auto q=qlist[ik];
                            const double arg = (q * (R * latvec)) * TWO_PI;
                            const complex<double> kphase = complex <double> ( cos(arg),  sin(arg) );
                            const auto &Cs_mat=*R_pair.second;
                            for(int i=0;i!=atom_nw[I];i++)
                                for(int mu=0;mu!=atom_mu[I];mu++)
                                    for(int j=0;j!=atom_nw[J];j++)
                                    {
                                        // if(ik==1)
                                        // {
                                        //     cout<<"in aims_Cc :"<<i<<j<<mu<<"   "<<wfc_k[iter_nks+ik](i_state,atom_iw_loc2glo(J,j))<<Cs_mat(i*atom_nw[J]+j,mu)<<endl;
                                        //     cout<<Cc_map[is][ik][I][i_state](i,mu)<<endl;
                                        // }
                                        Cc_map[is][ik][I][i_state](i,mu)+=mf.get_eigenvectors()[is][ik](i_state,atom_iw_loc2glo(J,j))*Cs_mat(i*atom_nw[J]+j,mu)*kphase;
                                    }
                        }
                    }
                    
                    print_complex_matrix(" aims_Cc ",Cc_map[is][ik][I][i_state]);
                }
            } 
        }
        iter_nks+=mf.get_n_kpoints();   
    }
    return Cc_map;
}
