#include "aperiodic_chi0.h"
#include "Gauss_Quadrature.h"
#include "lapack_connector.h"
#include "global_class.h"
//#include "src_global/matrix.h"
#include "input.h"
#include<fstream>
#include<sstream>
#include<iomanip>
#include<cmath>
#include<cstring>
#include<omp.h>
//#include<iostream>
using namespace std;

void Aperiodic_Chi0::main()
{
    init();
    double cc_tb=omp_get_wtime();
    double cc_te=omp_get_wtime();
    aims_C_tilde();
    double ct_te=omp_get_wtime();
    //map<size_t,map<size_t,ComplexMatrix>> C_tilde_map=C_band_aux(); //C_tilde[is][iat](mn,mu)
    vector<int> part_range=construct_part_range();
    //ph_cRPA(C_tilde_map[0]);
    //map<double,double> freq_grid(read_file_grid("/home/shirong/python_test/freq_grid.txt",'F')); 
    //map<double,double> freq_grid(construct_gauss_grid(20)); 
    map<double,double> freq_grid=cal_chi0.read_local_grid( "local_"+to_string(cal_chi0.grid_N)+"_freq_points.dat",'F');
    first_freq=freq_grid.begin()->first;
    vector<pair<size_t,size_t>> tot_pair(get_atom_pair(Vq));
    complex<double> tot_cRPA(0.0,0.0);
    cout<<"chi0_freq:"<<endl;
    double chi0_tb=omp_get_wtime();
    tmp_count=1;
    for(const auto &freq_pair:freq_grid)
    {
        double freq=freq_pair.first;
        double weight=freq_pair.second;
        map<size_t,map<size_t,map<size_t,ComplexMatrix>>> chi_omega_mat;
        for(size_t atom_pair=0;atom_pair!=tot_pair.size();atom_pair++)
        {
            auto I=tot_pair[atom_pair].first;
            auto J=tot_pair[atom_pair].second;
            map<size_t,ComplexMatrix> chi_omega_mat_tmp(cal_ap_chi0_mat(freq,I,J));
            chi_omega_mat[0][I][J]=std::move(chi_omega_mat_tmp[0]);
        }
        //cout<<"i_freq: "<<freq<<"     "<<chi_omega_mat[0][0](0,0)<<"     "<<chi_omega_mat[0][0](0,1)<<"    "<<chi_omega_mat[0][0](1,1)<<"    "<<chi_omega_mat[0][0](2,2)<<endl;

        map<size_t,map<size_t,ComplexMatrix>> pi_mat=cal_pi_mat(chi_omega_mat[0],freq);
        tmp_count++;
        complex<double> part_cRPA=cal_cRPA_omega(pi_mat,freq,part_range);
        tot_cRPA+=part_cRPA*weight/TWO_PI;
    }
    double chi0_te=omp_get_wtime();
    cout<<"Cc_TIME_USED: "<<cc_te-cc_tb<<endl;
    cout<<"Ct_TIME_USED: "<<ct_te-cc_te<<endl;
    cout<<"Chi0_TIME_USED: "<<chi0_te-chi0_tb<<endl;
    cout<<"ap_cRPA:  "<<setprecision(8)<<tot_cRPA<<endl;
    //cout<<"E_hf+cRPA:"<<tot_cRPA.real()+en.etot<<endl;
}

void Aperiodic_Chi0::init()
{

    cal_chi0.grid_N=6;
    cal_chi0.init();
    cout<<"grid_N: "<<cal_chi0.grid_N<<endl;
    for(size_t i=0;i!=n_kpoints;i++)
    {
        cout<<"i_k:"<<i<<endl;
        //print_complex_matrix( "wfc_k", wfc_k[i]);
        //print_complex_matrix("dm_2D: ", LOC.wfc_dm_2d.dm_k[is]);
    }
    range_all=0;
    cout<<"  atom_mu.size:  "<<atom_mu.size()<<endl;
    for(size_t i_at=0;i_at<atom_mu.size();i_at++)
	{
		range_all+=atom_mu[i_at];
	    cout<<"index_abfs["<<i_at<<"].count_size="<<atom_mu[i_at]<<endl;
	}
    Tot_band_num=NBANDS;
    cout<<"Tot_band_num: "<<Tot_band_num<<endl;
    cout<<" NSPIN= "<<NSPIN<<endl;
    for(int is=0;is!=NSPIN;is++)
    {
        cout<<"tag_0"<<endl;
        int temp_occ=0;
        cout<<"1"<<endl;
        for(int i_band=0;i_band!=Tot_band_num;i_band++)
        {
            cout<<"2"<<endl;
            if(wg(is,i_band)>=1e-6)    //kv.nks == NSPIN ?
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
}

vector<int> Aperiodic_Chi0::construct_part_range()
{
    int range_all;
    vector<int> part_range;
    cout<<" part_range.size(): "<<part_range.size()<<endl;
    part_range.clear();
    cout<<" tag_1"<<endl;
	part_range.push_back(0);
    cout<<" tag_2"<<endl;
	int count_range=0;
	for (int I=0;I!=atom_mu.size()-1;I++)
	{
        cout<<" tag_3"<<endl;
		count_range+=atom_mu[I];
		cout<<" tag_4"<<endl;
        part_range.push_back(count_range);
	}
	
	cout<<"part_range:"<<endl;
    cout<<" part_range.size(): "<<part_range.size()<<endl;
	for(int I=0;I!=atom_mu.size();I++)
	{
		cout<<part_range[I]<<endl;
	}
	cout<<"part_range over"<<endl;
    return part_range;
}

complex<double> Aperiodic_Chi0::cal_cRPA_omega(map<size_t,map<size_t,ComplexMatrix>> &pi_mat, const double &freq,vector<int> &part_range)
{
    //cout<<" cal_cRPA_omega"<<endl;
    ComplexMatrix pi_2D(range_all,range_all);
    //cout<<"range_all: "<<range_all<<endl;
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
	ComplexMatrix identity(range_all,range_all);
	ComplexMatrix identity_minus_pi(range_all,range_all);

	identity.set_as_identity_matrix();

	identity_minus_pi=identity-pi_2D;
    
	complex<double> det_for_rpa(1.0,0.0);
//	det_for_rpa=complexmatrix_det(range_all, identity_minus_pi);
	int info_LU =0;
	int *ipiv=new int[range_all];
	LapackConnector::zgetrf(range_all,range_all,identity_minus_pi,range_all,ipiv,&info_LU);
	for(int ib=0;ib!=range_all;ib++)
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
map<size_t,ComplexMatrix> Aperiodic_Chi0::cal_ap_chi0_mat(const double &omega ,const size_t &I, const size_t &J)
{
    //cout<<"cal_chi0_mat"<<endl;
    const size_t I_aux_num=atom_mu[I];
    const size_t J_aux_num=atom_mu[J];

    map<size_t,ComplexMatrix> chi0_k_mat; 
    for(int ik=0;ik!=n_kpoints/NSPIN;ik++)
        chi0_k_mat[ik].create(I_aux_num,J_aux_num);               
    
    assert(n_kpoints==1);
    for(size_t is=0;is!=NSPIN;is++)
    {   
        if (tot_ia[is]==0) continue;
        //cout<<" i_spin: "<<is<<endl;
        for(int ik=0;ik!=n_kpoints/NSPIN;ik++)
        {
            for(int iq=0;iq!=n_kpoints/NSPIN;iq++)
            {
                ComplexMatrix cx_omega_band_mat(construct_omega_band_mat(is,ik,iq,omega));
                if(omega== first_freq)
                {
                    print_complex_real_matrix("cx_omega_band_mat",cx_omega_band_mat);
                }
                
                if(ik==iq)
                {
                    ComplexMatrix C_omega=(transpose(C_tilde_map.at(is).at(ik).at(iq).at(I),false)) * cx_omega_band_mat;
                    chi0_k_mat.at(ik)+=conj(C_omega)*C_tilde_map.at(is).at(iq).at(ik).at(J);
                    //chi0_mat[I][J]+=C_omega*reshape_complexmat(occ_num[is],unocc_num[is],atom_mu[J],C_tilde_map.at(is).at(J));//
                   // cout<<"  Cs_step: I J: "<<I<<J<<endl;
                    //chi0_mat[I][J]=C_omega*reshape_complexmat(unocc_num,occ_num,C_tilde_map.at(J).nc,C_tilde_map.at(J));  //careabout second C_tilde_map (mn,mu)  band_index !!!
                          
                }
            }
        }
    }
    return chi0_k_mat;
}

map<size_t,map<size_t,ComplexMatrix>> Aperiodic_Chi0::cal_pi_mat(map<size_t,map<size_t,ComplexMatrix>> &chi0_omega_mat, const double &freq)
{
    //cout<<"cal_pi_mat"<<endl;
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
                        if(freq==first_freq)
                        {
                            cout<<" IJQ: "<<I<<" "<<J<<" "<<Q<<endl;
                            print_complex_real_matrix("vq:",Vq_mat);
                            print_complex_real_matrix("chi0:",chi0_omega_mat.at(J).at(Q));
                            print_complex_real_matrix("pi_mat:",pi_mat.at(I).at(Q));
                        }
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
    // cout<<" cap_pi_1"<<endl;
    // for(const auto &Q_pair:chi0_omega_mat)
    // {
    //     const auto Q=Q_pair.first;
    //     //const auto adj_Q=Abfs::get_adjs(Q);
    //     for(const auto &I_pair:Q_pair.second)
    //     {
    //         const auto I=I_pair.first;
    //         ComplexMatrix &chi0_mat=chi0_omega_mat.at(I).at(Q);
    //         for(const auto &J_pair:Q_pair.second)
    //         {
    //             const auto J=J_pair.first;
    //             cout<<"  cal_pi  IJQ:  "<<I<<J<<Q<<endl;
    //             ComplexMatrix temp_Vps(*Vq.at(Q).at(J).at({0,0,0}));
    //             pi_mat[I][J]+=chi0_mat*temp_Vps;
              
    //             // if(tmp_count==1)
    //             // {
    //             //     cout<<" cal_pi_mat_step: I J Q: "<<I<<J<<Q<<endl;
    //             //     print_complex_real_matrix("chi0_mat:",chi0_mat);
    //             //     print_complex_real_matrix("Vps_mat:",temp_Vps);
    //             //     print_complex_real_matrix("pi_tot",pi_mat[I][J]);
    //             // }  
                
    //         }
    //     }
    // }
    return pi_mat;
}

map<size_t,map<size_t,ComplexMatrix>> Aperiodic_Chi0::C_band_aux()
{
    cout<<"C_band_aux"<<endl;
    map<size_t,map<size_t,ComplexMatrix>> C_tilde_map;
    for(int is=0;is!=NSPIN;is++)
    {
        for(const auto &I_pair:Cs)
        {
            const auto I=I_pair.first;
            C_tilde_map[is][I].create(tot_ia[is],atom_mu[I]);
            for(const auto &J_pair:I_pair.second)
            {
                const auto J=J_pair.first;
                cout<<"  i_spin: "<<is<<"    C_band_aux:   I J: "<<I<<J<<endl;
                ComplexMatrix cc_plus_cc=construct_cc_plus_cc(is,I,J);
                print_complex_real_matrix( "cc_plus_cc ", cc_plus_cc );  
                for(const auto &R_pair:J_pair.second)
                {
                    const auto &Cs_mat=R_pair.second;
                    ComplexMatrix Cs_cx(*Cs_mat);
                    assert(Cs_cx.nr==cc_plus_cc.nc);
                    //print_complex_real_matrix( "Cs_mat ", Cs_cx );
                    ComplexMatrix temp_C_tilde=cc_plus_cc*Cs_cx;
                    //print_complex_real_matrix("temp_C_tilde",temp_C_tilde); 
                    C_tilde_map[is][I]+=cc_plus_cc*Cs_cx;
                }
            }
        }    
    }

    return C_tilde_map;
}

map<size_t,map<size_t,map<size_t,map<size_t,ComplexMatrix>>>> Aperiodic_Chi0::aims_Cc()
{
    cout<<"aims_Cc"<<endl;
    map<size_t,map<size_t,map<size_t,map<size_t,ComplexMatrix>>>> Cc_map;
        
    int iter_nks=0;
    for(int is=0;is!=NSPIN;is++)
    {
        for(int ik=0; ik!=n_kpoints/NSPIN; ++ik)
        {
            for(const auto &I_pair:Cs)
            {
                const auto I=I_pair.first;
                for(int i_state=0;i_state!=NBANDS;i_state++)
                {
                    Cc_map[is][ik][I][i_state].create(atom_nw[I],atom_mu[I]);
                    cout<<"ik   I  i_state: "<<ik<<"    "<<I<<"    "<<i_state<<endl;
                    for(const auto &J_pair:I_pair.second)
                    {
                        const auto J=J_pair.first;
                        for(const auto &R_pair:J_pair.second)
                        {
                            const auto R=R_pair.first;
                            const double arg = ( kvec_c[ik]* (R*latvec)) * TWO_PI;
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
                                        Cc_map[is][ik][I][i_state](i,mu)+=wfc_k[iter_nks+ik](i_state,atom_iw_loc2glo(J,j))*Cs_mat(i*atom_nw[J]+j,mu)*kphase;
                                    }
                        }
                    }
                    
                    print_complex_matrix(" aims_Cc ",Cc_map[is][ik][I][i_state]);
                }
            } 
        }
        iter_nks+=n_kpoints/NSPIN;   
    }
    return Cc_map;
}

void Aperiodic_Chi0::aims_C_tilde()
{
    map<size_t,map<size_t,map<size_t,map<size_t,ComplexMatrix>>>> aims_Cc_tmp=aims_Cc();
    cout<<" Begin aims_C_tilde"<<endl;
    int iter_nks=0;
    for(int is=0;is!=NSPIN;is++)
    {
        for(int ik=0;ik!=n_kpoints/NSPIN; ++ik)
        {
            for(int iq=0;iq!=n_kpoints/NSPIN; ++iq)
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
                                    
                                    C_tilde_map.at(is).at(ik).at(iq).at(I)(i_occ*unocc_num[is]+i_unocc-occ_num[is],mu)+=wfc_k[iter_nks+iq](i_unocc,iwt1)*aims_Cc_tmp.at(is).at(ik).at(I).at(i_occ)(i_orb,mu) \
                                                                                    + aims_Cc_tmp.at(is).at(iq).at(I).at(i_unocc)(i_orb,mu)*wfc_k[iter_nks+ik](i_occ,iwt1);
                                    // if(i_occ==0 && i_unocc==occ_num[is] && mu==3)
                                    // {
                                    //     cout<<"  C_tilde_map sum rounting: "<<i_occ<<"   "<<i_unocc<<endl;
                                    //     cout<<iwt1<<"   "<<fixed<<setprecision(8)<<wfc_k[iter_nks+iq](i_unocc,iwt1).real()<<"   "<<aims_Cc_tmp.at(is).at(ik).at(I).at(i_occ)(i_orb,mu).real() <<"   "<<aims_Cc_tmp.at(is).at(iq).at(I).at(i_unocc)(iwt1,mu).real()<<"  "<<wfc_k[iter_nks+ik](i_occ,iwt1).real()<<endl;
                                    //     cout<<C_tilde_map.at(is).at(ik).at(iq).at(I)(i_occ*unocc_num[is]+i_unocc-occ_num[is],mu)<<endl;
                                    // }
                                }
                            }
                        }
                    }
                    
                    print_complex_real_matrix("aims_C_tilde",C_tilde_map.at(is).at(ik).at(iq).at(I));
                }
            }
        }
        iter_nks+=n_kpoints/NSPIN;  
    }
}

ComplexMatrix Aperiodic_Chi0::construct_cc_plus_cc(const size_t &is, const size_t &I, const size_t &J)
{
    cout<<"construct_cc_plus_cc"<<endl;
  /*   const size_t it1=ucell.iat2it[I];
	const size_t ia1=ucell.iat2ia[I];
    
    const size_t it2=ucell.iat2it[J];
	const size_t ia2=ucell.iat2ia[J]; */
    
    
    int i_num=atom_nw[I];
    int j_num=atom_nw[J];
    ComplexMatrix cc_plus_cc(tot_ia[is], i_num*j_num);
    for(int i_occ=0;i_occ!=occ_num[is];i_occ++)
    {
        for(int i_unocc=occ_num[is];i_unocc!=Tot_band_num;i_unocc++)
        {
            for(int i_orb=0;i_orb!=i_num;i_orb++)
            {
                const size_t iwt1=atom_iw_loc2glo(I,i_orb);
                for(int j_orb=0;j_orb!=j_num;j_orb++)
                {
                    const size_t iwt2=atom_iw_loc2glo(J,j_orb);
                    cc_plus_cc(i_occ*unocc_num[is]+i_unocc-occ_num[is],i_orb*j_num+j_orb)=conj(wfc_k[is](i_occ,iwt1))*wfc_k[is](i_unocc,iwt2) \
                                                                                    +conj(wfc_k[is](i_occ,iwt2))*wfc_k[is](i_unocc,iwt1);
                }
            }
            
        }
    }
    return cc_plus_cc;
}     

matrix Aperiodic_Chi0::construct_omega_band_mat(const size_t &is, const size_t ik, const size_t iq, const double &omega)
{
    //cout<<"construct_omega_band_mat"<<endl;
    matrix omega_band_mat(tot_ia[is],tot_ia[is]);
    //complex<double> i_omega=complex<double>(0.0,omega);
    for(int i_occ=0;i_occ!=occ_num[is];i_occ++)
    {
        for(int i_unocc=occ_num[is];i_unocc!=Tot_band_num;i_unocc++)
        {
            double tran_energy=0.5*(ekb[ik][i_occ]-ekb[iq][i_unocc]);
            //omega_band_mat(unocc_num[is]*i_occ+i_unocc-occ_num[is],unocc_num[is]*i_occ+i_unocc-occ_num[is])= \
                            (2*(wg(is,i_occ)-wg(is,i_unocc))/(tran_energy-i_omega)).real(); 
            omega_band_mat(unocc_num[is]*i_occ+i_unocc-occ_num[is],unocc_num[is]*i_occ+i_unocc-occ_num[is])= \
                            2*(wg(ik,i_occ)-wg(iq,i_unocc))*tran_energy/(tran_energy*tran_energy+omega*omega);   //!!! consider the spin!
        }
    }
    return omega_band_mat;
}

ComplexMatrix Aperiodic_Chi0::reshape_complexmat(const size_t n1, const size_t n2, const size_t n3, const ComplexMatrix &mat) //(n1,n2,n3) -> (n2,n1,n3)
{
    const auto length = sizeof(complex<double>)*n3;
	const auto n13 = n1*n3;
	const complex<double> *m_ptr = mat.c;
	ComplexMatrix m_new( n1*n2, n3, false );
	for( size_t i1=0; i1!=n1; ++i1 )
	{
		complex<double> *m_new_ptr = m_new.c+i1*n3;
		for( size_t i2=0; i2!=n2; ++i2, m_ptr+=n3, m_new_ptr+=n13 )
			memcpy( m_new_ptr, m_ptr, length );
	}
    return m_new;
}
void Aperiodic_Chi0::ph_cRPA(const map<size_t,map<size_t,ComplexMatrix>> &C_tilde_map)
{
    vector<ComplexMatrix> ph_CVC(NSPIN);
    for(size_t is=0;is!=NSPIN;is++)
    {
        if (tot_ia[is]==0) continue;
        ph_CVC[is].create(tot_ia[is],tot_ia[is]);
        ComplexMatrix delta_occ(construct_delta_occ_mat(is));
        
        for(const auto &Q_pair:C_tilde_map.at(is))
        {
            const auto Q=Q_pair.first;
            ComplexMatrix ph_CV(tot_ia[is],atom_mu[Q]);
            for(const auto &I_pair:C_tilde_map.at(is))
            {
                const auto I=I_pair.first;
                const ComplexMatrix complex_Vps((*Vq.at(I).at(Q).at({0,0,0})));
                ph_CV+=C_tilde_map.at(is).at(I)*complex_Vps;
            }
            ph_CVC[is]+=2*ph_CV*transpose(C_tilde_map.at(is).at(Q),false);
        }
        ph_CVC[is]=ph_CVC[is]*delta_occ;
        print_complex_matrix("ph_CVC",ph_CVC[is]);
    }
    
    vector<ComplexMatrix> ph_CVC_plus_band(NSPIN);
    for(size_t is=0;is!=NSPIN;is++)
    {
        if (tot_ia[is]==0) continue;
        ph_CVC_plus_band[is].create(tot_ia[is],tot_ia[is]);
        ComplexMatrix band_mat(construct_band_mat(is));
        ph_CVC_plus_band[is]=ph_CVC[is]+band_mat;
        print_complex_matrix("ph_CVC_plus_band",ph_CVC_plus_band[is]);
    }
    
    int max_occ_num=0;
    int min_occ_num=occ_num[0];
    int max_unocc_num=0;
    for(size_t is=0;is!=NSPIN;is++)
    {
        if(occ_num[is]>max_occ_num)
            max_occ_num=occ_num[is];
        if(unocc_num[is]>max_unocc_num)
            max_unocc_num=unocc_num[is];
        if(occ_num[is]<min_occ_num)
            min_occ_num=occ_num[is];
    }
    int max_ia=max_occ_num*max_unocc_num;
    cout<<"  max_occ , min_occ, max_unocc num, : "<<max_occ_num<<"    "<<min_occ_num<<"   "<<max_unocc_num<<endl;
    
    ComplexMatrix SP_B(max_ia,max_ia);
    ComplexMatrix SP_A(max_ia,max_ia);
    for(size_t is=0;is!=NSPIN;is++)
    {
        if (tot_ia[is]==0) continue;
        for(int i_occ=0;i_occ!=occ_num[is];i_occ++)
        {
            for(int i_unocc=occ_num[is];i_unocc!=Tot_band_num;i_unocc++)
            {
                for(int j_occ=0;j_occ!=occ_num[is];j_occ++)
                {
                    for(int j_unocc=occ_num[is];j_unocc!=Tot_band_num;j_unocc++)
                    {
                        cout<<" is  max_ia,max_jb   tot_ia,tot_jb"<<is<<"   "<<max_unocc_num*i_occ+i_unocc-min_occ_num<<max_unocc_num*j_occ+j_unocc-min_occ_num<< \
                                "     "<<unocc_num[is]*i_occ+i_unocc-occ_num[is]<<unocc_num[is]*j_occ+j_unocc-occ_num[is]<<endl;
                        SP_B(max_unocc_num*i_occ+i_unocc-min_occ_num,max_unocc_num*j_occ+j_unocc-min_occ_num)+= \
                                ph_CVC[is](unocc_num[is]*i_occ+i_unocc-occ_num[is],unocc_num[is]*j_occ+j_unocc-occ_num[is]);
                                
                        SP_A(max_unocc_num*i_occ+i_unocc-min_occ_num,max_unocc_num*j_occ+j_unocc-min_occ_num)+= \
                                ph_CVC_plus_band[is](unocc_num[is]*i_occ+i_unocc-occ_num[is],unocc_num[is]*j_occ+j_unocc-occ_num[is]);
                    }
                }
            }
        }
    }
    
    print_complex_matrix("SP_B",SP_B);
    print_complex_matrix("SP_A",SP_A);
    
    matrix SP_mat(2*max_ia,2*max_ia);
    for(int i=0;i!=max_ia;i++)
        for(int j=0;j!=max_ia;j++)
            SP_mat(i,j)=SP_A(i,j).real();
    for(int i=0;i!=max_ia;i++)
        for(int j=0;j!=max_ia;j++)
            SP_mat(max_ia+i,j)=-1.0*SP_B(i,j).real();
    for(int i=0;i!=max_ia;i++)
        for(int j=0;j!=max_ia;j++)
            SP_mat(i,max_ia+j)=SP_B(i,j).real();
    for(int i=0;i!=max_ia;i++)
        for(int j=0;j!=max_ia;j++)
            SP_mat(max_ia+i,max_ia+j)=-1.0*SP_A(i,j).real();   
    print_matrix("  SP_mat ",SP_mat);
    complex<double> tr_A=trace(SP_A);
    int info;	
	vector<complex<double>> eigen_value;
    matrix r_eg_vec(SP_mat.nr,SP_mat.nc);
    dgeev(SP_mat,eigen_value,r_eg_vec,info);
	//LapackConnector::dsyev('V','L',SP_mat,VECTOR_TO_PTR(eigen_value),info);
    print_matrix("  SP_mat after diagonal ",SP_mat);
    double omega_ph=0.0;
    cout<<" SP_mat eigen_value:"<<endl;
    for(int i=0;i!=eigen_value.size();i++)
    {
        cout<<"    "<<eigen_value[i];
        if(eigen_value[i].real()>0)
        {
            omega_ph+=eigen_value[i].real();
        }
    }
    print_matrix("r_eg_vec",r_eg_vec);
    cout<<endl;
    cout<<"  omega_ph "<<omega_ph<<"     tr_A:"<<tr_A<<endl;
    double ph_cRPA=0.5*(omega_ph-tr_A.real());
    cout<<"  ph_cRPA:  "<<setprecision(8)<<ph_cRPA<<endl;
}
matrix Aperiodic_Chi0::construct_band_mat(const size_t &is)
{
    //cout<<"construct_omega_band_mat"<<endl;
    matrix band_mat(tot_ia[is],tot_ia[is]);
    
    for(int i_occ=0;i_occ!=occ_num[is];i_occ++)
    {
        for(int i_unocc=occ_num[is];i_unocc!=Tot_band_num;i_unocc++)
        {
            double tran_energy=0.5*(ekb[is][i_unocc]-ekb[is][i_occ]);
            band_mat(unocc_num[is]*i_occ+i_unocc-occ_num[is],unocc_num[is]*i_occ+i_unocc-occ_num[is])= tran_energy;
        }
    }
    return band_mat;
}
matrix Aperiodic_Chi0::construct_delta_occ_mat(const size_t &is)
{
    matrix delta_occ_mat(tot_ia[is],tot_ia[is]);
    
    for(int i_occ=0;i_occ!=occ_num[is];i_occ++)
    {
        for(int i_unocc=occ_num[is];i_unocc!=Tot_band_num;i_unocc++)
        {
            double delta_occ=wg(is,i_occ)-wg(is,i_unocc);
            delta_occ_mat(unocc_num[is]*i_occ+i_unocc-occ_num[is],unocc_num[is]*i_occ+i_unocc-occ_num[is])= 0.5*NSPIN*delta_occ;
        }
    }
    return delta_occ_mat;
}
map<double,double> Aperiodic_Chi0::construct_gauss_grid(const int Npoints)
{
    map<double,double> grid;
    grid.clear();
    cout<<"grid.size_1:"<<grid.size()<<endl;
    Gauss_Quadrature gauss_FT( 'L', Npoints);
    gauss_FT.open_Interval( 0.0, 'R' );	
    cout<<"finish gauss"<<endl;
    for(int i=0;i!=Npoints;++i)
    {
        //std::cout<<"gauss_FT.nodes:"<<dec<<gauss_FT.nodes[i]<<"    "<<dec<<gauss_FT.weights[i]<<std::endl; 
        grid.insert({gauss_FT.nodes[i],gauss_FT.weights[i]});
    }
    cout<<"grid.size_2:"<<grid.size()<<endl;
    cout<<"grid insert"<<endl;
    return grid;
}

map<double, double> Aperiodic_Chi0::get_minimax_grid(const size_t &npoints)
{
    vector<double> omega,womega;
    switch(npoints)
    {
        case 10:{
                omega={ 0.3481615248557796  ,
1.2313725909827136  ,
2.810913357339209  , 
6.15080805896027  ,  
13.877541481835147  ,
33.22069104542858  , 
86.44672242760085  , 
253.20259292285937  ,
891.1843910704623  , 
4451.847505887757 };
                womega={     0.7255994080655348  ,
   1.119228495374018  ,
  2.2036033180789114  ,
 4.872326808305027  ,
   11.650266939579032  ,
  30.235333999867677  ,
  87.14811713700114  ,
   291.510370299167  ,
  1237.6763436699034  ,
  8510.868573356214  }; break;
                }
        case 24:{
            omega={0.26698215681967413,
                   0.892223906310287  ,
                   1.8317378189492126 ,
                   3.448540624873043  ,
                   6.412332341704699  ,
                   12.004764760267175 ,
                   22.59296241728836  ,
                   41.726076718742085 ,
                   72.58076851129877  ,
                   121.78239886144758 ,
                   196.90750359477803 ,
                   308.9975125445375  ,
                   388.8935803445273  ,
                   688.4707600323343  ,
                   920.8807873273335  ,
                   2015.166433612731  ,
                   2356.565126502859  ,
                   5523.102564973685  ,
                   7655.084845919751  ,
                   14613.87657199115  ,
                   28273.763715937504 ,
                   62195.263439263705 ,
                   134547.49387371115 ,
            187295.53091852312 };
            womega={   0.1746132987395890,
                     0.23476218050354886 ,
                      0.3817729398008738 ,
                     0.6814853128853732  ,
                     1.2706200979578386  ,
                      2.41288518115062   ,
                     4.527058405056739   ,
                      7.815575677663938  ,
                     12.129166857287595  ,
                      19.678750472569515 ,
                      29.443140410052337 ,
                     27.635501922848125  ,
                     55.73445271426729   ,
                     84.65794107126393   ,
                     157.38596798021848  ,
                     469.9087693972736   ,
                     169.76980051838893  ,
                     1527.8095913852105  ,
                     548.9780030592839   ,
                     2779.8412053235775  ,
                      8705.670002053863  ,
                      6978.774383166463  ,
                      34967.87966785437  ,
            69917.60629628583 }; break;   
            }
/*            
                omega={       0.216947863738504304,
       0.696143263930536826,
       1.32393538227083329 ,
       2.24493289065218526 ,
       3.6857244350672298  ,
       6.0224506584234998  ,
       9.9008803388491895  ,
      16.4565791110765716  ,
      27.7282619261081713  ,
      47.4519548172217114  ,
      82.6245046141143717  ,
     146.664980422799914   ,
     265.999541189142519   ,
     494.243126716021436   ,
     943.9272183571494     ,
    1860.66714144249931    ,
    3805.67254521028053    ,
    8133.19404068376571    ,
   18336.1169518854222     ,
   44209.6699195113397     ,
  116406.949619594176      ,
  346689.530268761911      ,
 1250711.43098220462       ,
 6544186.61639788095};
                womega={      0.140451214016967085,
      0.169948669567292127,
      0.237182977881646939,
      0.360822736407188649,
      0.575764869953818259,
      0.944997311385078032,
      1.58226868303556811 ,
      2.69543067725147889 ,
      4.67081982406709617 ,
      8.24133262271182154 ,
     14.8299750745752927  ,
     27.2727544002122144  ,
     51.3902319089963342  ,
     99.5327191150612123  ,
    198.923242962995772   ,
    412.27689279712223    ,
    891.781785149740472   ,
   2030.51168298919492    ,
   4924.61235045882768    ,
  12943.4348138283331     ,
  37872.8556917399255     ,
 129237.556631576008      ,
 565357.867716293666      ,
4147374.97977405787}; break; 

        } 
      */  
    }
    double gap;
    gap=get_band_gap();
    map<double,double> minimax_grid;
    minimax_grid.clear();
    for(size_t i=0;i!=omega.size();++i)
        minimax_grid.insert({gap*omega[i],PI*gap*womega[i]});
    cout<<" MAXMINI_GRID "<<endl;
    for(const auto &m:minimax_grid)
        cout<<m.first<<"      "<<m.second<<endl;
    return minimax_grid;
} 

map<double,double> Aperiodic_Chi0::read_file_grid(const string &file_path, const char type)
{
    map<double,double> grid;
    ifstream infile;
    //string file_path="/home/shirong/python_test/freq_grid.txt";
    infile.open(file_path);
    
    vector<double> tran;
    string s0;
    int n_freq=0;
    //stringstream ss;
    //double ss_d;
    string pattern="          ";
    while(getline(infile,s0))
    {
        std::string::size_type pos;
        std::vector<std::string> result;
        s0+=pattern;//扩展字符串以方便操
    
        for(int i=0; i!=s0.size(); i++)
        {
            pos=s0.find(pattern,i);
            if(pos<s0.size())
            {
                std::string s=s0.substr(i,pos-i);
                result.push_back(s);
                i=pos+pattern.size()-1;
            }
        }
        grid.insert(pair<double,double>(std::stod(result[0]),std::stod(result[1])));
    }    
    infile.close();
    
    double gap;
    gap=get_band_gap();
    map<double,double> minimax_grid;
    minimax_grid.clear();
    switch(type)
    {
        case 'F':
        {
            for(auto &i_pair:grid)
                minimax_grid.insert({i_pair.first*gap,i_pair.second*gap*PI});
            
            cout<<" MINIMAX_GRID_Freq "<<endl;
            for(const auto &m:minimax_grid)
                cout<<m.first<<"      "<<m.second<<endl;
            break;
        }
        
        case 'T':
        {
            for(auto i_pair:grid)
                minimax_grid.insert({i_pair.first/(gap),i_pair.second/(gap)});
            
            cout<<" MINIMAX_GRID_Tau "<<endl;
            for(const auto &m:minimax_grid)
                cout<<m.first<<"      "<<m.second<<endl;
            break;
        }
    }
    
    return minimax_grid;
}
double Aperiodic_Chi0::get_band_gap()
{
    int HOMO_level;
	HOMO_level=0;
	for(int n=0;n!=NLOCAL;n++)
	{
		if(wg(0,n)>=1e-6)
		{
			HOMO_level+=1;
		}
		
	}
    double gap;
    gap=0.5*(ekb[0][HOMO_level]-ekb[0][HOMO_level-1]);
    
	cout<<"HOMO_level="<<HOMO_level<<"     gap:"<<gap<<endl;
    return gap;
}

void Aperiodic_Chi0::print_matrix( char* desc, const matrix &mat )
{
    int nr=mat.nr;
    int nc=mat.nc;
    printf( "\n %s\n", desc );
    for( int i = 0; i < nr; i++ ) 
    {
        for( int j = 0; j < nc; j++ ) printf( " %12.8f", mat.c[i*nc+j] );
            printf( "\n" );
    }
} 
void Aperiodic_Chi0::print_complex_matrix( char* desc, const ComplexMatrix &mat )
{
    int nr=mat.nr;
    int nc=mat.nc;
    printf( "\n %s\n", desc );
    for( int i = 0; i < nr; i++ ) 
    {
        for( int j = 0; j < nc; j++ ) printf("%8.4f ,%6.4f ",mat.c[i*nc+j].real(),mat.c[i*nc+j].imag());
            printf( "\n" );
    }
} 

void Aperiodic_Chi0::print_complex_real_matrix( char* desc, const ComplexMatrix &mat )
{
    int nr=mat.nr;
    int nc=mat.nc;
    printf( "\n %s\n", desc );
    for( int i = 0; i < nr; i++ ) 
    {
        for( int j = 0; j < nc; j++ ) printf( " %12.8f", mat.c[i*nc+j].real());
            printf( "\n" );
    }
} 

vector<pair<size_t,size_t>> Aperiodic_Chi0::get_atom_pair(const map<size_t,map<size_t,map<Vector3_Order<double>,std::shared_ptr<ComplexMatrix>>>> &m)
{
    vector<pair<size_t,size_t>> atom_pairs;
    for(const auto &mA : m)
		for(const auto &mB : mA.second)
			atom_pairs.push_back({mA.first,mB.first});
	printf("  tot atom_pairs tasks: %d\n",atom_pairs.size());
    return atom_pairs;
}

extern "C"
{
    void dgeev_(const char* jobvl,const char* jobvr,const int* n,double *a,
                const int* lda,double* wr,double* wi,double* vl, const int* ldvl, 
                double* vr, const int* ldvr,  double* work,const int* lwork, int* info);
}

void Aperiodic_Chi0::dgeev(matrix &mat, vector<complex<double>> &egValue, matrix &vr, int info)
{
    assert(mat.nr==mat.nc);
    vector<double> eg_r(mat.nr);
    vector<double> eg_i(mat.nr);
    matrix vl(mat.nr,mat.nc);
    vr.create(mat.nr,mat.nc);
    int lwork=16*mat.nr;
    vector<double> work(std::max(1,lwork));
    dgeev_("N","V",&mat.nr,mat.c,&mat.nr,VECTOR_TO_PTR(eg_r),VECTOR_TO_PTR(eg_i),vl.c,&mat.nr,vr.c,&mat.nr,VECTOR_TO_PTR(work),&lwork,&info);
    for(int i=0;i!=mat.nr;i++)
        egValue.push_back(complex<double>(eg_r[i],eg_i[i]));
    vr=transpose(vr);
}
