#include "cal_periodic_chi0.h"
#include "lapack_connector.h"
#include <omp.h>
#include <iomanip>
#include "global_class.h"
#include <iostream>
/*
#include "scalapack_connector.h"
#include "src_lcao/global_fp.h"
#include "src_pw/global.h"
#include "src_lcao/exx_lcao.h"
#include "src_lcao/abfs.h"
#include "src_external/src_test/test_function.h"
#include "Gauss_Quadrature.h"
*/

using namespace std;

void Cal_Periodic_Chi0::chi0_main() 
{
    init();
    double t_begin=omp_get_wtime();
    
    int p_id_local=para_mpi.get_myid();
    //LOC.wfc_dm_2d.cal_dm();
    cout<<"finish init"<<endl;
	set<Vector3_Order<int>> R_grid(construct_R_grid());
	cout<<"R_grid_num= "<<R_grid.size()<<endl; 
	cout<< minimax_grid_path+"/time_grid.txt"<<endl;
    map<double,double> time_grid(read_file_grid( minimax_grid_path+"/time_grid.txt",'T'));    
     
    map<double,double> freq_grid(read_file_grid( minimax_grid_path+"/freq_grid.txt",'F'));

    vector<double> tau_vec;
    for(auto tp:time_grid)
        tau_vec.push_back(tp.first);
    //copy(freq_tmp.begin(),freq_tmp.end(),inserter(freq_grid,freq_grid.begin()));
    //cout<<"freq_grid size: "<<freq_grid.size();
    auto iter_tau=time_grid.begin();
    first_tau=iter_tau->first;
    auto iter_freq=freq_grid.begin();
    //iter_freq++; 
    first_freq=iter_freq->first;
    //first_freq=freq_grid.begin()->first;
    cout<<"time_grid.size() "<<time_grid.size()<<endl;
    //construct_gauss_grid(20);
    cout<<"finish time_grid init"<<endl;

    
	for(const auto &time_tau:time_grid)
    {
        //cout<<"time_tau:"<<time_tau.first<<endl;
		for(const Vector3_Order<int> &R:R_grid)
        {
            cal_Green_func_R_tau( time_tau.first, R, kvec_c);
            cal_Green_func_R_tau( -1*time_tau.first, R, kvec_c);
        }
	}
    
    cout<<" FINISH GREEN FUN"<<endl;
    // for(auto &R:R_grid)
    // {
    //     //cout<<R<<endl;
    //     for(int I=0;I!=natom;I++)
    //     {
    //         for (int J=0;J!=natom;J++)
    //         {
    //             cout<<"   I  J  R   "<<I<<"  "<<J<<"  "<<R<<endl;
    //            // print_matrix(" Green_atom", Green_atom[0][I][J][time_grid.begin()->first][R]);
    //            // print_matrix(" Green_atom  -tau", Green_atom[0][I][J][-time_grid.begin()->first][R]);
    //         }
    //     }
    // }
    cout<<"  time_grid_size   "<<time_grid.size()<<endl;
    //vector<vector<pair<double,Vector3_Order<int>>>> p_task(process_task_tau_R(time_grid,R_grid));
    //vector<pair<size_t,size_t>> tot_pair(get_atom_pair(Vq));
    vector<pair<size_t,size_t>> tot_pair(process_task_ap_local());
    cout<<"  tot_pair.size :  "<<tot_pair.size()<<endl;
    double t_chi0_begin=omp_get_wtime();
    omp_lock_t chi0_lock;
    omp_init_lock(&chi0_lock);
    omp_lock_t t_lock;
    omp_init_lock(&t_lock);
    double time_task_tot=0;
   // print_matrix("Cs_last",*Cs.at(natom-1).at(natom-1).at(Vector3_Order<int>(0,0,0)));
    // for(auto &R:R_grid)
    // {
    //     cout<<"R:  "<<R<<endl;
    //     print_matrix("Green_last",Green_atom.at(0).at(natom-1).at(natom-1).at(first_tau).at(R));
    // }
    //matrix tmp_chi0(cal_chi0_element(first_tau,Vector3_Order<int>(0,0,0),natom-1,natom-1));

    #pragma omp parallel for 
   
    // for(size_t tau_R_index=0;tau_R_index!=p_task[p_id_local].size();tau_R_index++)
    // {
    //     auto tau=p_task[p_id_local][tau_R_index].first;
    //     auto R=p_task[p_id_local][tau_R_index].second;
    //     //cout<<" CHI0  thread, tau, R:   "<<omp_get_thread_num()<<"   "<<tau<<R<<endl;
    //     printf("CHI0 p_id:   %d,   thread: %d,  tau: %f,   R: ( %d,  %d, %d )\n ",p_id_local,omp_get_thread_num(),tau,R.x,R.y,R.z);
    //     for(auto &I_p:Vq)
    //     {
    //         const size_t I=I_p.first;
    //         for(auto &J_p:I_p.second)
    //         {
    //             const size_t J=J_p.first;
    //             cal_chi0_element(tau,R,I,J);
    //         }
    //     }
    // }
    for(size_t atom_pair=0;atom_pair!=tot_pair.size();atom_pair++)
    {
        auto I=tot_pair[atom_pair].first;
        auto J=tot_pair[atom_pair].second;
        //cout<<" CHI0  thread, tau, R:   "<<omp_get_thread_num()<<"   "<<tau<<R<<endl;
       // printf("CHI0 p_id:   %d,   thread: %d,  I: %d,   J: %d\n ",p_id_local,omp_get_thread_num(),I,J);
        //map<double,double> time_grid_local(read_file_grid( minimax_grid_path+"/time_grid.txt",'T'));    
        double task_begin=omp_get_wtime();

        for(auto it=0;it!=tau_vec.size();it++)
        {
            for(auto &R:R_grid)
            {
               // printf("   out thread: %d, IJ:   %d,%d,  R:(  %d,%d,%d  )  tau:%f\n",omp_get_thread_num(),I,J,R.x,R.y,R.z,tau_vec[it]);
                matrix tmp_chi0(cal_chi0_element(tau_vec[it],R,I,J));
                omp_set_lock(&chi0_lock);
                chi0[tau_vec[it]][R][I][J]=std::move(tmp_chi0);
                omp_unset_lock(&chi0_lock);
            }
        }
        double task_end=omp_get_wtime();
        double time_used=task_end-task_begin;
        omp_set_lock(&t_lock);
        time_task_tot+=time_used;
        printf("CHI0 p_id: %d,   thread: %d,     I: %d,   J: %d,      TIME_USED: %f \n ",p_id_local,omp_get_thread_num(),I,J,time_used);
        omp_unset_lock(&t_lock);
    }
    omp_destroy_lock(&chi0_lock);
    omp_destroy_lock(&t_lock);
   // int mid_atom=(natom-1)/2;
    //cout<<"  !!!  mid_atom="<<mid_atom;
    #pragma omp barrier
    double t_chi0_end=omp_get_wtime();
    #pragma omp single
    for(const auto &tau:time_grid)
    {
        cout<<" tau  "<<tau.first<<"   chi0";
        rt_m_max(chi0[tau.first][{0,0,0}][natom-1][natom-1]);
    }
 
    Cosine_to_chi0_freq(time_grid,freq_grid);
    for(const auto &freq:freq_grid)
    {
        cout<<"freq  "<<freq.first<<"   chi0_freq";
        rt_m_max(chi0_freq[freq.first][{0,0,0}][natom-1][natom-1]);
    }
    //cout<<"finish FT_to_chi0_freq"<<endl;
    

    // for(auto &R:R_grid)
    // {
    //     cout<<R<<endl;
    //     print_matrix("chi0_freq_mat  real",chi0_freq[freq_grid.begin()->first][R][0][0]);
    // }

    FT_R_to_K_chi0();
    for(const auto &freq:freq_grid)
    {
        cout<<"freq   "<<freq.first<<"   chi0_k";
        rt_cm_max(chi0_k[freq.first][{0,0,0}][natom-1][natom-1]);
    }
    // for(int i=0;i!=n_kpoints;i++)
    // {
    //     cout<<kvec_c[i];
    //     print_complex_matrix("chi0_k",chi0_k[freq_grid.begin()->first][kvec_c[i]][0][0]);
    // }
    //print_complex_matrix("last vq",(*Vq.at(natom-1).at(natom-1).at(Vector3_Order<double>{0,0,0})));
    cal_pi_k_use_aims_vq();
 
    for(auto &freq_pair:pi_k)
        cout<<dec<<freq_pair.first<<"     "<<pi_k[freq_pair.first][{0,0,0}][0][0](0,0)<<"     "<<pi_k[freq_pair.first][{0,0,0}][0][0](0,1)<<"     "<<pi_k[freq_pair.first][{0,0,0}][0][0](1,0)<<"     "<<pi_k[freq_pair.first][{0,0,0}][0][0](1,1)<<"     "<<pi_k[freq_pair.first][{0,0,0}][0][0](2,2)<<endl;
    
    cal_MP2_energy_pi_k(freq_grid);
    //cal_MP2_energy_pi_k_tau(time_grid,cal_pi_k_tau());
    RPA_correlation_energy(freq_grid);

    double t_end=omp_get_wtime();
    cout<<"TOTAL TIME USED :  "<<t_end-t_begin<<endl;
    cout<<"ONLY CHI0 TIME USED:  "<<t_chi0_end-t_chi0_begin<<endl;
    cout<<"  Average per_task time:"<<time_task_tot/(tot_pair.size());
   // cout<<" periodic_exx_lcao.enegy()"<<exx_lcao.energy<<endl;
    //cout<<" periodic_en.etot"<<en.etot<<endl;

    
}

vector<vector<pair<double,Vector3_Order<int>>>> Cal_Periodic_Chi0::process_task_tau_R(const map<double,double> &time_grid, const set<Vector3_Order<int>> &R_grid)
{
    int tot_num=time_grid.size()*R_grid.size();
    int p_num=para_mpi.get_size();
    vector<vector<pair<double,Vector3_Order<int>>>> p_task(p_num);
    int itR=0;
    for(auto &tau:time_grid)
        for(auto &R:R_grid)
        {   
            int p_id=itR%p_num;
            printf(" itR:  %d,     p_id:   %d \n",itR,p_id);
            p_task[p_id].push_back(pair<double,Vector3_Order<int>>(tau.first,R));
            itR++;
        }
    return p_task;
}

vector<pair<size_t,size_t>> Cal_Periodic_Chi0::process_task_ap_local()
{
    vector<pair<size_t,size_t>> tot_pair(get_atom_pair(Vq));
    int tot_num=tot_pair.size();
    int p_num=para_mpi.get_size();
    int p_local=para_mpi.get_myid();
    vector<pair<size_t,size_t>> p_task_ap;
    int itp=0;
    for(const auto &ap:tot_pair)
    {
        int p_id=itp%p_num;
        if(p_id==p_local) 
        {
            cout<<"  p_id:  "<<p_id<<"   I J :  "<<ap.first<<"   "<<ap.second<<endl;
            p_task_ap.push_back(ap);
        }
        itp++;
    }

    return p_task_ap;
}


vector<pair<size_t,size_t>> Cal_Periodic_Chi0::get_atom_pair(const map<size_t,map<size_t,map<Vector3_Order<double>,std::shared_ptr<ComplexMatrix>>>> &m)
{
    vector<pair<size_t,size_t>> atom_pairs;
    for(const auto &mA : m)
		for(const auto &mB : mA.second)
			atom_pairs.push_back({mA.first,mB.first});
	printf("  tot atom_pairs tasks: %d\n",atom_pairs.size());
    return atom_pairs;
}


double Cal_Periodic_Chi0::get_band_gap()
{
    double E_homo=-1e6;
    double E_lumo=1e6;
    for(int ik=0;ik!=n_kpoints;ik++)
    {
        int LUMO_level=0;
        for(int n=0;n!=NBANDS;n++)
        {
            if(wg(ik,n)>=1e-12)
            {
                LUMO_level+=1;
            }
		}
        if(LUMO_level>=1)
        {
            if(E_homo<ekb[ik][LUMO_level-1])
                E_homo=ekb[ik][LUMO_level-1];
        }
        if(E_lumo>ekb[ik][LUMO_level])
            E_lumo=ekb[ik][LUMO_level];
        //cout<<"    ik:"<<ik<<"     LUMO_level:"<<LUMO_level<<"     E_homo:<<"E_homo<<"      E_lumo:"<<E_lumo<<endl;
	}
    double gap;
    //gap=0.5*(wf.ekb[0][LUMO_level]-wf.ekb[0][LUMO_level-1]);
    gap=0.5*(E_lumo-E_homo);
	std::cout<<"E_homo(Ha):"<<E_homo*0.5<<"     E_lumo(Ha):"<<E_lumo*0.5<<"     gap(Ha): "<<gap<<endl;
    return gap;
    
}
void Cal_Periodic_Chi0::init()
{
    //para_mpi.mpi_init(argc,argv);
    temp_Cs_count=0;
    //cout<<" knpoints_num= "<<n_kpoints<<endl;
    //green_k.resize(n_kpoints);
/*
    wfc_k.resize(n_kpoints);
    //cout<<"begin to init wfc_k"<<endl;
    for(int i=0;i!=wfc_k.size();++i)
    {
        wfc_k[i].create(LOC.wfc_dm_2d.wfc_k[i].nr,LOC.wfc_dm_2d.wfc_k[i].nc);
        //cout<<"nc="<<LOC.wfc_dm_2d.wfc_k[i].nr<<"   nc="<<LOC.wfc_dm_2d.wfc_k[i].nc<<endl;
        //cout<<"after create "<<endl;
        //cout<<"wfc_k[i].size="<<wfc_k[i].size<<"   LOC.wfc_dm_2d.wfc_k[i].size="<<LOC.wfc_dm_2d.wfc_k[i].size<<endl;
        for(int j=0;j!=LOC.wfc_dm_2d.wfc_k[i].size;++j)
        {
            wfc_k[i].c[j]=LOC.wfc_dm_2d.wfc_k[i].c[j];
        }
        cout<<"   ik="<<i<<"  ";
        //print_complex_matrix("  wfc_k_mat ",LOC.wfc_dm_2d.wfc_k[i].nc,LOC.wfc_dm_2d.wfc_k[i].nr,LOC.wfc_dm_2d.wfc_k[i].c,LOC.wfc_dm_2d.wfc_k[i].nc);
    }
    //cout<<"end init wfc_k"<<endl;
    
    */

}

set<Vector3_Order<int>> Cal_Periodic_Chi0::construct_R_grid()
{
	//cout<<" begin to construct_R_grid"<<endl;
    set<Vector3_Order<int>> R_grid;
    R_grid.clear();
    
    R_grid.insert({0,0,0});
	for( int x=-(kv_nmp[0])/2; x<=(kv_nmp[0]-1)/2; ++x )
		for( int y=-(kv_nmp[1])/2; y<=(kv_nmp[1]-1)/2; ++y )
			for( int z=-(kv_nmp[2])/2; z<=(kv_nmp[2]-1)/2; ++z )
				R_grid.insert({x,y,z});
 
    
 
    R_periodic = Vector3<int>{kv_nmp[0],kv_nmp[1],kv_nmp[2]};
    for(const auto &R:R_grid)
        cout<<R<<endl;
    cout<<" R_periodic: "<<R_periodic<<endl; 
	return R_grid;
}



std::vector<ComplexMatrix> Cal_Periodic_Chi0::cal_Green_func_element_k( const matrix &wg, const double time_tau )
{
	//cout<<wg.nr<<"   "<<wg.nc<<"  "<<time_tau<<endl;
    // dm = wfc.T * wg * wfc.conj()
	// dm[ik](iw1,iw2) = \sum_{ib} wfc[ik](ib,iw1).T * wg(ik,ib) * wfc[ik](ib,iw2).conj()
    //assert(wg.nc<=NLOCAL);
	//assert(wg.nr==n_kpoints);
    std::vector<ComplexMatrix> green_k(n_kpoints);
	for(int ik=0; ik!=n_kpoints; ++ik)
	{			
        //cout<<" ik= "<<ik<<endl;
		std::vector<double> wg_local(wg.nc,0.0);
		for(int ib_global=0; ib_global!=wg.nc; ++ib_global)
		{
			//const int ib_local = ParaO.trace_loc_col[ib_global];
            double temp_e_up=-1*0.5*(ekb[ik][ib_global]-efermi)*time_tau;
            double ephase;
            if(temp_e_up >= 0)
            {
                ephase=1.0;
            }
            else
            {
                ephase = std::exp(temp_e_up);
            }
            wg_local[ib_global] = wg(ik,ib_global)*ephase;
          
        }
		// wg_wfc(ib,iw) = wg[ib] * wfc(ib,iw).conj();
       // cout<<wfc_k.size()<<endl;
       // cout<<"wfc_k  dim "<<wfc_k[ik].nr<<wfc_k[ik].nc<<endl;
		ComplexMatrix wg_wfc = conj(wfc_k[ik]);
       // cout<<" wg_wfc "<<wg_wfc.nr<<"  "<<wg_wfc.nc<<endl;
		for(int ir=0; ir!=wg_wfc.nr; ++ir)
			LapackConnector::scal( wg_wfc.nc, wg_local[ir], wg_wfc.c+ir*wg_wfc.nc, 1 );
		// C++: dm(iw1,iw2) = wfc(ib,iw1).T * wg_wfc(ib,iw2)
		green_k[ik].create( wfc_k[ik].nr, wfc_k[ik].nc );
        green_k[ik]=transpose(wfc_k[ik],0)*wg_wfc;
	}
    return green_k;
}

void Cal_Periodic_Chi0::cal_Green_func_R_tau( const double &time_tau, const Vector3_Order<int> &R,  const Vector3<double>* kvec_c)
{
    Vector3_Order<double> R_2(R.x,R.y,R.z);
   
    Vector3_Order<int> R_0(0,0,0);
    //cout<<"R_0 init"<<endl;
    //cout<<"R= "<<R<<"     time_tau= "<<time_tau<<endl;
    
    //cout<<"   zero_out";
    //cout<<"Green_function_zero_out"<<endl;
    std::vector<ComplexMatrix> green_k;
    if( time_tau>0 ) 
    {
        matrix i_mat(wg.nr,wg.nc);
        for(int i=0;i!=i_mat.nr*i_mat.nc;++i)
            i_mat.c[i]=1.0/n_kpoints*NSPIN; 
        //cal_Green_func_element_k(i_mat-wf.wg, time_tau);
        matrix wg_unocc=i_mat-wg;
        for(int i=0;i!=wg_unocc.nr*wg_unocc.nc;++i)
        {
            if(wg_unocc.c[i]<0)
                wg_unocc.c[i]=0;
        }
        //print_matrix("wg_unocc_mat:",wg_unocc);
        //cout<<"cal green ele"<<endl;
        green_k=cal_Green_func_element_k(wg_unocc, time_tau);
    } 
    else
    {
        green_k=cal_Green_func_element_k( wg, time_tau);
    }
	
    int iter_nks=0;
    for(int is=0;is!=NSPIN;is++)
    {
        ComplexMatrix Green_glo_tmp(green_k[0].nr,green_k[0].nc);
        Green_glo_tmp.zero_out();
        for(int ik=0; ik!=n_kpoints/NSPIN; ++ik)
        {
            const double arg =-1*( kvec_c[ik+iter_nks]* (R_2*latvec)) * TWO_PI; //latvec
            const complex<double> kphase = complex <double> ( cos(arg),  sin(arg) );
            Green_glo_tmp += green_k[ik+iter_nks] * kphase ;
        }

        if(time_tau<=0){Green_glo_tmp*=-1;}
        iter_nks+=n_kpoints/NSPIN;
        for (int I=0;I!=natom;I++)
        {
            const size_t I_num=atom_nw[I];
            for(int J=0;J!=natom;J++)
            {
                const size_t J_num=atom_nw[J];
                Green_atom[is][I][J][time_tau][R].create(I_num,J_num);
                for(size_t i=0;i!=I_num;i++)
                {
                    size_t i_glo=atom_iw_loc2glo(I,i);
                    for(size_t j=0;j!=J_num;j++)
                    {
                        size_t j_glo=atom_iw_loc2glo(J,j);
                        Green_atom[is][I][J][time_tau][R](i,j)=Green_glo_tmp(i_glo,j_glo).real();
                    }
                }
            }
        }
    }
}


void Cal_Periodic_Chi0::Cosine_to_chi0_freq(map<double,double> &time_grid, map<double,double> &freq_grid )
{
    //map<double,double> freq_grid(construct_gauss_grid(120));  
    int n_freq=freq_grid.size();
    int n_tau=time_grid.size();
    assert(n_freq==n_tau);
    vector<double> tran_gamma=read_cosine_trans_grid( minimax_grid_path+"/time2freq_transform_grid.txt");
    for(const auto &time_pair:chi0)
    {
        //const double time=time_pair.first;
        const auto chi0_time=time_pair.second; 
        for(const auto &R_pair:chi0_time)
        {
            const auto R=R_pair.first;
            const auto chi0_time_R=R_pair.second;
            for(const auto &I_pair:chi0_time_R)
            {
                const size_t I=I_pair.first;
                const auto chi0_time_R_I=I_pair.second;
                for(const auto &J_pair:chi0_time_R_I)
                {
                    const size_t J=J_pair.first;
                    const auto mat=J_pair.second;
                    for(const auto &freq_pair:freq_grid)
                    {
                        const double freq=freq_pair.first;
                        const double freq_weight=freq_pair.second;
                        chi0_freq[freq][R][I][J].create(mat.nr,mat.nc);
                        
                    }
                }
            }
        }
        break;
    }
    
    int iter_tau=0;    
    for(const auto &time_pair:chi0)
    {
        const double time=time_pair.first;
        const auto chi0_time=time_pair.second;
        for(const auto &R_pair:chi0_time)
        {
            const auto R=R_pair.first;
            const auto chi0_time_R=R_pair.second;
            for(const auto &I_pair:chi0_time_R)
            {
                const size_t I=I_pair.first;
                const auto chi0_time_R_I=I_pair.second;
                for(const auto &J_pair:chi0_time_R_I)
                {
                    const size_t J=J_pair.first;
                    const auto mat=J_pair.second;
                    int iter_freq=0;
                    for(const auto &freq_pair:freq_grid)
                    {
                        
                        const double freq=freq_pair.first;
                        const double cos_weight=tran_gamma[iter_freq*n_tau+iter_tau];
                        chi0_freq.at(freq).at(R).at(I).at(J)+=mat*cos_weight;  
                        //cout<<time<<"   "<<time_weight<<"    freq:"<<freq<<"   cos_weight"<<cos_weight  \
                            <<"      mat:  "<<mat(0,0)<<"    "<<chi0_freq[freq][R][I][J](0,0)<<endl;
                        iter_freq++;
                    }
                    
                }
            }
        }
        iter_tau++;
    }
    cout<<"finish cosine tran chi0 !"<<endl;
}

matrix Cal_Periodic_Chi0:: cal_chi0_element(const double &time_tau, const Vector3_Order<int> &R, const size_t &I_index, const size_t &J_index)
{
    //printf("     begin chi0  thread: %d,  I: %d, J: %d\n",omp_get_thread_num(),I_index,J_index);
    //for(const auto &K_pair:Cs[I_index])
    //Vector3_Order<int> R_0(0,0,0);
    
    const size_t i_num=atom_nw[I_index];
    const size_t mu_num=atom_mu[I_index];

    const size_t j_num=atom_nw[J_index];
    const size_t nu_num=atom_mu[J_index];

    chi0[time_tau][R][I_index][J_index].create(mu_num,nu_num);
    matrix O_sum(mu_num,nu_num);
    for(const auto &K_pair:Cs[I_index])
    {
        const auto K_index=K_pair.first;
        const size_t k_num=atom_nw[K_index];
        
        for(const auto &R1_index:K_pair.second)
        {
            const auto R1=R1_index.first; 
            const auto &Cs_mat1=R1_index.second;
           // cout<<"R1:  begin  "<<R1<<endl;
            for(int is=0;is!=NSPIN;is++)
            {
                matrix N_R2(k_num,j_num*nu_num);
                matrix N_conj_R2(k_num,j_num*nu_num);
            
                matrix X_R2(i_num,j_num*nu_num);
                matrix X_conj_R2(i_num,j_num*nu_num);
            
                for(const auto &L_pair:Cs[J_index])
                {
                    const auto L_index=L_pair.first;
                    const size_t l_num=atom_nw[L_index];
                    
                    for(const auto &R2_index:L_pair.second)
                    {
                        const auto R2=R2_index.first; 
                        const auto &Cs_mat2=R2_index.second;
                
                        Vector3_Order<int> R_temp_1(R+R2-R1);
                        Vector3_Order<int> R_temp_2(R2+R); 
                    
                        assert(j_num*l_num==(*Cs_mat2).nr); 
                        //printf("          thread: %d, IJKL:   %d,%d,%d,%d  R:(  %d,%d,%d  )  tau:%f\n",omp_get_thread_num(),I_index,J_index,K_index,L_index,R.x,R.y,R.z,time_tau);
                        matrix Cs2_reshape(reshape_Cs(j_num,l_num,nu_num,Cs_mat2));
                        //printf("          thread: %d, Green_K %d,%d,   Cs  %d,%d\n",omp_get_thread_num(),Green_atom[is][K_index][L_index][time_tau][R_temp_1%R_periodic].nr,Green_atom[is][K_index][L_index][time_tau][R_temp_1%R_periodic].nc,Cs2_reshape.nr,Cs2_reshape.nc);
                        N_R2+=Green_atom.at(is).at(K_index).at(L_index).at(time_tau).at(R_temp_1%R_periodic)*Cs2_reshape;
                        //printf("          thread: %d, 1.2\n",omp_get_thread_num());
                        N_conj_R2+=Green_atom.at(is).at(K_index).at(L_index).at(-time_tau).at(R_temp_1%R_periodic)*Cs2_reshape;
                        //printf("          thread: %d, 1.3\n",omp_get_thread_num());
                        X_R2+=Green_atom.at(is).at(I_index).at(L_index).at(time_tau).at(R_temp_2%R_periodic)*Cs2_reshape;
                        //printf("          thread: %d, 1.4\n",omp_get_thread_num());
                        X_conj_R2+=Green_atom.at(is).at(I_index).at(L_index).at(-time_tau).at(R_temp_2%R_periodic)*Cs2_reshape;
                       
                    }
                }
                //printf("          thread: %d, 2\n",omp_get_thread_num());
                matrix N_R2_rs(reshape_mat(k_num,j_num,nu_num,N_R2));
                matrix N_conj_R2_rs(reshape_mat(k_num,j_num,nu_num,N_conj_R2));
                matrix X_R2_rs(reshape_mat(i_num,j_num,nu_num,X_R2));
                matrix X_conj_R2_rs(reshape_mat(i_num,j_num,nu_num,X_conj_R2));
    
                Vector3_Order<int> R_temp_3(R-R1);
                //cal and sum O
                matrix O(i_num,k_num*nu_num);
                matrix Z(k_num,i_num*nu_num);
                //printf("          thread: %d, 3\n",omp_get_thread_num());
                O+=Green_atom.at(is).at(I_index).at(J_index).at(-time_tau).at(R)*N_R2_rs;
                O+=Green_atom.at(is).at(I_index).at(J_index).at(time_tau).at(R)*N_conj_R2_rs;
                Z+=Green_atom.at(is).at(K_index).at(J_index).at(-time_tau).at(R_temp_3%R_periodic)*X_R2_rs;
                Z+=Green_atom.at(is).at(K_index).at(J_index).at(time_tau).at(R_temp_3%R_periodic)*X_conj_R2_rs;
             
                matrix Z_rs(reshape_mat(k_num,i_num,nu_num,Z));
                
                O+=Z_rs;
                matrix OZ(reshape_mat_21(i_num,k_num,nu_num,O));
                matrix Cs1_tran(transpose(*Cs_mat1));
                O_sum+=Cs1_tran*OZ;
             //   cout<<"   K, R1:   "<<K_index<<"   "<<R1;
                //rt_m_max(O_sum);
            }
        }
    }
    //chi0[time_tau][R][I_index][J_index]=O_sum;
    //printf("     finish chi0  thread: %d,  I: %d, J: %d",omp_get_thread_num(),I_index,J_index);
    return O_sum;
} 

matrix Cal_Periodic_Chi0::reshape_Cs(const size_t n1, const size_t n2, const size_t n3, const shared_ptr<matrix> &Cs) //(n1*n2,n3) -> (n2,n1*n3)
{
    const auto length = sizeof(double)*n3;
	const auto n13 = n1*n3;
	const double *m_ptr = (*Cs).c;
	matrix m_new( n2, n1*n3, false );
	for( size_t i1=0; i1!=n1; ++i1 )
	{
		double *m_new_ptr = m_new.c+i1*n3;
		for( size_t i2=0; i2!=n2; ++i2, m_ptr+=n3, m_new_ptr+=n13 )
			memcpy( m_new_ptr, m_ptr, length );
	}
    return m_new;
}

matrix Cal_Periodic_Chi0::reshape_dim_Cs(const size_t n1, const size_t n2, const size_t n3, const shared_ptr<matrix> &Cs) //(n1*n2,n3) -> (n1,n2*n3)
{
    const auto length = sizeof(double)*n1*n2*n3;
	const auto n13 = n1*n3;
	const double *m_ptr = (*Cs).c;
	matrix m_new( n1, n2*n3, false );
    double *m_new_ptr = m_new.c;
    memcpy( m_new_ptr, m_ptr, length );
    return m_new;
}


matrix Cal_Periodic_Chi0::reshape_mat(const size_t n1, const size_t n2, const size_t n3, const matrix &mat) //(n1,n2*n3) -> (n2,n1*n3)
{
    const auto length = sizeof(double)*n3;
	const auto n13 = n1*n3;
	const double *m_ptr = mat.c;
	matrix m_new( n2, n1*n3, false );
	for( size_t i1=0; i1!=n1; ++i1 )
	{
		double *m_new_ptr = m_new.c+i1*n3;
		for( size_t i2=0; i2!=n2; ++i2, m_ptr+=n3, m_new_ptr+=n13 )
			memcpy( m_new_ptr, m_ptr, length );
	}
    return m_new;
}
matrix Cal_Periodic_Chi0::reshape_mat_21(const size_t n1, const size_t n2, const size_t n3, const matrix &mat) //(n1,n2*n3) -> (n1*n2,n3)
{
    const auto length = sizeof(double)*n1*n2*n3;
	const auto n13 = n1*n2*n3;
	const double *m_ptr = mat.c;
	matrix m_new( n1*n2 ,n3, false );
	double *m_new_ptr = m_new.c;
    memcpy( m_new_ptr, m_ptr, length );
    return m_new;
}


std::map<double,map<Vector3_Order<double>,map<size_t,map<size_t,ComplexMatrix>>>> Cal_Periodic_Chi0::cal_pi_k_use_aims_vq()
{
    cout<<" bigin to cal_pi_k"<<endl;
    for(auto &freq_pair:chi0_freq)
    {
        const double freq=freq_pair.first;
        const auto chi0_freq=freq_pair.second;
        for(auto &R_pair:chi0_freq)
        {
            const Vector3_Order<int> R=R_pair.first;
            auto chi0_freq_R=R_pair.second; 
            for(auto &I_pair:chi0_freq_R)
            {
                const size_t I=I_pair.first;
                const size_t I_mu=atom_mu[I];
                auto chi0_freq_R_I=I_pair.second;
                //const size_t mu_num=exx_lcao.index_abfs[ucell.iat2it[I]].count_size;
                for(auto &J_pair:chi0_freq_R_I)
                {
                    size_t J=J_pair.first;
                    const size_t J_mu=atom_mu[J];
                    //const auto chi0_mat=J_pair.second;
                    //(*Vps[I][J][R]).c
                    //cout<<"2"<<endl;
                   // const size_t nu_num=exx_lcao.index_abfs[ucell.iat2it[J]].count_size;
                    for(auto &q_pair:irk_weight)
                    {
                        //cout<<"3"<<endl;
                        Vector3_Order<double> ik_vec=q_pair.first;
                        //auto &vq_mat=*q_pair.second;
                        pi_k[freq][ik_vec][I][J].create(I_mu,J_mu);
                        if(I!=J)
                            pi_k[freq][ik_vec][J][I].create(J_mu,I_mu);
                        //cout<<freq<<"    "<<I<<"   "<<J<<"    ik:"<<ik<<"   ik_vec "<<ik_vec<<endl;
                    }
                }
            }
            break;
        }
    }

    for(auto &freq_p:chi0_k)
    {
        const double freq=freq_p.first;
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
                        
                        // if(freq==first_freq)
                        // {
                        //     //cout<<I<<"  "<<Q<<"  "<<J<<"   "<<ik_vec;
                        //     //rt_cm_max(Vq_mat);
                        //     // print_complex_matrix("vq",Vq_mat);
                        //     // if(J<=Q)
                        //     // {
                        //     //     print_complex_matrix("chi0_mat",chi0_k[freq][ik_vec][J][Q]);
                        //     //   //  print_complex_matrix("pi_mat",Vq_mat*chi0_k[freq][ik_vec][J][Q]);
                        //     // }
                        //     // else
                        //     // {
                        //     //     print_complex_matrix("chi0_mat",transpose(chi0_k[freq][ik_vec][Q][J],1));
                        //     //     //print_complex_matrix("pi_mat",Vq_mat*transpose(chi0_k[freq][ik_vec][Q][J],1));
                        //     // }
                        // }
                        if(J<=Q)
                        {
                            //cout<<"   1";
                            pi_k.at(freq).at(ik_vec).at(I).at(Q)+=Vq_mat*chi0_k.at(freq).at(ik_vec).at(J).at(Q);
                            // if(freq==first_freq)
                            // {
                            //     // cout<<" pi ";
                            //     // rt_cm_max(pi_k.at(freq).at(ik_vec).at(I).at(Q));
                            //     // cout<<" chi0 ";
                            //     // rt_cm_max(chi0_k.at(freq).at(ik_vec).at(J).at(Q));
                            //     // cout<<" Vq ";
                            //     // rt_cm_max(Vq_mat);
                            // }
                            // if(freq==first_freq && I==natom-1 && Q==natom-1)
                            // {
                            //     cout<<"IQJ: "<<I<<"   "<<Q<<"   "<<J<<"    pi_ele: "<<pi_k.at(freq).at(ik_vec).at(I).at(Q)(0,0)<<endl;
                            //     rt_cm_max(chi0_k.at(freq).at(ik_vec).at(J).at(Q));
                            //     rt_cm_max(Vq_mat);
                            // }
                            
                        }
                        else
                        {
                            //cout<<"   2";
                            pi_k.at(freq).at(ik_vec).at(I).at(Q)+=Vq_mat*transpose(chi0_k.at(freq).at(ik_vec).at(Q).at(J),1);
                            // if(freq==first_freq && I==natom-1 && Q==natom-1)
                            // {
                            //     cout<<"IQJ: "<<I<<"   "<<Q<<"   "<<J<<"    pi_ele: "<<pi_k.at(freq).at(ik_vec).at(I).at(Q)(0,0)<<endl;
                            //     rt_cm_max(chi0_k.at(freq).at(ik_vec).at(Q).at(J));
                            //     rt_cm_max(Vq_mat);
                            // }
                        }
                        
                        if(I!=J)
                        {
                            if(I<=Q)
                            {
                               // cout<<"   3"<<"  pi: "<<pi_k[freq][ik_vec][J][Q].nr<<pi_k[freq][ik_vec][J][Q].nc<<"  vq nr nc: "<<(transpose(Vq_mat,1)).nr<<"  "<<(transpose(Vq_mat,1)).nc<<"   chi0:"<<chi0_k[freq][ik_vec][I][Q].nr<<"  "<<chi0_k[freq][ik_vec][I][Q].nc<<endl;
                                pi_k.at(freq).at(ik_vec).at(J).at(Q)+=transpose(Vq_mat,1)*chi0_k.at(freq).at(ik_vec).at(I).at(Q);
                                // if(freq==first_freq && J==natom-1 && Q==natom-1)
                                // {
                                //     cout<<"IQJ: "<<I<<"   "<<Q<<"   "<<J<<"    pi_ele: "<<pi_k.at(freq).at(ik_vec).at(J).at(Q)(0,0)<<endl;
                                //     rt_cm_max(chi0_k.at(freq).at(ik_vec).at(I).at(Q));
                                //     rt_cm_max(Vq_mat);
                                // }
                            }
                            else
                            {
                                //cout<<"   4";
                                pi_k.at(freq).at(ik_vec).at(J).at(Q)+=transpose(Vq_mat,1)*transpose(chi0_k.at(freq).at(ik_vec).at(Q).at(I),1);
                                // if(freq==first_freq && J==natom-1 && Q==natom-1)
                                // {
                                //     cout<<"IQJ: "<<I<<"   "<<Q<<"   "<<J<<"    pi_ele: "<<pi_k.at(freq).at(ik_vec).at(J).at(Q)(0,0)<<endl;
                                //     rt_cm_max(chi0_k.at(freq).at(ik_vec).at(Q).at(I));
                                //     rt_cm_max(Vq_mat);
                                // }
                            }
                        }   
                             
                    }
                }
            }
        }
    }
  //  print_complex_matrix("  last_pi_mat:",pi_k.at(first_freq).at({0,0,0}).at(natom-1).at(natom-1));
    cout<<"Finish cal pi_k"<<endl;
    return pi_k;
}

void Cal_Periodic_Chi0::FT_R_to_K_chi0()
{
    cout<<" bigin to FT_R_to_K_chi0"<<endl;
    for(auto &freq_pair:chi0_freq)
    {
        const double freq=freq_pair.first;
        const auto chi0_freq=freq_pair.second;
        for(auto &R_pair:chi0_freq)
        {
            const Vector3_Order<int> R=R_pair.first;
            auto chi0_freq_R=R_pair.second; 
            for(auto &I_pair:chi0_freq_R)
            {
                const size_t I=I_pair.first;
                auto chi0_freq_R_I=I_pair.second;
                //const size_t mu_num=exx_lcao.index_abfs[ucell.iat2it[I]].count_size;
                for(auto &J_pair:chi0_freq_R_I)
                {
                    size_t J=J_pair.first;
                    //const auto chi0_mat=J_pair.second;
                    //(*Vps[I][J][R]).c
                    //cout<<"2"<<endl;
                   // const size_t nu_num=exx_lcao.index_abfs[ucell.iat2it[J]].count_size;
                    for(auto &q_pair:irk_weight)
                    {
                        //cout<<"3"<<endl;
                        Vector3_Order<double> ik_vec=q_pair.first;
                        //auto &vq_mat=*q_pair.second;
                        chi0_k[freq][ik_vec][I][J].create(atom_mu[I],atom_mu[J]);
                        //cout<<freq<<"    "<<I<<"   "<<J<<"    ik:"<<ik<<"   ik_vec "<<ik_vec<<endl;
                    }
                }
            }
            break;
        }
    }
    cout<<"  finish init chi0_k"<<endl;
    temp_chi0_k_count=0;
    for(auto &freq_pair:chi0_freq)
    {
        const double freq=freq_pair.first;
        const auto chi0_freq=freq_pair.second;

        for(auto &R_pair:chi0_freq)
        {
            const Vector3_Order<int> R=R_pair.first;
            auto chi0_freq_R=R_pair.second; 
            for(auto &I_pair:chi0_freq_R)
            {
                const size_t I=I_pair.first;
                auto chi0_freq_R_I=I_pair.second;
                const size_t mu_num=atom_mu[I];
                for(auto &J_pair:chi0_freq_R_I)
                {
                    size_t J=J_pair.first;
                    const size_t nu_num=atom_mu[J];
                    auto &chi0_mat=J_pair.second;
                    size_t mat_size=chi0_mat.nr*chi0_mat.nc;
                    for(auto &ikp:chi0_k[freq])
                    {
                        auto ik_vec=ikp.first;
                        const double arg = ( ik_vec* (R*latvec)) * TWO_PI;
                        const complex<double> kphase = complex <double> ( cos(arg),  sin(arg) );
                        //cout<<"  I J:  "<<I<<"  "<<J<<endl;
                        //cout<<chi0_k[freq][ikp.first][I][J].nr<<"  "<<chi0_k[freq][ikp.first][I][J].nc<<"  "<<chi0_mat.nr<<"  "<<chi0_mat.nc<<endl;
                        for(size_t ci=0;ci!=mat_size;ci++)
                            chi0_k.at(freq).at(ikp.first).at(I).at(J).c[ci]+=kphase*chi0_mat.c[ci];
                    }
                }
            }
        }
        temp_chi0_k_count++;
    }
}


void Cal_Periodic_Chi0::cal_MP2_energy_pi_k( map<double,double> &freq_grid ) 
{
    //std::map<double,map<Abfs::Vector3_Order<double>,map<size_t,map<size_t,ComplexMatrix>>>> pi_k(cal_pi_k());
    std::map<double,map<Vector3_Order<double>,map<size_t,map<size_t,ComplexMatrix>>>> pi_k_square;
    for(auto &freq_pair:pi_k)
    {
        const double freq=freq_pair.first;
        const auto pi_omega=freq_pair.second;
        for(auto &k_pair:pi_omega)
        {
            const auto k=k_pair.first;
            auto pi_omega_k=k_pair.second; 
            for(auto &I_pair:pi_omega_k)
            {
                const size_t I=I_pair.first;
                auto pi_omega_k_I=I_pair.second;
                const size_t mu_num=atom_mu[I];
                pi_k_square[freq][k][I][I].create(mu_num,mu_num);
                for(auto &J_pair:pi_omega_k_I)
                {
                    const size_t J=J_pair.first;
                    const auto pi_mat=J_pair.second;
                    const auto pi2_mat=pi_omega_k[J][I];
                    const size_t nu_num=atom_mu[J];
                   // cout<<freq<<"    "<<R<<"   "<<I<<"    "<<J<<"    mu_num="<<mu_num<<"    nu_num="<<nu_num<<endl;
                    
                    
                    const size_t xi_num=atom_mu[J];
                    
                    LapackConnector::gemm(
                        'N','N',
                        mu_num,mu_num,xi_num,
                        1,
                        pi_mat.c,xi_num,
                        pi2_mat.c,mu_num,
                        1,
                        pi_k_square[freq][k][I][I].c, 
                        mu_num);
             
                }
                //print_complex_matrix("pi_square_mat",pi_square[freq][R][I][I].nc,pi_square[freq][R][I][I].nr,pi_square[freq][R][I][I].c,pi_square[freq][R][I][I].nc);
            }
        }
    }
    
    complex<double> tot_trace_pi(0.0,0.0);
    for(auto &freq_pair:pi_k_square) 
    {
        auto freq=freq_pair.first;
        auto pi_omega=freq_pair.second;
        complex<double> trace_pi(0.0,0.0);
        for(auto &k_pair:pi_omega)
        {
            auto pi_omega_k=k_pair.second;
            for(auto &I_pair:pi_omega_k)
            {
                const size_t I=I_pair.first;
                auto pi_omega_k_I=I_pair.second;
                trace_pi+=trace(pi_omega_k_I[I]);
            }
        }
        tot_trace_pi+=trace_pi*freq_grid[freq]*NSPIN/n_kpoints/TWO_PI/2*(-1);
        cout<<" pi_k   freq:"<<freq<<"   trace:"<<trace_pi<<endl;
    }
    cout<<" tot_MP2_energy_pi_k: "<<tot_trace_pi;
}

void Cal_Periodic_Chi0::RPA_correlation_energy( map<double,double> &freq_grid )
{
    int range_all=0;

    for(auto &iat:atom_mu)
    {
        range_all+=iat.second;
    }
    
    vector<int> part_range;
    part_range.resize(atom_mu.size());
	part_range[0]=0;
	int count_range=0;
	for (int I=0;I!=atom_mu.size()-1;I++)
	{
		count_range+=atom_mu[I];
		part_range[I+1]=count_range;
	}
	
	cout<<"part_range:"<<endl;
	for(int I=0;I!=atom_mu.size();I++)
	{
		cout<<part_range[I]<<endl;
	}
	cout<<"part_range over"<<endl;
    //std::map<double,map<Abfs::Vector3_Order<double>,map<size_t,map<size_t,ComplexMatrix>>>> pi_k(cal_pi_k());
    std::map<double,map<Vector3_Order<double>,ComplexMatrix>> pi_freq_2D;
    // std::map<double,map<Vector3_Order<double>,ComplexMatrix>> chi0_2D;
    // std::map<Vector3_Order<double>,ComplexMatrix> Vq_2D;
    for(const auto &freq_pair:pi_k)
    {
        const auto freq=freq_pair.first;
        const auto pi_freq=freq_pair.second;
        
        for(const auto &k_pair:pi_freq)
        {
            const auto kvec_c=k_pair.first;
            const auto pi_freq_k=k_pair.second;
           // chi0_2D[freq][kvec_c].create(range_all,range_all);
            pi_freq_2D[freq][kvec_c].create(range_all,range_all);
        
            for(const auto &I_pair:pi_freq_k)
            {
                const auto I=I_pair.first;
                const auto pi_freq_k_I=I_pair.second;
                const size_t mu_num=atom_mu[I];
                for(const auto &J_pair:pi_freq_k_I)
                {
                    const auto J=J_pair.first;
                    const auto mat=J_pair.second;
                    const size_t nu_num=atom_mu[J];
                    for(size_t mu=0;mu!=mu_num;++mu)
                    {
                        for(size_t nu=0;nu!=nu_num;++nu)
                        {
                            pi_freq_2D.at(freq).at(kvec_c)(part_range[I]+mu,part_range[J]+nu)+=mat(mu,nu); 
                            
                        }
                    }
                    // if(I<=J)
                    // {
                    //     for(size_t mu=0;mu!=mu_num;++mu)
                    //     {
                    //         for(size_t nu=0;nu!=nu_num;++nu)
                    //         {
                    //             chi0_2D.at(freq).at(kvec_c)(part_range[I]+mu,part_range[J]+nu)+=chi0_k.at(freq).at(kvec_c).at(I).at(J)(mu,nu);
                    //         }
                    //     }
                    // }
                    // else
                    // {
                    //     for(size_t mu=0;mu!=mu_num;++mu)
                    //     {
                    //         for(size_t nu=0;nu!=nu_num;++nu)
                    //         {
                    //             chi0_2D.at(freq).at(kvec_c)(part_range[I]+mu,part_range[J]+nu)+=conj(chi0_k.at(freq).at(kvec_c).at(J).at(I)(nu,mu));
                    //         }
                    //     }
                    // }
                }
            }
            //cout<<"  freq:   "<<freq;
            //rt_cm_max(pi_freq_2D.at(freq).at(kvec_c));
            // if(freq==first_freq)
            // {
            //     cout<<"  pi_2d mat , ik_vec:  "<<kvec_c<<endl;
            //     print_complex_matrix("pi_2d",pi_freq_2D[freq][kvec_c]);
            // }
        }
    } 
    // for(auto &qvec:irk_weight)
    // {
    //     auto kvec_c=qvec.first;
    //     Vq_2D[kvec_c].create(range_all,range_all); 
    //     for(auto &I_p:pi_k[first_freq][kvec_c])
    //     {
    //         const auto I=I_p.first;
    //         const size_t mu_num=atom_mu[I];
    //         for(auto &J_p:I_p.second)
    //         {
    //             const auto J=J_p.first;
    //             const size_t nu_num=atom_mu[J];
    //             if(I<=J)
    //             {
    //                 for(size_t mu=0;mu!=mu_num;++mu)
    //                 {
    //                     for(size_t nu=0;nu!=nu_num;++nu)
    //                     {
    //                         Vq_2D[kvec_c](part_range[I]+mu,part_range[J]+nu)+=(*Vq[I][J][kvec_c])(mu,nu);
    //                     }
    //                 }
    //             }
    //             else
    //             {
    //                 for(size_t mu=0;mu!=mu_num;++mu)
    //                 {
    //                     for(size_t nu=0;nu!=nu_num;++nu)
    //                     {
    //                         Vq_2D[kvec_c](part_range[I]+mu,part_range[J]+nu)+=conj((*Vq[J][I][kvec_c])(nu,mu));
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }

    // ofstream fs;
    // fs.open("librpa_pi.txt");
    // for(auto &irk:irk_weight)
    // {
    //     auto kvec_c=irk.first;
    //     fs<<kvec_c<<endl; 
    //     print_complex_matrix_file("Vq_2D",Vq_2D[kvec_c],fs);
    //     print_complex_matrix_file("chi0_2D",chi0_2D[first_freq][kvec_c],fs);
    //     print_complex_matrix_file("pi_2D",pi_freq_2D[first_freq][kvec_c],fs);

    // }
    // fs.close();

    
    complex<double> tot_RPA_energy(0.0,0.0);
    map<Vector3_Order<double>,complex<double>> cRPA_k;
    for(const auto &freq_pair:pi_freq_2D)
    {
        const auto freq=freq_pair.first;
        const auto freq_k=freq_pair.second;
        const double weight=freq_grid[freq];
        for(const auto &k_pair:freq_k)
        {
            const auto kvec_c=k_pair.first;
            const auto mat=k_pair.second;
            complex<double> rpa_for_omega_k(0.0,0.0);
            ComplexMatrix identity(range_all,range_all);
            ComplexMatrix identity_minus_pi(range_all,range_all);
        
            identity.set_as_identity_matrix();
        
            identity_minus_pi=identity-pi_freq_2D[freq][kvec_c];
            complex<double> det_for_rpa(1.0,0.0);
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
            trace_pi=trace(pi_freq_2D.at(freq).at(kvec_c));
            cout<<"PI trace vector:"<<endl;
            if(freq==first_freq)
            {
                for(int ir=0;ir!=pi_freq_2D.at(freq).at(kvec_c).nr;ir++)
                    cout<<pi_freq_2D.at(freq).at(kvec_c)(ir,ir).real()<<"    ";
            }
            cout<<endl;
            rpa_for_omega_k=ln_det+trace_pi; 
            cout<<" freq:"<<freq<<"      rpa_for_omega_k: "<<rpa_for_omega_k<<"      lnt_det: "<<ln_det<<"    trace_pi "<<trace_pi<<endl;
            cRPA_k[kvec_c]+=rpa_for_omega_k*weight*irk_weight[kvec_c]*NSPIN/TWO_PI;
            tot_RPA_energy+=rpa_for_omega_k*weight*irk_weight[kvec_c]*NSPIN/TWO_PI;
        }
    }
    for(auto &ikvec:cRPA_k)
    {
        cout<<ikvec.first<<ikvec.second<<endl;
    }
    cout<<"Grid_num_"<<freq_grid.size()<<"  tot_RPA_energy:  "<<setprecision(8)<<tot_RPA_energy<<endl; 
}

vector<double> Cal_Periodic_Chi0::read_cosine_trans_grid(const string &file_path)
{
    ifstream infile;
    infile.open(file_path);
    
    vector<double> tran;
    string s;
    int n_freq=0;
    //stringstream ss;
    //double ss_d;
    while(getline(infile,s))
    {
        if(s[0]=='F')
        {
            n_freq++;
        }
        else
        {
            stringstream ss(s);
            double ss_d;
            ss>>ss_d;
            tran.push_back(ss_d);
        }
    }    
    infile.close();
    
    double gap;
    gap=get_band_gap();
    cout<<"Cosine_tran_grid"<<endl;
    for(int i=0;i!=tran.size();i++)
    {
        tran[i]/=gap;
        //cout<<tran[i]<<endl;
    }
    return tran;
}

map<double,double> Cal_Periodic_Chi0::read_file_grid(const string &file_path, const char type)
{
    map<double,double> grid;
    ifstream infile;
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
        s0+=pattern;//扩展字符串以方便操作
    
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


void Cal_Periodic_Chi0::print_matrix( char* desc, const matrix &mat )
{
    int nr=mat.nr;
    int nc=mat.nc;
    printf( "\n %s\n", desc );
    for( int i = 0; i < nr; i++ ) 
    {
        for( int j = 0; j < nc; j++ ) printf( " %8.4f", mat.c[i*nc+j] );
            printf( "\n" );
    }
} 
void Cal_Periodic_Chi0::print_complex_matrix( char* desc, const ComplexMatrix &mat )
{
    int nr=mat.nr;
    int nc=mat.nc;
    printf( "\n %s\n", desc );
    for( int i = 0; i < nr; i++ ) 
    {
        for( int j = 0; j < nc; j++ ) printf("%10.6f,%8.6f ",mat.c[i*nc+j].real(),mat.c[i*nc+j].imag());
            printf( "\n" );
    }
} 
void Cal_Periodic_Chi0::print_complex_real_matrix( char* desc, const ComplexMatrix &mat )
{
    int nr=mat.nr;
    int nc=mat.nc;
    printf( "\n %s\n", desc );
    for( int i = 0; i < nr; i++ ) 
    {
        for( int j = 0; j < nc; j++ ) printf( " %10.6f", mat.c[i*nc+j].real());
            printf( "\n" );
    }
} 

void Cal_Periodic_Chi0::print_complex_matrix_file( char* desc, const ComplexMatrix &mat , ofstream &fs)
{
    int nr=mat.nr;
    int nc=mat.nc;
    fs<<desc<<endl;
    for( int i = 0; i < nr; i++ ) 
    {
        for( int j = 0; j < nc; j++ ) fs<<setw(15)<<showpoint<<fixed<<setprecision(8)<<mat.c[i*nc+j].real()<<setw(14)<<showpoint<<fixed<<setprecision(8)<<mat.c[i*nc+j].imag();
        fs<<"\n";
    }
} 


extern "C"
{
    void dgeev_(const char* jobvl,const char* jobvr,const int* n,double *a,
                const int* lda,double* wr,double* wi,double* vl, const int* ldvl, 
                double* vr, const int* ldvr,  double* work,const int* lwork, int* info);
}

void Cal_Periodic_Chi0::dgeev(matrix &mat, vector<complex<double>> &egValue, matrix &vr, int info)
{
    assert(mat.nr==mat.nc);
    vector<double> eg_r(mat.nr);
    vector<double> eg_i(mat.nr);
    matrix vl(mat.nr,mat.nc);
    vr.create(mat.nr,mat.nc);
    int lwork=16*mat.nr;
    vector<double> work(std::max(1,lwork));
    dgeev_("N","V",&mat.nr,mat.c,&mat.nr,&eg_r[0],&eg_i[0],vl.c,&mat.nr,vr.c,&mat.nr,&work[0],&lwork,&info);
    for(int i=0;i!=mat.nr;i++)
        egValue.push_back(complex<double>(eg_r[i],eg_i[i]));
    vr=transpose(vr);
}

void Cal_Periodic_Chi0::rt_cm_max(ComplexMatrix &m)
{
    double rv=m.c[0].real();
    int pos=0;
    for(int i=1;i!=m.size;i++)
    {
        if (m.c[i].real()>rv)
        {
            rv=m.c[i].real();
            pos=i;
        }
    }
    int ir=pos/m.nc;
    int ic=pos%m.nc;
    cout<<"  max_value: "<<rv<<"    position:  ( "<<ir<<"  , "<<ic<<" )"<<endl;
}

void Cal_Periodic_Chi0::rt_m_max(matrix &m)
{
    double rv=m.c[0];
    int pos=0;
    for(int i=1;i!=m.nr*m.nc;i++)
    {
        if (m.c[i]>rv)
        {
            rv=m.c[i];
            pos=i;
        }
    }
    int ir=pos/m.nc;
    int ic=pos%m.nc;
    cout<<"  max_value: "<<rv<<"    position:  ( "<<ir<<"  , "<<ic<<" )"<<endl;
}