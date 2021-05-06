#include "cal_periodic_chi0.h"
#include "lapack_connector.h"
#include <iomanip>
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
    //LOC.wfc_dm_2d.cal_dm();
    cout<<"finish init"<<endl;
	set<Vector3_Order<int>> R_grid(construct_R_grid());
	cout<<"R_grid_num= "<<R_grid.size()<<endl; 
	cout<< minimax_grid_path+"/time_grid.txt"<<endl;
    map<double,double> time_grid(read_file_grid( minimax_grid_path+"/time_grid.txt",'T'));    
     
    map<double,double> freq_grid(read_file_grid( minimax_grid_path+"/freq_grid.txt",'F'));    
    //copy(freq_tmp.begin(),freq_tmp.end(),inserter(freq_grid,freq_grid.begin()));
    //cout<<"freq_grid size: "<<freq_grid.size();
    first_tau=time_grid.begin()->first;
    first_freq=freq_grid.begin()->first;
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
    for(auto &R:R_grid)
    {
        cout<<R<<endl;
        print_complex_matrix("First Green_mat tau",Green_function[0][time_grid.begin()->first][R]); 
        print_complex_matrix("First Green_mat -tau",Green_function[0][-time_grid.begin()->first][R]); 
    }
/*    
    for(int is=0;is!=NSPIN;is++)
    {
        cout<<"  ****TEST Green_function  i_spin: "<<is<<endl;
        for(const auto &time_tau:time_grid)
        {   
            const auto t=time_tau.first;
            cout<<t<<"      "<<Green_function[is][t][{0,0,0}](0,0).real()<<"   "<<Green_function[is][t][{0,0,0}](0,1).real()<<"  "<<Green_function[is][t][{0,0,0}](1,0).real()<<"  "<<Green_function[is][t][{0,0,0}](1,1).real()<<endl;
        }
        cout<<"  ****TEST Green_function[-tau]  i_spin: "<<is<<endl;
        for(const auto &time_tau:time_grid)
        {   
            const auto t=-1*time_tau.first;
            cout<<t<<"      "<<Green_function[is][t][{0,0,0}](0,0).real()<<"   "<<Green_function[is][t][{0,0,0}](0,1).real()<<"  "<<Green_function[is][t][{0,0,0}](1,0).real()<<"  "<<Green_function[is][t][{0,0,0}](1,1).real()<<endl;
        }
        cout<<"finish Green_R_tau"<<endl;
    }
*/
//    READ_AIMS_Cs("FHI-aims-out/Cs_data.txt");
    // for(auto &Ip:Cs)
    //     for(auto &Jp:Ip.second)
    //         for(auto &Rp:Jp.second)
    //         {
    //             auto &cs_mat=*Rp.second;
    //             print_matrix("Cs_mat",cs_mat);
    //         }
   // READ_AIMS_Vq("FHI-aims-out/Coulomb_q_mat");
   /*  for(auto &Ip:Vq)
        for(auto &Jp:Ip.second)
            for(auto &qp:Jp.second)
            {
                auto &vq_mat=*qp.second;
                print_complex_matrix("Vq_mat ", vq_mat);
            } */
   
    
   // cal_local_Cs();
    //cout<<"finish cal_local_Cs"<<endl;
    
    for(const auto &time_tau:time_grid)
        for(const Vector3_Order<int> &R:R_grid)
            cal_chi0_element(time_tau.first,R);
    
    for(const Vector3_Order<int> &R:R_grid)
        print_complex_matrix("chi0_tau",chi0[time_grid.begin()->first][R][0][0]);
 /*   cout<<"  ***TEST chi0_tau"<<endl;
    cout<<"  mu_num:"<<chi0[time_grid.begin()->first][{0,0,0}][0][0].nr<<"    nu_num:"<<chi0[time_grid.begin()->first][{0,0,0}][0][0].nc<<endl;
    for(const auto &time_tau:time_grid)
        cout<<time_tau.first<<"     "<<chi0[time_tau.first][{0,0,0}][0][0](0,0)<<"    "<<chi0[time_tau.first][{0,0,0}][0][0](0,1)<<"    "<<chi0[time_tau.first][{0,0,0}][0][0](1,0) \
            <<"    "<<chi0[time_tau.first][{0,0,0}][0][0](1,1)<<"    "<<chi0[time_tau.first][{0,0,0}][0][0](2,2)<<endl;
    cout<<"finish cal_chi0_element"<<endl;
*/
    Cosine_to_chi0_freq(R_grid,time_grid,freq_grid);
    cout<<"finish FT_to_chi0_freq"<<endl;
    
/*    cout<<"  ***TEST chi0_for_freq:"<<endl;
    for(auto &freq_pair:chi0_freq)
        cout<<dec<<freq_pair.first<<"     "<<chi0_freq[freq_pair.first][{0,0,0}][0][0](0,0)<<"     "<<chi0_freq[freq_pair.first][{0,0,0}][0][0](0,1) \
            <<"     "<<chi0_freq[freq_pair.first][{0,0,0}][0][0](1,0)<<"     "<<chi0_freq[freq_pair.first][{0,0,0}][0][0](1,1)<<"     "<<chi0_freq[freq_pair.first][{0,0,0}][0][0](2,2)<<endl;
*/
    for(auto &R:R_grid)
    {
        cout<<R<<endl;
        print_complex_matrix("chi0_freq_mat",chi0_freq[freq_grid.begin()->first][R][0][0]);
    }
/*
    for(auto &freq_pair:chi0_freq)
    {
        auto freq=freq_pair.first;
        for(auto &box_pair:freq_pair.second)
        {
            auto box=box_pair.first;
            for(auto &I_pair:box_pair.second)
            {
                auto I=I_pair.first;
                for(auto &J_pair:I_pair.second)
                {
                    auto J=J_pair.first;
                    auto &mat=J_pair.second;
                    cout<<freq<<"   "<<box<<"   "<<I<<"   "<<J<<endl;
                    print_complex_matrix("chi0_freq mat", mat);
                }
            }
        }
    }  
    //cal_pi_freq();
    //cout<<"finish cal_pi_freq"<<endl;
    // for(auto &freq_pair:pi_freq)
        // cout<<dec<<freq_pair.first<<"     "<<pi_freq[freq_pair.first][{0,0,0}][0][0](0,0)<<"     "<<pi_freq[freq_pair.first][{0,0,0}][0][0](0,1)<<"     "<<pi_freq[freq_pair.first][{0,0,0}][0][0](1,0)<<"     "<<pi_freq[freq_pair.first][{0,0,0}][0][0](1,1)<<"     "<<pi_freq[freq_pair.first][{0,0,0}][0][0](2,2)<<endl;
    //cal_MP2_energy( freq_grid );
    //cout<<"finish cal_MP2_energy"<<endl;

    cal_pi_k();
    */
    FT_R_to_K_chi0();
    for(int i=0;i!=n_kpoints;i++)
    {
        cout<<kvec_c[i];
        print_complex_matrix("chi0_k",chi0_k[freq_grid.begin()->first][kvec_c[i]][0][0]);
    }
    
    cal_pi_k_use_aims_vq();
    /*
    for(int i=0;i!=n_kpoints;i++)
    {
        cout<<kvec_c[i];
        print_complex_matrix("PI mat",pi_k[freq_grid.begin()->first][kvec_c[i]][0][0]);
    }
    */
    for(auto &freq_pair:pi_k)
        cout<<dec<<freq_pair.first<<"     "<<pi_k[freq_pair.first][{0,0,0}][0][0](0,0)<<"     "<<pi_k[freq_pair.first][{0,0,0}][0][0](0,1)<<"     "<<pi_k[freq_pair.first][{0,0,0}][0][0](1,0)<<"     "<<pi_k[freq_pair.first][{0,0,0}][0][0](1,1)<<"     "<<pi_k[freq_pair.first][{0,0,0}][0][0](2,2)<<endl;
    
    cal_MP2_energy_pi_k(freq_grid);
    //cal_MP2_energy_pi_k_tau(time_grid,cal_pi_k_tau());
    RPA_correlation_energy(freq_grid);
   // cout<<" periodic_exx_lcao.enegy()"<<exx_lcao.energy<<endl;
    //cout<<" periodic_en.etot"<<en.etot<<endl;

    
}
/*
void Cal_Periodic_Chi0::adj_atom_search()
{
    cout<<"begin adj_atom_search"<<endl;
    //adj_atom_pair.clear();
 
    for(const auto &I_pair:exx_lcao.Cps)
    {
        const auto I=I_pair.first;
        const auto Cs_I=I_pair.second;
        cout<<"  adj_2"<<endl;
        for(const auto &J_pair:Cs_I)
        {
            cout<<"   adj_2.5"<<endl;
            const auto J=J_pair.first;
            adj_atom_pair[I][J].create(1,1);
            cout<<"   adj_2.6"<<endl;
            if(I<J)
            {  
                adj_atom_pair[J][I].create(1,1);
            }
        }
    }
    
    cout<<"  adj_atom_pair"<<endl;
    for(const auto &I_pair:adj_atom_pair)
    {   
        const auto I=I_pair.first;
        const auto I_set=I_pair.second;
        cout<<I<<"  :  ";
        for(const auto &J_pair:I_set)
            cout<<"   "<<J_pair.first;
        cout<<endl;
    }
}
*/
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
    temp_Cs_count=0;
    //cout<<" knpoints_num= "<<n_kpoints<<endl;
    green_k.resize(n_kpoints);
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



void Cal_Periodic_Chi0::cal_Green_func_element_k( const matrix &wg, const double time_tau )
{
	//cout<<wg.nr<<"   "<<wg.nc<<"  "<<time_tau<<endl;
    // dm = wfc.T * wg * wfc.conj()
	// dm[ik](iw1,iw2) = \sum_{ib} wfc[ik](ib,iw1).T * wg(ik,ib) * wfc[ik](ib,iw2).conj()
    //assert(wg.nc<=NLOCAL);
	//assert(wg.nr==n_kpoints);
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
}

void Cal_Periodic_Chi0::cal_Green_func_R_tau( const double &time_tau, const Vector3_Order<int> &R,  const Vector3<double>* kvec_c)
{
    Vector3_Order<double> R_2(R.x,R.y,R.z);
   
    Vector3_Order<int> R_0(0,0,0);
    //cout<<"R_0 init"<<endl;
    //cout<<"R= "<<R<<"     time_tau= "<<time_tau<<endl;
    
    //cout<<"   zero_out";
    //cout<<"Green_function_zero_out"<<endl;
    const double OCC_CUT=0.1;
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
        cal_Green_func_element_k(wg_unocc, time_tau);
    } 
    else
    {
        cal_Green_func_element_k( wg, time_tau);
    }
	
    //cout<<"finish cal_Green_func_element_k"<<endl;
    //cout<<"green_k[0].nr="<<green_k[0].nr<<"   green_k[0].nc="<<green_k[0].nc<<endl;
    
    //cout<<"Green_function_size="<<Green_function[time_tau][R].size<<endl;
    for(int is=0;is!=NSPIN;is++)
    {
        Green_function[is][time_tau][R].create(green_k[0].nr,green_k[0].nc);
        Green_function[is][time_tau][R].zero_out();
    }
    //const complex<double> n_kpoints_Rec(1.0/n_kpoints,0.0);
    //cout<<"  n_kpoints_Rec= "<<n_kpoints_Rec<<endl;
    int iter_nks=0;
    for(int is=0;is!=NSPIN;is++)
    {
        for(int ik=0; ik!=n_kpoints/NSPIN; ++ik)
        {
            const double arg =-1*( kvec_c[ik+iter_nks]* (R_2*latvec)) * TWO_PI; //latvec
            const complex<double> kphase = complex <double> ( cos(arg),  sin(arg) );
            Green_function[is][time_tau][R] += green_k[ik+iter_nks] * kphase ;
            // if(time_tau==first_tau)
            // {
            //     cout<<" IN cal_Green   R, kvec_d   kphase: "<<R_2<<kvec_c[ik+iter_nks]<<kphase<<endl;
            //     print_complex_matrix(" green_k_mat",green_k[ik+iter_nks]);
            //     print_complex_matrix(" green_k_mat_plus_kphase",green_k[ik+iter_nks]*kphase);
            // }
            // for(int nr=0;nr!=Green_function[is][time_tau][R].nr;nr++)
            //     for(int nc=0;nc!=Green_function[is][time_tau][R].nc;nc++)
            //     {
            //         if(nc<nr) Green_function[is][time_tau][R](nr,nc)=Green_function[is][time_tau][R](nc,nr);
            //     }
        }

        if(time_tau<=0){Green_function[is][time_tau][R]*=-1;}
        iter_nks+=n_kpoints/NSPIN;
    }
}


/*
void Cal_Periodic_Chi0::FT_to_chi0_freq(const set<Abfs::Vector3_Order<int>> &R_grid, map<double,double> &time_grid, map<double,double> &freq_grid )
{
    //map<double,double> freq_grid(construct_gauss_grid(120));  

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
                    for(const auto &freq_pair:freq_grid)
                    {
                        const double freq=freq_pair.first;
                        const double freq_weight=freq_pair.second;
                        const double time_weight=time_grid[time];
                        const double cos_weight=2*cos(freq*time)*time_weight;
                        chi0_freq[freq][R][I][J]+=mat*cos_weight;  
                        //cout<<time<<"   "<<time_weight<<"    freq:"<<freq<<"   cos_weight"<<cos_weight  \
                            <<"      mat:  "<<mat(0,0)<<"    "<<chi0_freq[freq][R][I][J](0,0)<<endl;
                        
                    }
                }
            }
        } 
    }
    
    
    // for(const auto &freq_pair:freq_grid)
    // {
        // const double freq=freq_pair.first;
        // for(const auto &R:R_grid)
        // {
            // for(const auto &I_pair:chi0[(*time_grid.begin()).first][R])
            // {
                // const size_t I=I_pair.first;
                // const auto J_map=I_pair.second;
                // const size_t mu_num=exx_lcao.index_abfs[ucell.iat2it[I]].count_size;
                // for(const auto &J_pair:J_map)
                // {
                    // const size_t J=J_pair.first;
                    // const auto mat=J_pair.second;
                    // const size_t nu_num=exx_lcao.index_abfs[ucell.iat2it[J]].count_size;
                    // chi0_freq[freq][R][I][J].create(mu_num,nu_num);
                    // for(const auto &time:time_grid)
                    // {
                        // double cos_weight=2*std::cos(freq*time.first)*time.second;
                        // chi0_freq[freq][R][I][J]+=mat*cos_weight;
                    // }
                // }
            // }
        // }
    // }
} 
*/

void Cal_Periodic_Chi0::Cosine_to_chi0_freq(const set<Vector3_Order<int>> &R_grid, map<double,double> &time_grid, map<double,double> &freq_grid )
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
                        chi0_freq[freq][R][I][J]+=mat*cos_weight;  
                        //cout<<time<<"   "<<time_weight<<"    freq:"<<freq<<"   cos_weight"<<cos_weight  \
                            <<"      mat:  "<<mat(0,0)<<"    "<<chi0_freq[freq][R][I][J](0,0)<<endl;
                        iter_freq++;
                    }
                    
                }
            }
        }
        iter_tau++;
    }
    
}
void Cal_Periodic_Chi0::cal_chi0_element_origin(const double &time_tau, const Vector3_Order<int> &R)
{
    //cout<<"In cal_chi0_element::time_tau="<<time_tau<<"     R="<<R<<"    Green_function"<<Green_function[time_tau][R].size<<endl;
    //complex<double> *chi0_ptr=chi0[time_tau][R].c;
    //size_t i_shift=0;
    Vector3_Order<int> R_0(0,0,0);
    for(auto &I_pair:Cs)
    {
        const size_t I_index=I_pair.first;
        const size_t i_num=atom_nw[I_index];
        const size_t mu_num=atom_mu[I_index];
        for(auto &J_pair:Cs)
        {
            const size_t J_index=J_pair.first;
            const size_t j_num=atom_nw[J_index];
            const size_t nu_num=atom_mu[J_index];

            chi0[time_tau][R][I_index][J_index].create(mu_num,nu_num);

            ComplexMatrix O_sum(mu_num,nu_num);

            for(const auto &K_pair:I_pair.second)
            {
                const auto K_index=K_pair.first;
                const size_t k_num=atom_nw[K_index];
                
                for(const auto &R1_index:K_pair.second)
                {
                    const auto R1=R1_index.first; 
                    const auto &Cs_mat1=R1_index.second;
                    
                    for(int is=0;is!=NSPIN;is++)
                    {
                        ComplexMatrix N_R2(k_num,j_num*nu_num);
                        ComplexMatrix N_conj_R2(k_num,j_num*nu_num);
                        
                        ComplexMatrix X_R2(i_num,j_num*nu_num);
                        ComplexMatrix X_conj_R2(i_num,j_num*nu_num);
                    
                        for(const auto &L_pair:J_pair.second)
                        {
                            const auto L_index=L_pair.first;
                            const size_t l_num=atom_nw[L_index];
                            
                            for(const auto &R2_index:L_pair.second)
                            {
                                const auto R2=R2_index.first; 
                                const auto &Cs_mat2=R2_index.second;
                                // if(temp_Cs_count == 0)
                                // {
                                //     cout<<"R IJKL:"<<R<<"   "<<I_index<<J_index<<K_index<<L_index<<"    R1,R2:"<<R1<<R2<<endl;
                                // }
                               
                                ComplexMatrix Green_kl(k_num,l_num);
                                ComplexMatrix Green_kl_conj(k_num,l_num);
                                //cout<<"R2_debug_0"<<endl;
                                Vector3_Order<int> R_temp_1(R1-R2-R);
                            
                                generate_Green_part(is,time_tau,R_temp_1%R_periodic,K_index,L_index,Green_kl,Green_kl_conj);
                            
                                ComplexMatrix Green_il(i_num,l_num);
                                ComplexMatrix Green_il_conj(i_num,l_num); 
                            
                                Vector3_Order<int> R_temp_2(R2+R); 
                                generate_Green_part(is,-1*time_tau,R_temp_2%R_periodic,I_index,L_index,Green_il,Green_il_conj);
                            
                                assert(j_num*l_num==(*Cs_mat2).nr); 
                                
                                ComplexMatrix Cs2_reshape(reshape_Cs(j_num,l_num,nu_num,Cs_mat2));
                            
                                LapackConnector::gemm(
                                    'N','N',
                                    k_num,j_num*nu_num,l_num,
                                    1,
                                    Green_kl.c,l_num,
                                    Cs2_reshape.c,j_num*nu_num,
                                    1,
                                    N_R2.c,
                                    j_num*nu_num);
                                    
                                LapackConnector::gemm(
                                    'N','N',
                                    k_num,j_num*nu_num,l_num,
                                    1,
                                    Green_kl_conj.c,l_num,
                                    Cs2_reshape.c,j_num*nu_num,
                                    1,
                                    N_conj_R2.c,
                                    j_num*nu_num);
                                
                                LapackConnector::gemm(
                                    'N','N',
                                    i_num,j_num*nu_num,l_num,
                                    1,
                                    Green_il.c,l_num,
                                    Cs2_reshape.c,j_num*nu_num,
                                    1,
                                    X_R2.c,
                                    j_num*nu_num);
                                    
                                LapackConnector::gemm(
                                    'N','N',
                                    i_num,j_num*nu_num,l_num,
                                    1,
                                    Green_il_conj.c,l_num,
                                    Cs2_reshape.c,j_num*nu_num,
                                    1,
                                    X_conj_R2.c,
                                    j_num*nu_num);
                            }
                        }
                        
                        ComplexMatrix N_R2_rs(reshape_complexmat(k_num,j_num,nu_num,N_R2));
                        ComplexMatrix N_conj_R2_rs(reshape_complexmat(k_num,j_num,nu_num,N_conj_R2));
                        ComplexMatrix X_R2_rs(reshape_complexmat(i_num,j_num,nu_num,X_R2));
                        ComplexMatrix X_conj_R2_rs(reshape_complexmat(i_num,j_num,nu_num,X_conj_R2));
                        
                        ComplexMatrix Green_ij(i_num,j_num);
                        ComplexMatrix Green_ij_conj(i_num,j_num);
                        generate_Green_part(is,-1*time_tau,R,I_index,J_index,Green_ij,Green_ij_conj);
                        
                        ComplexMatrix Green_kj(k_num,j_num);
                        ComplexMatrix Green_kj_conj(k_num,j_num);
                        Vector3_Order<int> R_temp_3(R1-R);
                        generate_Green_part(is,time_tau,R_temp_3%R_periodic,K_index,J_index,Green_kj,Green_kj_conj);
                    
                        
                        //cal and sum O
                        ComplexMatrix O(i_num*k_num,nu_num);
                        ComplexMatrix Z(k_num*i_num,nu_num);
                        
                        LapackConnector::gemm(
                                'N','N',
                                i_num,k_num*nu_num,j_num,
                                1,
                                Green_ij.c,j_num,
                                N_R2_rs.c,k_num*nu_num,
                                1,
                                O.c,
                                k_num*nu_num);
                        
                        LapackConnector::gemm(
                                'N','N',
                                i_num,k_num*nu_num,j_num,
                                1,
                                Green_ij_conj.c,j_num,
                                N_conj_R2_rs.c,k_num*nu_num,
                                1,
                                O.c,
                                k_num*nu_num);
                        
                        LapackConnector::gemm(
                                'N','N',
                                k_num,i_num*nu_num,j_num,
                                1,
                                Green_kj.c,j_num,
                                X_R2_rs.c,i_num*nu_num,
                                1,
                                Z.c,
                                i_num*nu_num);
                        
                        LapackConnector::gemm(
                                'N','N',
                                k_num,i_num*nu_num,j_num,
                                1,
                                Green_kj_conj.c,j_num,
                                X_conj_R2_rs.c,i_num*nu_num,
                                1,
                                Z.c,
                                i_num*nu_num);
                        
                        ComplexMatrix Z_rs(reshape_complexmat(k_num,i_num,nu_num,Z));
                        
                        O+=Z_rs;
                        
                        ComplexMatrix Cs1_tran(transpose(*Cs_mat1));
                        
                        LapackConnector::gemm(
                                'N','N',
                                mu_num,nu_num,i_num*k_num,
                                1,
                                Cs1_tran.c,i_num*k_num,
                                O.c,nu_num,
                                1,
                                O_sum.c,
                                nu_num);
                            //l_shift+=l_num;
                            //cout<<"   the_end  l_shift"<<l_shift<<endl;
                    }
                }
            }
            //j_shift+=j_num;
            chi0[time_tau][R][I_index][J_index]=O_sum;
            
            //cblas_zcopy(chi0_sum_KL.size, chi0_sum_KL.c, 1, chi0_ptr, 1);
            //cout<<"chi0_sum_KL.c[0]:"<<chi0_sum_KL.c[0]<<"     chi0.c[ptr]:"<<*chi0_ptr<<endl;
            //chi0_ptr+=chi0_sum_KL.size;
            //joint to big mat 
        }
        //i_shift+=i_num;
    }
    ++temp_Cs_count;
} 

/*

map<double,double> Cal_Periodic_Chi0::construct_gauss_grid(const int Npoints)
{
    // int npoints=5;
    // double dis_points=0.1;  
    //set<double> time_tau_grid;
    // for(int i=0;i!=npoints;++i)
        // time_tau_grid.insert(-i*dis_points);
    //Npoints=40;
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

void Cal_Periodic_Chi0::cal_local_Cs()
{   
    
    //{i,j}-blocks -> {I,J}-atoms 
    std::vector<int> iwt_row(ParaO.nrow,0);
    cout<<"cal_Cs 1"<<endl;
    for(int iw_global=0;iw_global!=green_k[0].nr;++iw_global)
    {
        const int iw_local = ParaO.trace_loc_row[iw_global];
        if(iw_local>=0)
            iwt_row[iw_local]=iw_global;
    }
    
    iwt2iat_row.clear();
    for(int i=0;i!=iwt_row.size();++i)
    {
        iwt2iat_row.insert(ucell.iwt2iat[iwt_row[i]]);
    }
    cout<<"cal_Cs 2"<<endl;
    cout<<"ParaO.ncol="<<ParaO.ncol<<endl;
    std::vector<int> iwt_col(ParaO.ncol,0);
    cout<<"cal_Cs 2.5"<<endl;
    for(int iw_global=0;iw_global!=green_k[0].nc;++iw_global)
    {
        const int iw_local = ParaO.trace_loc_col[iw_global];
        if(iw_local>=0)
            iwt_col[iw_local]=iw_global;
    }
    
    cout<<"cal_Cs 2.6"<<endl;
    //cout<<"iwt2iat_col.size="<<iwt_col.size()<<endl;
    iwt2iat_col.clear();
    cout<<iwt2iat_col.size()<<endl;
    for(int i=0;i!=iwt_col.size();++i)
    {
        //cout<<i<<"   "<<ucell.iwt2iat[iwt_col[i]]<<endl;
        iwt2iat_col.insert(ucell.iwt2iat[iwt_col[i]]);
    }
    // vector<int> m(iwt_col.size());
    // for(int i=0;i!=iwt_col.size();++i)
    // {
        // m[i]=ucell.iwt2iat[iwt_col[i]];
    // }
    // cout<<"cal_Cs 2.7"<<endl;
    // iwt2iat_col.insert(std::begin(m),std::end(m));    
    //iwt2iat_col={0,1};   
    all_mu=0;
    for(set<size_t>::iterator I_index=iwt2iat_row.begin();I_index!=iwt2iat_row.end();++I_index)
    {
        //const size_t i_num=ucell.atoms[ucell.iat2it[*I_index]].nw;
        const size_t mu_num=exx_lcao.index_abfs[ucell.iat2it[*I_index]].count_size;
        all_mu+=mu_num;
    }
    
    all_nu=0;
    for(set<size_t>::iterator J_index=iwt2iat_col.begin();J_index!=iwt2iat_col.end();++J_index)
    {
        //const size_t j_num=ucell.atoms[ucell.iat2it[*J_index]].nw;
        const size_t nu_num=exx_lcao.index_abfs[ucell.iat2it[*J_index]].count_size;
        all_nu+=nu_num;
    }
    
    std::set<size_t> iat_row_col(iwt2iat_row.begin(),iwt2iat_row.end());
    iat_row_col.insert(iwt2iat_col.begin(),iwt2iat_col.end());
    cout<<"cal_Cs 4"<<endl;
    // vector<map<size_t,vector<Abfs::Vector3_Order<int>>>> adjs_row(iwt2iat_row.size());
    // int vec_iter=0;
    // for(set<size_t>::iterator ad=iwt2iat_row.begin();ad!=iwt2iat_row.end();++ad)
    // {
        // adjs_row[vec_iter]=Abfs::get_adjs( *ad );
        // vec_iter++;
        // cout<<"vec_iter: "<<vec_iter<<"    ad: "<<*ad<<endl;
    // }
    
    cout<<"loop1_3"<<endl;
    cout<<"iwt2iat_vec.size()"<<iwt2iat_row.size()<<endl;
    
    //Pre-Cal Cs_m
    map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,shared_ptr<matrix>>>> Cs_temp;
    Cs_temp = Abfs::cal_Cs( iat_row_col, exx_lcao.m_abfs_abfs,exx_lcao.m_abfslcaos_lcaos, exx_lcao.index_abfs,exx_lcao.index_lcaos, exx_lcao.info.c_threshold, exx_lcao.Cws,exx_lcao.Vws );
    Cs_m=Abfs::cal_mps( R_periodic, Cs_temp );
    cout<<"loop1_4"<<endl;
    //cout<<"Cs_m first term: "<<(*Cs_m[0][0][{0,0,0}])(0,0)<<endl; 
    //print_matrix(" (*Cs_m[0][0][{0,0,0}]) ",(*Cs_m[0][0][{0,0,0}]).nc,(*Cs_m[0][0][{0,0,0}]).nr,(*Cs_m[0][0][{0,0,0}]).c,(*Cs_m[0][0][{0,0,0}]).nc);
    for(auto I=iat_row_col.begin();I!=iat_row_col.end();++I)
    {
        cout<<"  I: "<<*I;
        auto adj_I=Abfs::get_adjs(*I);
        for(const auto &J_pair:adj_I)
        {
            cout<<"    "<<J_pair.first;
        }
        cout<<endl;
    }
}

*/
void Cal_Periodic_Chi0::cal_chi0_element(const double &time_tau, const Vector3_Order<int> &R)
{
    //cout<<"In cal_chi0_element::time_tau="<<time_tau<<"     R="<<R<<"    Green_function"<<Green_function[time_tau][R].size<<endl;
    //complex<double> *chi0_ptr=chi0[time_tau][R].c;
    //size_t i_shift=0;
    Vector3_Order<int> R_0(0,0,0);
    for(auto &I_pair:Cs)
    {
        const size_t I_index=I_pair.first;
        const size_t i_num=atom_nw[I_index];
        const size_t mu_num=atom_mu[I_index];
        for(auto &J_pair:Cs)
        {
            const size_t J_index=J_pair.first;
            const size_t j_num=atom_nw[J_index];
            const size_t nu_num=atom_mu[J_index];

            chi0[time_tau][R][I_index][J_index].create(mu_num,nu_num);

            ComplexMatrix O_sum(mu_num,nu_num);

            for(const auto &K_pair:I_pair.second)
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
                        ComplexMatrix N_R2(k_num,j_num*nu_num);
                        ComplexMatrix N_conj_R2(k_num,j_num*nu_num);
                        
                        ComplexMatrix X_R2(i_num,j_num*nu_num);
                        ComplexMatrix X_conj_R2(i_num,j_num*nu_num);
                    
                        for(const auto &L_pair:J_pair.second)
                        {
                            const auto L_index=L_pair.first;
                            const size_t l_num=atom_nw[L_index];
                            
                            for(const auto &R2_index:L_pair.second)
                            {
                                const auto R2=R2_index.first; 
                                const auto &Cs_mat2=R2_index.second;
                                if(time_tau==first_tau && R==Vector3_Order<int>(0,0,0) && R1==Vector3_Order<int>(0,0,0))
                                {
                                    cout<<"R IJKL:"<<R<<"   "<<I_index<<J_index<<K_index<<L_index<<"    R1,R2:"<<R1<<R2<<endl;
                                }
                               
                                ComplexMatrix Green_kl(k_num,l_num);
                                ComplexMatrix Green_kl_conj(k_num,l_num);
                                //cout<<"R2_debug_0"<<endl;
                               // Vector3_Order<int> R_temp_1(R1-R2-R);
                                Vector3_Order<int> R_temp_1(R+R2-R1);
                            
                                generate_Green_part(is,time_tau,R_temp_1%R_periodic,K_index,L_index,Green_kl,Green_kl_conj);//attention R
                            
                                ComplexMatrix Green_il(i_num,l_num);
                                ComplexMatrix Green_il_conj(i_num,l_num); 
                            
                                Vector3_Order<int> R_temp_2(R2+R); 
                                generate_Green_part(is,time_tau,R_temp_2%R_periodic,I_index,L_index,Green_il,Green_il_conj);
                            
                                assert(j_num*l_num==(*Cs_mat2).nr); 
                                
                                matrix Cs2_reshape_temp(reshape_Cs(j_num,l_num,nu_num,Cs_mat2));
                                ComplexMatrix Cs2_reshape(reshape_dim_Cs(l_num,j_num,mu_num,make_shared<matrix>(Cs2_reshape_temp)));
                              //  Green_kl=transpose(Green_kl,0);
                               // Green_kl_conj=transpose(Green_kl_conj,0);    
                                //Green_il_conj=transpose(Green_il_conj,0);
                                N_R2+=Green_kl*Cs2_reshape;

                                N_conj_R2+=Green_kl_conj*Cs2_reshape;

                                X_R2+=Green_il*Cs2_reshape;

                                X_conj_R2+=Green_il_conj*Cs2_reshape;
                                // LapackConnector::gemm(
                                //     'N','N',
                                //     k_num,j_num*nu_num,l_num,
                                //     1,
                                //     Green_kl.c,l_num,
                                //     Cs2_reshape.c,j_num*nu_num,
                                //     1,
                                //     N_R2.c,
                                //     j_num*nu_num);
                                    
                                // LapackConnector::gemm(
                                //     'N','N',
                                //     k_num,j_num*nu_num,l_num,
                                //     1,
                                //     Green_kl_conj.c,l_num,
                                //     Cs2_reshape.c,j_num*nu_num,
                                //     1,
                                //     N_conj_R2.c,
                                //     j_num*nu_num);
                                
                                // LapackConnector::gemm(
                                //     'N','N',
                                //     i_num,j_num*nu_num,l_num,
                                //     1,
                                //     Green_il.c,l_num,
                                //     Cs2_reshape.c,j_num*nu_num,
                                //     1,
                                //     X_R2.c,
                                //     j_num*nu_num);
                                    
                                // LapackConnector::gemm(
                                //     'N','N',
                                //     i_num,j_num*nu_num,l_num,
                                //     1,
                                //     Green_il_conj.c,l_num,
                                //     Cs2_reshape.c,j_num*nu_num,
                                //     1,
                                //     X_conj_R2.c,
                                //     j_num*nu_num);
                                // if(time_tau==first_tau && R==Vector3_Order<int>(0,0,0) && R1==Vector3_Order<int>(0,0,0))
                                // {
                                //     print_complex_matrix("Green_kl",Green_kl);
                                //     print_complex_matrix("Green_il_conj",Green_il_conj);
                                //     print_complex_matrix("Green_kl_conj",Green_kl_conj);
                                //     print_complex_matrix("Green_il",Green_il);
                                //    // print_complex_matrix("Cs2_reshape",Cs2_reshape);
                                //    // print_complex_matrix("N_conj_R2",N_R2);
                                //     //print_complex_matrix("X_R2",X_conj_R2);
                                // }
                            }
                        }
                        
                        ComplexMatrix N_R2_rs(reshape_complexmat(k_num,j_num,nu_num,N_R2));
                        ComplexMatrix N_conj_R2_rs(reshape_complexmat(k_num,j_num,nu_num,N_conj_R2));
                        ComplexMatrix X_R2_rs(reshape_complexmat(i_num,j_num,nu_num,X_R2));
                        ComplexMatrix X_conj_R2_rs(reshape_complexmat(i_num,j_num,nu_num,X_conj_R2));
                        
                        ComplexMatrix Green_ij(i_num,j_num);
                        ComplexMatrix Green_ij_conj(i_num,j_num);
                        generate_Green_part(is,time_tau,R,I_index,J_index,Green_ij,Green_ij_conj);
                        
                        ComplexMatrix Green_kj(k_num,j_num);
                        ComplexMatrix Green_kj_conj(k_num,j_num);
                        Vector3_Order<int> R_temp_3(R-R1);
                        generate_Green_part(is,time_tau,R_temp_3%R_periodic,K_index,J_index,Green_kj,Green_kj_conj);
                    
                        
                        //cal and sum O
                        ComplexMatrix O(i_num,k_num*nu_num);
                        ComplexMatrix Z(k_num,i_num*nu_num);
                        
                        O+=Green_ij_conj*N_R2_rs;
                        O+=Green_ij*N_conj_R2_rs;

                        Z+=Green_kj_conj*X_R2_rs;
                        Z+=Green_kj*X_conj_R2_rs;
                        
                        // LapackConnector::gemm(
                        //         'N','N',
                        //         i_num,k_num*nu_num,j_num,
                        //         1,
                        //         Green_ij.c,j_num,
                        //         N_R2_rs.c,k_num*nu_num,
                        //         1,
                        //         O.c,
                        //         k_num*nu_num);
                        // if(time_tau==first_tau && R==Vector3_Order<int>(0,0,0) && R1==Vector3_Order<int>(0,0,0))
                        // {
                        //     print_complex_matrix("Green_ij",Green_ij);
                        //     print_complex_matrix("Green_ij_conj",Green_ij_conj);
                        //     print_complex_matrix("N_R2_rs",N_R2_rs);
                        //     print_complex_matrix("N_conj_R2_rs",N_conj_R2_rs);
                        //     print_complex_matrix("O_first",O);
                        // }
                        // LapackConnector::gemm(
                        //         'N','N',
                        //         i_num,k_num*nu_num,j_num,
                        //         1,
                        //         Green_ij_conj.c,j_num,
                        //         N_conj_R2_rs.c,k_num*nu_num,
                        //         1,
                        //         O.c,
                        //         k_num*nu_num);
                        
                        // LapackConnector::gemm(
                        //         'N','N',
                        //         k_num,i_num*nu_num,j_num,
                        //         1,
                        //         Green_kj.c,j_num,
                        //         X_R2_rs.c,i_num*nu_num,
                        //         1,
                        //         Z.c,
                        //         i_num*nu_num);
                        // if(time_tau==first_tau && R==Vector3_Order<int>(0,0,0) && R1==Vector3_Order<int>(0,0,0))
                        // {
                        //     // print_complex_matrix("Green_kj",Green_kj);
                        //     // print_complex_matrix("Green_kj_conj",Green_kj_conj);
                        //     // print_complex_matrix("X_R2_rs",X_R2_rs);
                        //     // print_complex_matrix("X_conj_R2_rs",X_conj_R2_rs);
                        //     // print_complex_matrix("Z_first",Z);
                        // }
                        // LapackConnector::gemm(
                        //         'N','N',
                        //         k_num,i_num*nu_num,j_num,
                        //         1,
                        //         Green_kj_conj.c,j_num,
                        //         X_conj_R2_rs.c,i_num*nu_num,
                        //         1,
                        //         Z.c,
                        //         i_num*nu_num);
                        // if(time_tau==first_tau && R==Vector3_Order<int>(0,0,0))
                        // {
                        //     print_complex_matrix("O_mat_temp",O);
                        //     print_complex_matrix("Z_mat_temp",Z);
                        // }
                        ComplexMatrix Z_rs(reshape_complexmat(k_num,i_num,nu_num,Z));
                        
                        O+=Z_rs;
                        ComplexMatrix OZ(reshape_complexmat_21(i_num,k_num,nu_num,O));
                        ComplexMatrix Cs1_tran(transpose(*Cs_mat1));
                        O_sum+=Cs1_tran*OZ;
                        // LapackConnector::gemm(
                        //         'N','N',
                        //         mu_num,nu_num,i_num*k_num,
                        //         1,
                        //         Cs1_tran.c,i_num*k_num,
                        //         O.c,nu_num,
                        //         1,
                        //         O_sum.c,
                        //         nu_num);
                        // if((time_tau==-first_tau || time_tau==first_tau) && R==Vector3_Order<int>(0,0,0) )
                        // {
                        //     cout<<"R1:  "<<R1<<endl;
                        //     print_complex_matrix("Cs1_tran",Cs1_tran);
                        //     print_complex_matrix("O_last",O);
                        //     print_complex_matrix("O_sum_temp",Cs1_tran*O);
                        // }
                            //l_shift+=l_num;
                            //cout<<"   the_end  l_shift"<<l_shift<<endl;
                    }
                }
            }
            //j_shift+=j_num;
            chi0[time_tau][R][I_index][J_index]=O_sum;
            
            //cblas_zcopy(chi0_sum_KL.size, chi0_sum_KL.c, 1, chi0_ptr, 1);
            //cout<<"chi0_sum_KL.c[0]:"<<chi0_sum_KL.c[0]<<"     chi0.c[ptr]:"<<*chi0_ptr<<endl;
            //chi0_ptr+=chi0_sum_KL.size;
            //joint to big mat 
        }
        //i_shift+=i_num;
    }
    ++temp_Cs_count;
} 

matrix Cal_Periodic_Chi0::reshape_Cs(const size_t n1, const size_t n2, const size_t n3, const shared_ptr<matrix> &Cs) //(n1,n2,n3) -> (n2,n1,n3)
{
    const auto length = sizeof(double)*n3;
	const auto n13 = n1*n3;
	const double *m_ptr = (*Cs).c;
	matrix m_new( n1*n2, n3, false );
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


ComplexMatrix Cal_Periodic_Chi0::reshape_complexmat(const size_t n1, const size_t n2, const size_t n3, const ComplexMatrix &mat) //(n1,n2*n3) -> (n2,n1*n3)
{
    const auto length = sizeof(complex<double>)*n3;
	const auto n13 = n1*n3;
	const complex<double> *m_ptr = mat.c;
	ComplexMatrix m_new( n2, n1*n3, false );
	for( size_t i1=0; i1!=n1; ++i1 )
	{
		complex<double> *m_new_ptr = m_new.c+i1*n3;
		for( size_t i2=0; i2!=n2; ++i2, m_ptr+=n3, m_new_ptr+=n13 )
			memcpy( m_new_ptr, m_ptr, length );
	}
    return m_new;
}
ComplexMatrix Cal_Periodic_Chi0::reshape_complexmat_21(const size_t n1, const size_t n2, const size_t n3, const ComplexMatrix &mat) //(n1,n2*n3) -> (n1*n2,n3)
{
    const auto length = sizeof(complex<double>)*n1*n2*n3;
	const auto n13 = n1*n2*n3;
	const complex<double> *m_ptr = mat.c;
	ComplexMatrix m_new( n1*n2 ,n3, false );
	complex<double> *m_new_ptr = m_new.c;
    memcpy( m_new_ptr, m_ptr, length );
    return m_new;
}


void Cal_Periodic_Chi0::generate_Green_part(const size_t &is, const double &time_tau, const Vector3_Order<int> &R, const size_t K_index, const size_t L_index, ComplexMatrix &Green_kl, ComplexMatrix &Green_kl_conj)
{
    //cout<<"    generate_Green_part   R:"<<R<<endl;
    const size_t k_num=atom_nw[K_index];
    const size_t l_num=atom_nw[L_index];
    //cout<<"time_tau:"<<time_tau<<"    R:"<<R<<endl;
    //cout<<"   generate_green_part:  k,l num:"<<k_num<<"   "<<l_num<<endl; 
    //cout<<" Green_kl.size:"<<Green_kl.size<<"    "<<Green_kl_conj.size<<endl;
    for(size_t k_iw=0;k_iw!=k_num;++k_iw)
    {
        size_t k_global=atom_iw_loc2glo(K_index,k_iw);
        //cout<<"loop_debug_1.7"<<endl;
        for(size_t l_iw=0;l_iw!=l_num;++l_iw)
        {
            //cout<<"loop_debug_2"<<endl;
            size_t l_global=atom_iw_loc2glo(L_index,l_iw);
            Green_kl(k_iw,l_iw)=Green_function[is][time_tau][R](k_global,l_global);
            Green_kl_conj(k_iw,l_iw)=conj(Green_function[is][-1*time_tau][R](k_global,l_global));
        
        } 
    }
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
                auto chi0_freq_R_I=I_pair.second;
                //const size_t mu_num=exx_lcao.index_abfs[ucell.iat2it[I]].count_size;
                for(auto &J_pair:chi0_freq_R_I)
                {
                    size_t J=J_pair.first;
                    //const auto chi0_mat=J_pair.second;
                    //(*Vps[I][J][R]).c
                    //cout<<"2"<<endl;
                   // const size_t nu_num=exx_lcao.index_abfs[ucell.iat2it[J]].count_size;
                    for(auto &q_pair:Vq[I][J])
                    {
                        //cout<<"3"<<endl;
                        Vector3_Order<double> ik_vec=q_pair.first;
                        auto &vq_mat=*q_pair.second;
                        pi_k[freq][ik_vec][I][J].create(vq_mat.nr,vq_mat.nc);
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
                        pi_k[freq][ik_vec][I][Q]+=transpose(Vq_mat,0)*chi0_k[freq][ik_vec][J][Q];
                        if(freq==first_freq)
                        {
                            cout<<I<<"  "<<Q<<"  "<<J<<"   "<<ik_vec<<endl;
                            print_complex_matrix("chi0_mat",chi0_k[freq][ik_vec][I][J]);
                            print_complex_matrix("vq",transpose(Vq_mat,0));
                            print_complex_matrix("pi_mat",pi_k[freq][ik_vec][I][J]);
                        }
                    }
                }
            }
        }
    }

    // temp_pi_count=0;                               
    // for(auto &freq_pair:chi0_freq)
    // {
    //     const double freq=freq_pair.first;
    //     const auto chi0_freq=freq_pair.second;
    //     for(auto &R_pair:chi0_freq)
    //     {
    //         const Vector3_Order<int> R=R_pair.first;
    //         auto chi0_freq_R=R_pair.second; 
    //         for(auto &I_pair:chi0_freq_R)
    //         {
    //             const size_t I=I_pair.first;
    //             auto chi0_freq_R_I=I_pair.second;
    //             const size_t mu_num=atom_mu[I];
    //             for(auto &J_pair:chi0_freq_R_I)
    //             {
    //                 size_t J=J_pair.first;
    //                 //const auto chi0_mat=J_pair.second;
    //                 //(*Vps[I][J][R]).c
                    
    //                 const size_t nu_num=atom_mu[J];
    //                 //cout<<tau<<"    "<<R<<"   "<<I<<"    "<<J<<"    mu_num="<<mu_num<<"    nu_num="<<nu_num<<endl;
    //                 //cout<<pi_freq[freq][R][I][J].size<<endl;
    //                 //pi_freq[freq][R][I][J].create(mu_num,nu_num);
    //                 //cout<<"  omega:"<<freq<<"    R"<<R<<endl;
    //                 for(auto &Q_pair:chi0_freq_R)
    //                 {
    //                     const auto Q=Q_pair.first;
    //                     const size_t xi_num=atom_mu[Q]; 
    //                     auto Vps_Q=Q_pair.second;
                        
    //                     auto chi0_mat=chi0_freq_R_I[Q];
                        

    //                     for(auto &Vq_pair:Vq[Q][J])
    //                     {
    //                         auto qvec=Vq_pair.first;
    //                         auto Vq_mat=*Vq_pair.second;
    //                         //cout<<"   Q="<<Q<<"    P="<<P<<"   R_Vps:"<<Vps_R<<"    xi_num="<<xi_num<<endl;
    //                        // ComplexMatrix Vs_complex(*Vps_mat_ptr);  //R in Vps_Q != R in chi0  
        
    //                         Vector3_Order<double> ik_vec(qvec); //latvec
    //                         const double arg = ( ik_vec* (R*latvec)) * TWO_PI;
    //                         const complex<double> kphase = complex <double> ( cos(arg),  sin(arg) );
                            
    //                         auto pi_tmp=chi0_mat*Vq_mat;
    //                         pi_tmp*=kphase;
    //                         LapackConnector::gemm(
    //                                 'N','N',
    //                                 mu_num,nu_num,xi_num,
    //                                 kphase,
    //                                 chi0_mat.c,xi_num,
    //                                 Vq_mat.c,nu_num,
    //                                 1,
    //                                 pi_k[freq][ik_vec][I][J].c,
    //                                 nu_num);
                            
    //                         if(freq==first_freq)
    //                         {
    //                             cout<<ik_vec<<kphase<<endl;
    //                             print_complex_real_matrix("chi0_mat",chi0_mat);
    //                             print_complex_real_matrix("vq",Vq_mat);
    //                             print_complex_real_matrix("pi_mat",pi_tmp);

    //                         }
                            
    //                     }
                      
    //                 }
    //             }
    //         }
    //     }
    //     temp_pi_count++;
    // }
    cout<<"Finish cal pi_k"<<endl;
    return pi_k;
}

/*
std::map<double,map<Abfs::Vector3_Order<double>,map<size_t,map<size_t,ComplexMatrix>>>> Cal_Periodic_Chi0::cal_pi_k()
{
    cout<<" bigin to cal_pi_k"<<endl;
    map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,shared_ptr<matrix>>>> Vs_m;		// Vs[iat1][iat2][box2]
	map<size_t,map<size_t,map<Abfs::Vector3_Order<int>,shared_ptr<matrix>>>> Vps_m;
    Vs_m = Abfs::cal_Vs( exx_lcao.atom_pairs_core_origin, exx_lcao.m_abfs_abfs, exx_lcao.index_abfs, exx_lcao.info.ccp_rmesh_times, exx_lcao.info.v_threshold, exx_lcao.Vws );
    Vps_m = Abfs::cal_mps( exx_lcao.Born_von_Karman_period, Vs_m );
    cout<<"1"<<endl;
    //std::map<double,map<Abfs::Vector3_Order<double>,map<size_t,map<size_t,ComplexMatrix>>>> pi_k;
    cout<<" n_kpoints="<<n_kpoints<<endl;
    for(auto &freq_pair:chi0_freq)
    {
        const double freq=freq_pair.first;
        const auto chi0_freq=freq_pair.second;
        for(auto &R_pair:chi0_freq)
        {
            const Abfs::Vector3_Order<int> R=R_pair.first;
            auto chi0_freq_R=R_pair.second; 
            for(auto &I_pair:chi0_freq_R)
            {
                const size_t I=I_pair.first;
                auto chi0_freq_R_I=I_pair.second;
                const size_t mu_num=exx_lcao.index_abfs[ucell.iat2it[I]].count_size;
                for(auto &J_pair:chi0_freq_R_I)
                {
                    size_t J=J_pair.first;
                    //const auto chi0_mat=J_pair.second;
                    //(*Vps[I][J][R]).c
                    //cout<<"2"<<endl;
                    const size_t nu_num=exx_lcao.index_abfs[ucell.iat2it[J]].count_size;
                    for(int ik=0;ik!=n_kpoints/NSPIN;++ik)
                    {
                        //cout<<"3"<<endl;
                        Abfs::Vector3_Order<double> ik_vec(kv.kvec_c[ik]);
                        pi_k[freq][ik_vec][I][J].create(mu_num,nu_num);
                        //cout<<freq<<"    "<<I<<"   "<<J<<"    ik:"<<ik<<"   ik_vec "<<ik_vec<<endl;
                    }
                }
            }
            break;
        }
    }
                     
    int temp_count=0;                
    for(auto &freq_pair:chi0_freq)
    {
        const double freq=freq_pair.first;
        const auto chi0_freq=freq_pair.second;
        for(auto &R_pair:chi0_freq)
        {
            const Abfs::Vector3_Order<int> R=R_pair.first;
            auto chi0_freq_R=R_pair.second; 
            for(auto &I_pair:chi0_freq_R)
            {
                const size_t I=I_pair.first;
                auto chi0_freq_R_I=I_pair.second;
                const size_t mu_num=exx_lcao.index_abfs[ucell.iat2it[I]].count_size;
                for(auto &J_pair:chi0_freq_R_I)
                {
                    size_t J=J_pair.first;
                    //const auto chi0_mat=J_pair.second;
                    //(*Vps[I][J][R]).c
                    
                    const size_t nu_num=exx_lcao.index_abfs[ucell.iat2it[J]].count_size;
                    //cout<<tau<<"    "<<R<<"   "<<I<<"    "<<J<<"    mu_num="<<mu_num<<"    nu_num="<<nu_num<<endl;
                    //cout<<pi_freq[freq][R][I][J].size<<endl;
                    //pi_freq[freq][R][I][J].create(mu_num,nu_num);
                    //cout<<"  omega:"<<freq<<"    R"<<R<<endl;
                    for(auto &Q_pair:chi0_freq_R)
                    {
                        const auto Q=Q_pair.first;
                        const size_t xi_num=exx_lcao.index_abfs[ucell.iat2it[Q]].count_size; 
                        auto Vps_Q=Q_pair.second;
                        
                        auto chi0_mat=chi0_freq_R_I[Q];
                        
                        if(Q<=J)
                        {
                            for(auto &Vps_R_pair:Vps_m[Q][J])
                            {
                                auto Vps_R=Vps_R_pair.first;
                                auto Vps_mat_ptr=Vps_R_pair.second;
                                //cout<<"   Q="<<Q<<"    P="<<P<<"   R_Vps:"<<Vps_R<<"    xi_num="<<xi_num<<endl;
                                ComplexMatrix Vs_complex(*Vps_mat_ptr);  //R in Vps_Q != R in chi0  
                               
                                
                                Abfs::Vector3_Order<int> R_plus_R2(R+Vps_R);
                                Abfs::Vector3_Order<double> R_plus_dl(R_plus_R2.x,R_plus_R2.y,R_plus_R2.z);
                                for(size_t ik=0;ik!=n_kpoints/NSPIN;++ik)
                                {
                                    Abfs::Vector3_Order<double> ik_vec(kv.kvec_c[ik]);
                                    const double arg = -1*( ik_vec* R_plus_dl) * TWO_PI;
                                    const complex<double> kphase = complex <double> ( cos(arg),  sin(arg) );
                                    
                                    
                                    LapackConnector::gemm(
                                            'N','N',
                                            mu_num,nu_num,xi_num,
                                            kphase,
                                            chi0_mat.c,xi_num,
                                            Vs_complex.c,nu_num,
                                            1,
                                            pi_k[freq][ik_vec][I][J].c,
                                            nu_num);
                                    // if( ik_vec==Abfs::Vector3_Order<double>(0,0,0) && I==0 && J==0)
                                    // {
                                        // cout<<freq<<"  Q:"<<Q<<" R,R2:"<<R<<Vps_R<<"  R_plus_dl"<<R_plus_dl<<"  kphase"<<kphase<<" pi_k:"<<pi_k[freq][ik_vec][I][J](0,0)<<endl;
                                    // }   
                                    // if(temp_count == 0) 
                                    // {
                                        // cout<<"  freq="<<freq<<"   I J Q:"<<I<<J<<Q<<"     R:"<<R<<"   Vps_R:"<<Vps_R<<"  ik_vec:"<<ik_vec<<"   arg:"<<arg<<"  kphase:"<<kphase<<endl;
                                        // print_complex_matrix("chi0_mat",chi0_mat);
                                        // print_complex_matrix("Vps_mat",Vs_complex);
                                    // }
                                    
                                }
                            }
                        }
                        if(Q>J)
                        {
                            for(auto &Vps_R_pair:Vps_m[J][Q])
                            {
                                auto Vps_R=Vps_R_pair.first;
                                auto Vps_mat_ptr=Vps_R_pair.second;
                                //cout<<"   Q="<<Q<<"    P="<<P<<"   R_Vps:"<<Vps_R<<"    xi_num="<<xi_num<<endl;
                                ComplexMatrix Vs_complex(transpose(*Vps_mat_ptr));  //R in Vps_Q != R in chi0  
                                //print_complex_matrix("chi0_mat",chi0_mat.nc,chi0_mat.nr,chi0_mat.c,chi0_mat.nc);
                                //print_complex_matrix("Vps_mat",Vs_complex.nc,Vs_complex.nr,Vs_complex.c,Vs_complex.nc);
                                Abfs::Vector3_Order<int> R_plus_R2(R+Vps_R);
                                Abfs::Vector3_Order<double> R_plus_dl(R_plus_R2.x,R_plus_R2.y,R_plus_R2.z);
                                for(size_t ik=0;ik!=n_kpoints/NSPIN;++ik)
                                {
                                    Abfs::Vector3_Order<double> ik_vec(kv.kvec_c[ik]);
                                    const double arg =-1* ( ik_vec* R_plus_dl) * TWO_PI;
                                    const complex<double> kphase = complex <double> ( cos(arg),  sin(arg) );
                                    
                                    //cout<<freq<<"    "<<I<<"   "<<J<<"  pi_k  kphase"<<kphase<<endl;
                                    LapackConnector::gemm(
                                            'N','N',
                                            mu_num,nu_num,xi_num,
                                            kphase,
                                            chi0_mat.c,xi_num,
                                            Vs_complex.c,nu_num,
                                            1,
                                            pi_k[freq][ik_vec][I][J].c,
                                            nu_num);
                                            
                                    if( ik_vec==Abfs::Vector3_Order<double>(0,0,0) && I==0 && J==0)
                                    {
                                        cout<<freq<<"  Q:"<<Q<<" R,R2:"<<R<<Vps_R<<"  R_plus_dl"<<R_plus_dl<<"  kphase"<<kphase<<" pi_k:"<<pi_k[freq][ik_vec][I][J](0,0)<<endl;
                                    }
                                    // if(temp_count == 0) 
                                    // {
                                        // cout<<"  **freq="<<freq<<"   I J Q:"<<I<<J<<Q<<"     R:"<<R<<"   Vps_R:"<<Vps_R<<"  ik_vec:"<<ik_vec<<"   arg:"<<arg<<"  kphase:"<<kphase<<endl;
                                        // print_complex_matrix("chi0_mat",chi0_mat);
                                        // print_complex_matrix("Vps_mat",Vs_complex);
                                        // print_complex_matrix("pi_freq_mat",pi_k[freq][ik_vec][I][J]);
                                    // }
                                
                                }
                            }
                        }
                        
                    }
                    
                }
            }
        }
        temp_count++;
    }
    return pi_k;
}
*/
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
                    for(auto &q_pair:Vq[I][J])
                    {
                        //cout<<"3"<<endl;
                        Vector3_Order<double> ik_vec=q_pair.first;
                        auto &vq_mat=*q_pair.second;
                        chi0_k[freq][ik_vec][I][J].create(vq_mat.nr,vq_mat.nc);
                        //cout<<freq<<"    "<<I<<"   "<<J<<"    ik:"<<ik<<"   ik_vec "<<ik_vec<<endl;
                    }
                    
                    // for(int ik=0;ik!=n_kpoints;ik++)
                    // {
                    //     //cout<<"3"<<endl;
                    //     Vector3_Order<double> ik_vec=kvec_c[ik];
                    //     chi0_k[freq][ik_vec][I][J].create(atom_mu[I],atom_mu[J]);
                    //     //cout<<freq<<"    "<<I<<"   "<<J<<"    ik:"<<ik<<"   ik_vec "<<ik_vec<<endl;
                    // }
                }
            }
            break;
        }
    }
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
                    for(auto &ikp:chi0_k[freq])
                    {
                        auto ik_vec=ikp.first;
                        const double arg = ( ik_vec* (R*latvec)) * TWO_PI;
                        const complex<double> kphase = complex <double> ( cos(arg),  sin(arg) );
                        chi0_k[freq][ikp.first][I][J]+=kphase*chi0_mat;

                        //if(temp_chi0_k_count<=0)
                        //{
                            //cout<<ik_vec<<"  "<<R<<"  "<<kphase<<endl;
                            //print_complex_matrix("temp_chi0_k ",(kphase*chi0_mat));
                        //}
                        
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
/*

std::map<double,map<Abfs::Vector3_Order<double>,map<size_t,map<size_t,ComplexMatrix>>>> Cal_Periodic_Chi0::cal_pi_k_tau()
{
    cout<<" bigin to cal_pi_k"<<endl;
    cout<<"1"<<endl;
    //std::map<double,map<Abfs::Vector3_Order<double>,map<size_t,map<size_t,ComplexMatrix>>>> pi_k;
    cout<<" n_kpoints="<<n_kpoints<<endl;
    for(auto &tau_pair:chi0)
    {
        const double tau=tau_pair.first;
        const auto chi0_tau=tau_pair.second;
        for(auto &R_pair:chi0_tau)
        {
            const Abfs::Vector3_Order<int> R=R_pair.first;
            auto chi0_tau_R=R_pair.second; 
            for(auto &I_pair:chi0_tau_R)
            {
                const size_t I=I_pair.first;
                auto chi0_tau_R_I=I_pair.second;
                const size_t mu_num=atom_mu[I];
                for(auto &J_pair:chi0_tau_R_I)
                {
                    size_t J=J_pair.first;
                    //const auto chi0_mat=J_pair.second;
                    //(*Vps[I][J][R]).c
                    //cout<<"2"<<endl;
                    const size_t nu_num=atom_mu[J];
                    for(int ik=0;ik!=n_kpoints/NSPIN;++ik)
                    {
                        //cout<<"3"<<endl;
                        Vector3_Order<double> ik_vec(kv.kvec_c[ik]);
                        pi_k[tau][ik_vec][I][J].create(mu_num,nu_num);
                        //cout<<freq<<"    "<<I<<"   "<<J<<"    ik:"<<ik<<"   ik_vec "<<ik_vec<<endl;
                    }
                }
            }
            break;
        }
    }
                     
    int temp_count=0;                
    for(auto &tau_pair:chi0)
    {
        const double tau=tau_pair.first;
        const auto chi0_tau=tau_pair.second;
        for(auto &R_pair:chi0_tau)
        {
            const Vector3_Order<int> R=R_pair.first;
            auto chi0_tau_R=R_pair.second; 
            for(auto &I_pair:chi0_tau_R)
            {
                const size_t I=I_pair.first;
                auto chi0_tau_R_I=I_pair.second;
                const size_t mu_num=atom_m[I];
                for(auto &J_pair:chi0_tau_R_I)
                {
                    size_t J=J_pair.first;
                    //const auto chi0_mat=J_pair.second;
                    //(*Vps[I][J][R]).c
                    
                    const size_t nu_num=atom_m[J];
                   
                    for(auto &Q_pair:chi0_tau_R)
                    {
                        const auto Q=Q_pair.first;
                        const size_t xi_num=atom_mu[Q]; 
                        auto Vps_Q=Q_pair.second;
                        
                        auto chi0_mat=chi0_tau_R_I[Q];
                        
                        if(Q<=J)
                        {
                            for(auto &Vps_R_pair:Vq[Q][J])
                            {
                                auto Vps_R=Vps_R_pair.first;
                                auto Vps_mat_ptr=Vps_R_pair.second;
                                //cout<<"   I="<<I<<"    J="<<J<<"     Q:"<<Q<<"     R:"<<R<<"   R_Vps:"<<Vps_R<<endl;
                                ComplexMatrix Vs_complex(*Vps_mat_ptr);  //R in Vps_Q != R in chi0  
                               
                                
                                Abfs::Vector3_Order<int> R_plus_R2(R+Vps_R);
                                Abfs::Vector3_Order<double> R_plus_dl(R_plus_R2.x,R_plus_R2.y,R_plus_R2.z);
                                for(size_t ik=0;ik!=n_kpoints/NSPIN;++ik)
                                {
                                    Abfs::Vector3_Order<double> ik_vec(kv.kvec_c[ik]);
                                    const double arg = ( ik_vec* R_plus_dl) * TWO_PI;
                                    const complex<double> kphase = complex <double> ( cos(arg),  sin(arg) );
                                    
                                    //cout<<freq<<"    "<<I<<"   "<<J<<"  pi_k  kphase"<<kphase<<endl;
                                    LapackConnector::gemm(
                                            'N','N',
                                            mu_num,nu_num,xi_num,
                                            kphase,
                                            chi0_mat.c,xi_num,
                                            Vs_complex.c,nu_num,
                                            1,
                                            pi_k[tau][ik_vec][I][J].c,
                                            nu_num);
                                       
                                    // if(temp_count == 0) 
                                    // {
                                        // cout<<"  tau="<<tau<<"   I J Q:"<<I<<J<<Q<<"     R:"<<R<<"   Vps_R:"<<Vps_R<<"  ik_vec:"<<ik_vec<<"   arg:"<<arg<<"  kphase:"<<kphase<<endl;
                                        // print_complex_matrix("chi0_mat",chi0_mat);
                                        // print_complex_matrix("Vps_mat",Vs_complex);
                                    // }
                                    
                                }
                            }
                        }
                        if(Q>J)
                        {
                            for(auto &Vps_R_pair:Vps_m[J][Q])
                            {
                                auto Vps_R=Vps_R_pair.first;
                                auto Vps_mat_ptr=Vps_R_pair.second;
                                //cout<<"   I="<<I<<"    J="<<J<<"     Q:"<<Q<<"     R:"<<R<<"   R_Vps:"<<Vps_R<<endl;
                                ComplexMatrix Vs_complex(transpose(*Vps_mat_ptr));  //R in Vps_Q != R in chi0  
                                //print_complex_matrix("chi0_mat",chi0_mat.nc,chi0_mat.nr,chi0_mat.c,chi0_mat.nc);
                                //print_complex_matrix("Vps_mat",Vs_complex.nc,Vs_complex.nr,Vs_complex.c,Vs_complex.nc);
                                Abfs::Vector3_Order<int> R_plus_R2(R+Vps_R);
                                Abfs::Vector3_Order<double> R_plus_dl(R_plus_R2.x,R_plus_R2.y,R_plus_R2.z);
                                for(size_t ik=0;ik!=n_kpoints/NSPIN;++ik)
                                {
                                    Abfs::Vector3_Order<double> ik_vec(kv.kvec_c[ik]);
                                    const double arg = ( ik_vec* R_plus_dl) * TWO_PI;
                                    const complex<double> kphase = complex <double> ( cos(arg),  sin(arg) );
                                    
                                    //cout<<freq<<"    "<<I<<"   "<<J<<"  pi_k  kphase"<<kphase<<endl;
                                    LapackConnector::gemm(
                                            'N','N',
                                            mu_num,nu_num,xi_num,
                                            kphase,
                                            chi0_mat.c,xi_num,
                                            Vs_complex.c,nu_num,
                                            1,
                                            pi_k[tau][ik_vec][I][J].c,
                                            nu_num);
                                            
                                    // if(temp_count == 0) 
                                    // {
                                        // cout<<"  **tau="<<tau<<"   I J Q:"<<I<<J<<Q<<"     R:"<<R<<"   Vps_R:"<<Vps_R<<"  ik_vec:"<<ik_vec<<"   arg:"<<arg<<"  kphase:"<<kphase<<endl;
                                        // print_complex_matrix("chi0_mat",chi0_mat);
                                        // print_complex_matrix("Vps_mat",Vs_complex);
                                        // print_complex_matrix("pi_freq_mat",pi_k[tau][ik_vec][I][J]);
                                    // }
                                
                                }
                            }
                        }
                        
                    }
                    
                }
            }
        }
        temp_count++;
    }
    return pi_k;
}

void Cal_Periodic_Chi0::cal_MP2_energy_pi_k_tau( map<double,double> &tau_grid ,std::map<double,map<Abfs::Vector3_Order<double>,map<size_t,map<size_t,ComplexMatrix>>>> pi_k ) 
{
    //std::map<double,map<Abfs::Vector3_Order<double>,map<size_t,map<size_t,ComplexMatrix>>>> pi_k(cal_pi_k());
    std::map<double,map<Abfs::Vector3_Order<double>,map<size_t,map<size_t,ComplexMatrix>>>> pi_k_square;
    for(auto &tau_pair:pi_k)
    {
        const double tau=tau_pair.first;
        const auto pi_tau=tau_pair.second;
        for(auto &k_pair:pi_tau)
        {
            const auto k=k_pair.first;
            auto pi_tau_k=k_pair.second; 
            for(auto &I_pair:pi_tau_k)
            {
                const size_t I=I_pair.first;
                auto pi_tau_k_I=I_pair.second;
                const size_t mu_num=exx_lcao.index_abfs[ucell.iat2it[I]].count_size;
                pi_k_square[tau][k][I][I].create(mu_num,mu_num);
                for(auto &J_pair:pi_tau_k_I)
                {
                    const size_t J=J_pair.first;
                    const auto pi_mat=J_pair.second;
                    const auto pi2_mat=pi_tau_k[J][I];
                    const size_t nu_num=exx_lcao.index_abfs[ucell.iat2it[J]].count_size;
                   // cout<<freq<<"    "<<R<<"   "<<I<<"    "<<J<<"    mu_num="<<mu_num<<"    nu_num="<<nu_num<<endl;
                    
                    
                    const size_t xi_num=exx_lcao.index_abfs[ucell.iat2it[J]].count_size;
                    
                    LapackConnector::gemm(
                        'N','N',
                        mu_num,mu_num,xi_num,
                        1,
                        pi_mat.c,xi_num,
                        pi2_mat.c,mu_num,
                        1,
                        pi_k_square[tau][k][I][I].c, 
                        mu_num);
             
                }
                //print_complex_matrix("pi_square_mat",pi_square[freq][R][I][I].nc,pi_square[freq][R][I][I].nr,pi_square[freq][R][I][I].c,pi_square[freq][R][I][I].nc);
            }
        }
    }
    
    complex<double> tot_trace_pi(0.0,0.0);
    for(auto &tau_pair:pi_k_square) 
    {
        auto tau=tau_pair.first;
        auto pi_tau=tau_pair.second;
        complex<double> trace_pi(0.0,0.0);
        for(auto &k_pair:pi_tau)
        {
            auto pi_tau_k=k_pair.second;
            for(auto &I_pair:pi_tau_k)
            {
                const size_t I=I_pair.first;
                auto pi_tau_k_I=I_pair.second;
                trace_pi+=trace(pi_tau_k_I[I]);
            }
        }
        tot_trace_pi+=trace_pi*tau_grid[tau]*NSPIN/n_kpoints/2*(-1);
        cout<<" pi_k   tau:"<<tau<<"   trace:"<<trace_pi<<endl;
    }
    cout<<" tot_tau_MP2_energy_pi_k: "<<tot_trace_pi;
}
*/
void Cal_Periodic_Chi0::RPA_correlation_energy( map<double,double> &freq_grid )
{
    int range_all=0;
/* 	for(int Mu=0;Mu!=ucell.nat;Mu++)
	{
		//cout<<"tau["<<Mu<<"]="<<ucell.atoms[ucell.iat2it[Mu]].tau[ucell.iat2ia[Mu]]<<endl;
		range_all+=exx_lcao.index_abfs[ucell.iat2it[Mu]].count_size;
	    //cout<<"index_abfs["<<Mu<<"].count_size="<<exx_lcao.index_abfs[ucell.iat2it[Mu]].count_size<<endl;
	} */
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
    for(const auto &freq_pair:pi_k)
    {
        const auto freq=freq_pair.first;
        const auto pi_freq=freq_pair.second;
        
        for(const auto &k_pair:pi_freq)
        {
            const auto kvec_c=k_pair.first;
            const auto pi_freq_k=k_pair.second;
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
                            pi_freq_2D[freq][kvec_c](part_range[I]+mu,part_range[J]+nu)+=mat(mu,nu); 
                        }
                    }
                }
            }
        }
    } 
    
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
            trace_pi=trace(pi_freq_2D[freq][kvec_c]);
            rpa_for_omega_k=ln_det+trace_pi; 
            //cout<<" freq:"<<freq<<"      rpa_for_omega_k: "<<rpa_for_omega<<endl;
            cRPA_k[kvec_c]+=rpa_for_omega_k*weight*irk_wight[kvec_c]*NSPIN/TWO_PI;
            tot_RPA_energy+=rpa_for_omega_k*weight*irk_wight[kvec_c]*NSPIN/TWO_PI;
        }
    }
    for(auto &ikvec:cRPA_k)
    {
        cout<<ikvec.first<<ikvec.second<<endl;
    }
    // for(const auto &freq_pair:freq_grid)
    // {
        // const auto freq=freq_pair.first;
        // const auto weight=freq_pair.second;
        
        // complex<double> rpa_for_omega(0.0,0.0);
        // ComplexMatrix identity(range_all,range_all);
        // ComplexMatrix identity_minus_pi(range_all,range_all);
       
        // identity.set_as_identity_matrix();
       
        // identity_minus_pi=identity-pi_freq_2D[freq];
        // complex<double> det_for_rpa(1.0,0.0);
        // int info_LU =0;
        // int *ipiv=new int[range_all];
        // LapackConnector::zgetrf(range_all,range_all,identity_minus_pi,range_all,ipiv,&info_LU);
        // for(int ib=0;ib!=range_all;ib++)
        // {
            // if(ipiv[ib] !=(ib+1)) det_for_rpa=-det_for_rpa*identity_minus_pi(ib,ib);
            // else det_for_rpa=det_for_rpa*identity_minus_pi(ib,ib);
        // }
        // delete[] ipiv;
        
        // complex<double> trace_pi;
        // complex<double> ln_det;
        // ln_det=std::log(det_for_rpa);
        // trace_pi=trace(pi_freq_2D[freq]);
        // rpa_for_omega=ln_det+trace_pi; 
        // cout<<" freq:"<<freq<<"      rpa_for_omega: "<<rpa_for_omega<<endl;
        // tot_RPA_energy+=rpa_for_omega*weight/TWO_PI;
    // }
    cout<<"Grid_num_"<<freq_grid.size()<<"  tot_RPA_energy:  "<<setprecision(8)<<tot_RPA_energy<<endl; 
}
/*
Abfs::Vector3_Order<int> Cal_Periodic_Chi0::periodic_R(const Abfs::Vector3_Order<int> &R)
{
    Abfs::Vector3_Order<int> Rp;
    if(R.x>kv_nmp[0])
    {
        Rp.x=(R.x+kv_nmp[0])%(2*kv_nmp[0]+1)-kv_nmp[0];
    }
    else if(R.x<-kv_nmp[0])
    {
        Rp.x=(R.x-kv_nmp[0])%(2*kv_nmp[0]+1)+kv_nmp[0];
    }
    else
    {
        Rp.x=R.x;
    }
    
    if(R.y>kv_nmp[1])
    {
        Rp.y=(R.y+kv_nmp[1])%(2*kv_nmp[1]+1)-kv_nmp[1];
    }
    else if(R.y<-kv_nmp[1])
    {
        Rp.y=(R.y-kv_nmp[1])%(2*kv_nmp[1]+1)+kv_nmp[1];
    }
    else
    {
        Rp.y=R.y;
    }
    
    if(R.z>kv_nmp[2])
    {
        Rp.z=(R.z+kv_nmp[2])%(2*kv_nmp[2]+1)-kv_nmp[2];
    }
    else if(R.z<-kv_nmp[2])
    {
        Rp.z=(R.z-kv_nmp[2])%(2*kv_nmp[2]+1)+kv_nmp[2];
    }
    else
    {
        Rp.z=R.z;
    }
    
    return Rp;
}

map<double, double> Cal_Periodic_Chi0::get_minimax_grid(const size_t &npoints)
{
    vector<double> omega,womega;
    switch(npoints)
    {
        case 10:{
                omega={ 0.328674516168464881,
   1.15613158254563664 ,
   2.62856154125761776 ,
   5.74937217095778053 ,
   13.0248895967233107 ,
   31.4326650917899961 ,
   82.7766397660356006 ,
  246.534887108615123  ,
  889.397043040403901  ,
 4653.65622419843385   };
                womega={  0.217621574654979028,
  0.332538790770825932,
  0.653721363306432379,
  1.45317967871747711 ,
  3.50610171553335448 ,
  9.20598376115589545 ,
 26.9325507894137992  ,
 91.9027380611465077  ,
 402.033401296635873  ,
2949.25231903565054   }; break;
                }
        case 24:{
                omega={ 0.014965706746460873 , 
0.2947974622113442 , 
1.9739210092351913 , 
4.291146210471523 ,  
6.410668763699374 ,  
11.250490198626732 , 
18.744306431523714 , 
31.807131175591113 , 
53.30110495764749 ,  
88.6932230770902 ,   
147.12514987141506 , 
243.38417288969 ,    
402.13027063158114 , 
663.8814331260972 ,  
1095.3887464093016 , 
1806.8092044651203 , 
2979.729287006711 ,  
4913.559195806112 ,  
8101.8869253074645 , 
13358.535120824672 , 
22025.294475232186 , 
36314.344601809076 , 
59872.99694075114 ,  
98714.6400338346 };
                womega={      0.07587421985052269 ,
    0.6323913474664682 ,
    1.4266722495302258 ,
   190.63599836650423 ,
   24.91460585765506 ,
    8.977752919953629 ,
    19.200447274024235 ,
    32.78153950437948 ,
   53.198394619610305 ,
  88.2140478707364 ,
    147.12913350531045 ,
 243.57559549283005 ,
    402.43271093309784 ,
   663.914036821879 ,
    1095.7077076850487 ,
    1807.0111241097388 ,
   2980.430295674161 ,
   4914.052127586602 ,
    8102.114308232925 ,
    13358.765379845956 ,
    22025.259783235728 ,
    36314.706088716426 ,
   59873.433510239025 ,
  98714.49640280618 }; break;
        }
    }
    double gap;
    gap=get_band_gap();
    map<double,double> minimax_grid;
    minimax_grid.clear();
    for(size_t i=0;i!=omega.size();++i)
        minimax_grid.insert({gap*omega[i],PI*gap*womega[i]});
    cout<<" MINIMAX_GRID_Freq "<<endl;
    for(const auto &m:minimax_grid)
        cout<<m.first<<"      "<<m.second<<endl;
    return minimax_grid;
}

map<double, double> Cal_Periodic_Chi0::get_minimax_grid_tau(const size_t &npoints)
{
    vector<double> omega,womega;
    switch(npoints)
    {
        case 10:{
                omega={0.0056851434686149525 ,
0.03150111614201419 ,  
0.08533148884551203 ,  
0.18409078313643584 ,  
0.3605269668520229 ,   
0.6679993384336078 ,   
1.1067787325810754 ,   
1.37324861153843 ,     
2.265933442758568 ,    
4.064769311966836   };
                womega={         0.014756769924543464 ,
    0.03800071747714435 ,
    0.07237208887758975 ,
    0.13050050087459458 ,
   0.2313730498283273 ,
   0.3938993555249255 ,
   0.34560988616609883 ,
 0.4700665327672719 ,
  1.2704476680517849 ,
  2.52649597296431  }; break;
                }
        case 24:{
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
    }
    double gap;
    gap=get_band_gap();
    map<double,double> minimax_grid;
    minimax_grid.clear();
    for(size_t i=0;i!=omega.size();++i)
        minimax_grid.insert({omega[i]/(gap),womega[i]/(gap)});
    cout<<" MINIMAX_GRID_Tau "<<endl;
    for(const auto &m:minimax_grid)
        cout<<m.first<<"      "<<m.second<<endl;
    return minimax_grid;
}
*/

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
        for( int j = 0; j < nc; j++ ) printf("%10.6f,%6.4f ",mat.c[i*nc+j].real(),mat.c[i*nc+j].imag());
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