#include "Hamiltonian.h"
#include "fermi_energy_occupation.h"
#include "fermi_energy_occupation.h"
#include <iostream>
#include "pbc.h"
#include "constants.h"


// 定义复数类型
using cplxdb = std::complex<double>;

//G0, for scRPA
std::vector<std::vector<cplxdb>> build_G0(
    MeanField& meanfield,
    const std::vector<double>& freq_nodes,
    int ispin,
    int ikpt,
    int n_states) {
    // 获取虚频点列表和能量列表
    std::vector<double> eigenvals(n_states);

    for (int n = 0; n < n_states; ++n) {
        eigenvals[n] = meanfield.get_eigenvals()[ispin](ikpt, n);
    }

    // 创建 G0 矩阵，大小为 n_states x freq_nodes.size()
    std::vector<std::vector<cplxdb>> G0(n_states, std::vector<cplxdb>(freq_nodes.size()));

    for (int n = 0; n < n_states; ++n) {
        for (size_t w = 0; w < freq_nodes.size(); ++w) {
            // 计算 G_n^0(iω) = 1 / (iω - ε_n)
            cplxdb iw(0.0, freq_nodes[w]);  // 虚频 iω
            G0[n][w] = 1.0 / (iw - eigenvals[n]);
        }
    }
    return G0;
}
// mode A
Matz build_correlation_potential_spin_k_modeA(
    const std::vector<std::vector<std::vector<cplxdb>>>& sigc_spin_k,
    int n_states) {
    std::cout << "QSGW: mode A" << std::endl;
    Matz Vc_spin_k(n_states, n_states, MAJOR::COL);
    std::map<int, Matz> Re_sigma;
    std::map<int, Matz> sigma;
    std::cout << "check112 " << std::endl;
    for (int k = 0 ; k < n_states; ++k)
    {
        Re_sigma[k] = Matz(n_states, n_states, MAJOR::COL);
        sigma[k] = Matz(n_states, n_states, MAJOR::COL);
        for (int i = 0; i < n_states; ++i)
        {
            for (int j = 0; j < n_states; ++j)
            {              
                sigma[k](i,j) = sigc_spin_k[i][j][k] ;    
            }
        }
        Re_sigma[k] = 0.5 * (sigma[k] + transpose(sigma[k],true));
    }
    std::cout << "check113 " << std::endl;
    for (int i = 0; i < n_states; ++i)
    {
        for (int j = 0; j < n_states; ++j)
        { 
            // 构建关联势矩阵
            std::complex<double> Vc_ij = 0.5 * (Re_sigma[i](i,j) + Re_sigma[j](i,j));
            Vc_spin_k(i, j) = Vc_ij;  
        }
    }    

    return Vc_spin_k;
}

//mode B
Matz build_correlation_potential_spin_k(
    const std::vector<std::vector<std::vector<cplxdb>>>& sigc_spin_k,
    int n_states) {
    // std::cout << "QSGW: mode B" << std::endl;
    Matz Vc_spin_k(n_states, n_states, MAJOR::COL);
    std::map<int, Matz> Re_sigma;
    std::map<int, Matz> sigma;
    for (int k = 0 ; k < n_states+1; ++k)
    {
        Re_sigma[k] = Matz(n_states, n_states, MAJOR::COL);
        sigma[k] = Matz(n_states, n_states, MAJOR::COL);
        for (int i = 0; i < n_states; ++i)
        {
            for (int j = 0; j < n_states; ++j)
            {              
                sigma[k](i,j) = sigc_spin_k[i][j][k] ;    
            }
        }
        Re_sigma[k] = 0.5 * (sigma[k] + transpose(sigma[k],true));
    }

    for (int i = 0; i < n_states; ++i)
    {
        
        std::complex<double> Vc_ii = Re_sigma[i](i,i);
        Vc_spin_k(i, i) = Vc_ii;
        for (int j = 0; j < n_states; ++j)
        {
            if(i!=j){
                
                std::complex<double> Vc_ij = Re_sigma[n_states](i,j);
                Vc_spin_k(i, j) = Vc_ij;  
            }  
        }
    }
    return Vc_spin_k;
}


// 计算scRPA交换-关联能量矩阵
Matz calculate_scRPA_exchange_correlation(
    MeanField& meanfield,
    const std::vector<double>& freq_nodes,
    const std::vector<double>& freq_weights,
    const std::map<double, Matz>& sigc_spin_k, 
    const std::vector<std::vector<cplxdb>>& G0, 
    int ispin,
    int ikpt,
    int n_states, 
    double mu, 
    double temperature) 
{
    Matz V_rpa_ks(n_states, n_states, MAJOR::COL); // 交换-关联能量矩阵，列优先存储
    const double K_B = 3.16681e-6;  // Hartree/K
    const double threshold = 1e-5;  // 设定一个阈值，如果差值太小就不进行除法计算
    const double sigma = 0.0001;  // 用于计算 δfn / δεn 的高斯展宽
    for (int n = 0; n < n_states; ++n) {
        // 获取填充因子 fn
        double energy_n = meanfield.get_eigenvals()[ispin](ikpt, n);
        double f_n = fermi_dirac(energy_n, mu, temperature) * 2.0 / meanfield.get_n_spins();
        double efermi = meanfield.get_efermi();
        for (int m = 0; m < n_states; ++m) {
            cplxdb V_nm = 0.0;
            
            if (n == m) {
                // // 计算 δfn / δεn

                // double df_n_depsilon_n = 1/ (sqrt(2.0 * PI ) * sigma) * exp(-pow((energy_n - efermi),2)/(2 * sigma));

                // // 如果 df_n_depsilon_n 太小，则跳过该计算
                // if (std::abs(df_n_depsilon_n) < threshold) {
                //     V_rpa_ks(n, n) = 0.0;
                //     continue;  // 跳过这个循环的进一步计算
                // }

                // if(energy_n > efermi){
                //     V_rpa_ks(n, n) = 0.0;
                //     continue;
                // }
                
                // 累加自能项
                for (size_t w = 0; w < freq_weights.size(); ++w) {
                    cplxdb sigc_nn_iw = sigc_spin_k.at(freq_nodes[w])(n, n);
                    
                    // n = m 的对角元计算
                    // V_nm += freq_weights[w] * sigc_nn_iw * G0[n][w] * G0[n][w];
                    V_nm += freq_weights[w] * sigc_nn_iw * G0[n][w];
                }

                // 对 δfn / δεn 进行归一化处理
                // V_rpa_ks(n, n) = V_nm / df_n_depsilon_n;
                V_rpa_ks(n, n) = V_nm;
            } 
            else {
                // 获取填充因子 fm
                double energy_m = meanfield.get_eigenvals()[ispin](ikpt, m);
                double f_m = fermi_dirac(energy_m, mu, temperature) * 2.0 / meanfield.get_n_spins();
                // double delta_nm = (f_n - f_m)/(energy_n - energy_m);
                double delta_nm = f_n - f_m;
                // 如果 f_n - f_m 太小，则跳过该计算
                if (std::abs(delta_nm) < threshold) {
                    V_rpa_ks(n, m) = 0.0;
                    continue;  // 跳过这个循环的进一步计算
                }

                // 累加交换-关联项
                for (size_t w = 0; w < freq_weights.size(); ++w) {
                    cplxdb sigc_nm_iw = sigc_spin_k.at(freq_nodes[w])(n, m);
                    V_nm += freq_weights[w] * sigc_nm_iw * G0[n][w] * G0[m][w];
                }

                // 归一化并存储
                V_rpa_ks(n, m) = V_nm * (energy_n - energy_m) / delta_nm;
            }  
        }
    }
    V_rpa_ks = 0.5 * (V_rpa_ks + transpose(V_rpa_ks,true));

    return V_rpa_ks;
}



std::map<int, std::map<int, Matz>> construct_H0_GW(
    MeanField& meanfield,
    const std::map<int, std::map<int, Matz>> & H_KS_all,
    const std::map<int, std::map<int, Matz>> & vxc_all,
    const std::map<int, std::map<int, Matz>> & Hexx_all,
    const std::map<int, std::map<int, Matz>> & Vc_all,
    int n_spins, int n_kpoints, int n_states) {

    // 初始化 GW 哈密顿量矩阵
    std::map<int, std::map<int, Matz>> H0_GW_all;
    
    double efermi = meanfield.get_efermi();
    for (int ispin = 0; ispin < n_spins; ++ispin)
    {
        for (int ikpt = 0; ikpt < n_kpoints; ++ikpt)
        {
            Matz Hexx_ispin_ik = Hexx_all.at(ispin).at(ikpt);
            Matz Vxc_construct_ispin_ik = Hexx_ispin_ik + Vc_all.at(ispin).at(ikpt);
            // // realize in k-space
            // for (int i = 0; i < n_states; ++i){   
            //     for (int j = 0; j < n_states; ++j){
            //         Vxc_construct_ispin_ik(i,j) = std::real(Vxc_construct_ispin_ik(i,j));
            //         // Vxc_construct_ispin_ik(i,j) = std::real(Vc_all.at(ispin).at(ikpt)(i,j)) + Hexx_ispin_ik(i,j);
            //     }
            // }
            
            // cut if possible
            // Matz Vxc_diff_spin_k = Hexx_ispin_ik + Vc_all.at(ispin).at(ikpt) - vxc_all.at(ispin).at(ikpt);
            // for (int i = 0; i < n_states; ++i){
            //     double energy_i = meanfield.get_eigenvals()[ispin](ikpt, i);
            //     for (int j = 0; j < n_states; ++j){
            //         double energy_j = meanfield.get_eigenvals()[ispin](ikpt, j);
            //         if(energy_i > efermi+1.75||energy_j > efermi+1.75){
            //             Vxc_diff_spin_k(i, j)= 0.0 ;                    
            //         }
            //     }
            // }
            // 构建 GW 哈密顿量矩阵
            Matz H0_GW_spin_k = H_KS_all.at(ispin).at(ikpt) - vxc_all.at(ispin).at(ikpt) + Vxc_construct_ispin_ik;

            H0_GW_all[ispin][ikpt] = H0_GW_spin_k;   
        }
    }

    return H0_GW_all;
}

std::map<int, std::map<int, Matz>> construct_H0_HF(
    MeanField& meanfield,
    const std::map<int, std::map<int, Matz>> & H_KS_all,
    const std::map<int, std::map<int, Matz>> & vxc_all,
    const std::map<int, std::map<int, Matz>> & Hexx_all,
    int n_spins, int n_kpoints, int n_states) {

    // 初始化 HF 哈密顿量矩阵
    std::map<int, std::map<int, Matz>> H0_HF_all;
    for (int ispin = 0; ispin < n_spins; ++ispin)
    {
        for (int ikpt = 0; ikpt < n_kpoints; ++ikpt)
        {
            Matz Hexx_ispin_ik = Hexx_all.at(ispin).at(ikpt);
            Matz H0_HF_spin_k = H_KS_all.at(ispin).at(ikpt) - vxc_all.at(ispin).at(ikpt) + Hexx_ispin_ik ;
            H0_HF_all[ispin][ikpt] = H0_HF_spin_k;   
        }
    }

    return H0_HF_all;
}

void diagonalize_and_store(MeanField& meanfield, const std::map<int, std::map<int, Matz>>& H0_GW_all,
                           int n_spins, int n_kpoints, int dimension)
{
    int nao = meanfield.get_n_bands();

    for (int ispin = 0; ispin < n_spins; ++ispin)
    {
        for (int ikpt = 0; ikpt < n_kpoints; ++ikpt)
        {
            // 取出相应的哈密顿量矩阵
            const auto &h = H0_GW_all.at(ispin).at(ikpt).copy();
            const std::string final_banner(90, '-');
            // 对角化哈密顿量，得到 QP 波函数在 KS 表象下的表示
            std::vector<double> w;
            Matz eigvec_KS;

            // // 打印哈密顿量矩阵
            // printf("%77s\n", final_banner.c_str());
            // printf("Hamiltonian matrix (h):\n");
            // for (int i = 0; i < dimension; ++i) {
            //     for (int j = 0; j < nao; ++j) {
            //         printf("%20.16f ", h(i, j).real());
            //     }
            //     printf("\n"); // 换行
            // }
            // printf("%77s\n", final_banner.c_str());

            eigsh(h, w, eigvec_KS);

            // // 打印本征值 w
            // printf("Eigenvalues (w):\n");
            // for (int i = 0; i < dimension; ++i) {
            //     printf("%20.16f\n", w[i]);
            // }
            // printf("%77s\n", final_banner.c_str());

            // // 打印本征向量矩阵 eigvec_KS
            // printf("Eigenvectors (eigvec_KS):\n");
            // for (int i = 0; i < dimension; ++i) {
            //     for (int j = 0; j < nao; ++j) {
            //         printf("%20.16f ", eigvec_KS(i, j).real());
            //     }
            //     printf("\n"); // 换行
            // }
            // printf("%77s\n", final_banner.c_str());


            // 将本征值存储到 MeanField 的 eskb 矩阵
            for (int ib = 0; ib < dimension; ++ib)
            {
                meanfield.get_eigenvals()[ispin](ikpt, ib) = w[ib];
            }
            
            Matz wfc(dimension, nao, MAJOR::COL);
            
            for (int ib = 0; ib < dimension; ++ib)
            {
                for (int iao = 0; iao < nao; iao++)
                {
                    wfc(ib, iao) = meanfield.get_eigenvectors0()[ispin][ikpt](ib, iao);
                }
            }
            auto eigvec_NAO = transpose(eigvec_KS) * wfc;
            
            // printf("%77s\n", final_banner.c_str());
            // printf("Eigenvectors2:\n");
            // for (int i = 0; i < meanfield.get_n_bands(); i++) {
            //     for (int j = 0; j < meanfield.get_n_bands(); j++) {
            //         const auto &eigenvectors = meanfield.get_eigenvectors()[ispin][ikpt](i, j) ;
            //         printf("%20.16f ", eigenvectors.real());
            //     }
            //     printf("\n"); // 换行
            // }
            // printf("%77s\n", final_banner.c_str());
            // printf("\n");
            // 将 KS 表示旋转到 NAO 表示

            for (int ib = 0; ib < dimension; ++ib)
            {
                for (int iao = 0; iao < nao; iao++)
                {
                    meanfield.get_eigenvectors()[ispin][ikpt](ib, iao) = eigvec_NAO(ib, iao);
                }
            }
            // printf("%77s\n", final_banner.c_str());
            // printf("Eigenvectors3:\n");
            // for (int i = 0; i < meanfield.get_n_bands(); i++) {
            //     for (int j = 0; j < meanfield.get_n_bands(); j++) {
            //         const auto &eigenvectors = meanfield.get_eigenvectors()[ispin][ikpt](i, j) ;
            //         printf("%20.16f ", eigenvectors.real()); 
            //     }
            //     printf("\n"); // 换行
            // }
            // printf("%77s\n", final_banner.c_str());
            // printf("\n");
        }
    }
    std::cout << "所有本征值已存储到 MeanField 对象。" << std::endl;
}


// 
std::map<int, Matz> FT_K_TO_R(MeanField& meanfield, const std::map<int, Matz>& Vk_ispin, const std::vector<Vector3_Order<int>>& Rlist)
{    
    std::map<int, Matz> Vr_ispin;
    const auto n_aos = meanfield.get_n_aos();
    // 遍历实空间格点向量列表
    for (auto R : Rlist) 
    {  
        auto iteR = std::find(Rlist.cbegin(), Rlist.cend(), R);
        auto iR = std::distance(Rlist.cbegin(), iteR);   
        // std::cout << "iR 的值是: " << iR << std::endl;
        Vr_ispin[iR] = Matz(n_aos, n_aos, MAJOR::COL);  
        Vr_ispin[iR].zero_out();
        for (int i = 0; i < n_aos; i++){
            for (int j = 0; j < n_aos; j++){
                complex<double> element = 0.0;
                for (int i_kpoint = 0; i_kpoint < meanfield.get_n_kpoints(); i_kpoint++)
                {
                    auto Vk_element = Vk_ispin.at(i_kpoint)(i,j);
                    auto k = kfrac_list[i_kpoint];
                    double ang = - k * R * TWO_PI ;
                    // double ang_check = - k * R ;    
                    // std::cout << "ang 的值是: " << ang_check << std::endl;
                    complex<double> kphase = complex<double>(cos(ang), sin(ang));
                    element += Vk_element * kphase;
                    // element += 1.0 * kphase;
                    
                }
                Vr_ispin[iR](i,j)=element;
                // for (int i = 0; i < n_aos; ++i){
                //     for (int j = 0; j < n_aos; ++j){
                //         Vr_ispin[iR](i,j)=std::real(Vr_ispin[iR](i,j));
                //     }
                // }
            }
        }
    }
    return Vr_ispin;
}




std::map<int, Matz> FT_R_TO_K(MeanField& meanfield, const std::map<int, Matz>& Vr_ispin, const std::vector<Vector3_Order<int>>& Rlist)
{    
    std::map<int, Matz> Vk_ispin;
    const auto n_aos = meanfield.get_n_aos();
    // 遍历实空间格点向量列表
    for (int i_kpoint = 0; i_kpoint < meanfield.get_n_kpoints(); i_kpoint++)
    {  
        Vk_ispin[i_kpoint] = Matz(n_aos, n_aos, MAJOR::COL);  
        Vk_ispin[i_kpoint].zero_out();
        auto k = kfrac_list[i_kpoint];
        for (int i = 0; i < n_aos; i++){
            for (int j = 0; j < n_aos; j++){
                complex<double> element = 0.0;
                for (auto R : Rlist) 
                {                      
                    auto iteR = std::find(Rlist.cbegin(), Rlist.cend(), R);
                    auto iR = std::distance(Rlist.cbegin(), iteR);   
                    auto Vr_element = Vr_ispin.at(iR)(i,j);
                    double ang =  k * R * TWO_PI ;
                    complex<double> kphase = complex<double>(cos(ang), sin(ang));
                    element += Vr_element * kphase / meanfield.get_n_kpoints() ;
                    // element += 1.0 * kphase;
                }
                Vk_ispin[i_kpoint](i,j)=element;
            }
        }
    }
    return Vk_ispin;
}

//interation diag

Matz calculate_scRPA_xc_g_matrix(
    const std::vector<double>& freq_nodes,
    const std::vector<double>& freq_weights,
    const Matz& sigx_spin_k, 
    const std::map<double, Matz>& sigc_spin_k, 
    const std::vector<std::vector<cplxdb>>& G0, 
    int n_states ) 
{
    Matz V_rpa_ks(n_states, n_states, MAJOR::COL); 
    for (int n = 0; n < n_states; ++n) {
        for (int m = 0; m < n_states; ++m) {
            cplxdb V_nm = 0.0;
            for (size_t w = 0; w < freq_weights.size(); ++w) {
                cplxdb sig_nm_iw = sigc_spin_k.at(freq_nodes[w])(n, m)+sigx_spin_k(n, m);
                V_nm += freq_weights[w] * sig_nm_iw * G0[n][w] ;
            }
            V_rpa_ks(n, m) = V_nm ;            
        }
    }
    return V_rpa_ks;
}

Matz calculate_scRPA_lambda_matrix(
    MeanField& meanfield,
    const Matz& H_KS0,
    const Matz& vxc0,
    const Matz& gxc_matrix_diger,
    double mu, 
    double temperature,
    int ispin,
    int ikpt,
    int n_states)
{
    Matz lambda_matrix(n_states, n_states, MAJOR::COL);
    for (int i = 0; i < n_states; ++i)
    {
        for (int j = 0; j < n_states; ++j)
        {
            double energy_n = meanfield.get_eigenvals()[ispin](ikpt, j);
            double f_n = fermi_dirac(energy_n, mu, temperature) * 2.0 / meanfield.get_n_spins();
            std::complex<double> Vc_ij = f_n *(H_KS0(i, j) - vxc0(i, j)) + gxc_matrix_diger(i, j);
            lambda_matrix(i, j) = Vc_ij;  
        }
    }    
    return lambda_matrix;
}

Matz construct_F_matrix(const Matz& lambda_matrix,std::vector<double>& F_eigenvalues,int n_states,int iteration)
{
    Matz F_matrix(n_states, n_states, MAJOR::COL);
    Matz lambda_matrix_diger = transpose(lambda_matrix,true);
    std::vector<double> w;
    if(iteration==1){
        Matz Re_lambda_matrix(n_states, n_states, MAJOR::COL);        
        Matz eigvec;
        Re_lambda_matrix = 0.5 * (lambda_matrix + lambda_matrix_diger);
        eigsh(Re_lambda_matrix, w, eigvec);
    }
    else{
        w=F_eigenvalues;
    }
    
    
    for (int i = 0; i < n_states; ++i)
    {
        F_matrix(i, i) = w[i];  
        for (int j = 0; j < n_states; ++j)
        {
            if(i<j){
                std::complex<double> F_ij = lambda_matrix(i, j) - lambda_matrix_diger(i, j);
                F_matrix(i, j) = F_ij;  
            }
            if(i>j){
                std::complex<double> F_ij = lambda_matrix_diger(i, j) - lambda_matrix(i, j);
                F_matrix(i, j) = F_ij;  
            }
        }
    }   

    if(F_matrix==transpose(F_matrix,true)){
        std::cout << "F_matrix is hermitian" << std::endl;
    }
    else{
        std::cout << "F_matrix is not hermitian" << std::endl;
    }
    return F_matrix;
}

std::vector<double> solve_update_F_matrix(MeanField& meanfield,const Matz& F_matrix,int ispin,int ikpt,int n_states)
{
    std::vector<double> w;
    Matz eigvec_KS;
    eigsh(F_matrix, w, eigvec_KS);
    std::vector<double> F_eigenvalues;
    for (int i = 0; i < n_states; ++i)
    {
        F_eigenvalues.push_back(w[i]);
    }

    Matz wfc(n_states, n_states, MAJOR::COL);
            
    for (int ib = 0; ib < n_states; ++ib)
    {
        for (int iao = 0; iao < n_states; iao++)
        {
            wfc(ib, iao) = meanfield.get_eigenvectors0()[ispin][ikpt](ib, iao);
        }
    }
    auto eigvec_NAO = transpose(eigvec_KS) * wfc;

    for (int ib = 0; ib < n_states; ++ib)
    {
        for (int iao = 0; iao < n_states; iao++)
        {
            meanfield.get_eigenvectors()[ispin][ikpt](ib, iao) = eigvec_NAO(ib, iao);
        }
    }




    return F_eigenvalues;
}