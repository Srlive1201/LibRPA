#include "fermi_energy_occupation.h"

#include <cmath>
#include <vector>
#include "constants.h"

// 费米分布
double fermi_dirac(double energy, double mu, double temperature)
{
    // const double K_B = 3.16681e-6;  // Hartree/K
    // return 1.0 / (1.0 + exp((energy - mu) / (K_B * temperature)));
    if (energy <= mu) {
        return 1.0;
    } else {
        return 0.0;
    }
}

// 计算给定化学势下的总占据态
double calculate_total_occupation(const MeanField &mf, double mu, double temperature) {
    double total_occupation = 0.0;

    for (int ispin = 0; ispin < mf.get_n_spins(); ++ispin) {
        for (int ikpt = 0; ikpt < mf.get_n_kpoints(); ++ikpt) {
            for (int ib = 0; ib < mf.get_n_bands(); ++ib) {
                double energy = mf.get_eigenvals()[ispin](ikpt, ib);
                double occupation = fermi_dirac(energy, mu, temperature) * 2.0 / (mf.get_n_kpoints() * mf.get_n_spins());
                total_occupation += occupation;
            }
        }
    }

    return total_occupation;
}
//calculate_local_occupation_for each ispin-ikpoint_0-temperature
static double calculate_local_occupation(const MeanField &mf, double mu, double temperature, int ispin, int ikpt) {
    double local_occupation = 0.0;
    for (int ib = 0; ib < mf.get_n_bands(); ++ib) {
        double energy = mf.get_eigenvals()[ispin](ikpt, ib);
        double occupation = fermi_dirac(energy, mu, temperature) * 2.0 / mf.get_n_spins();
        local_occupation += occupation;
    }
    return local_occupation;
}

// //1
// double calculate_fermi_energy(const MeanField &mf, double temperature, double total_electrons) {
//     double tolerance = 1e-4;  
//     double total_occupation = 0.0;
//     double mu = 0.0;
//     double gap = 0.0;
//     double vbm = -10000.0;  // 比mu小的最大值
//     double cbm = 10000.0;  // 比mu大的最小值

//     for (int ispin = 0; ispin < mf.get_n_spins(); ++ispin) {
//         for (int ikpt = 0; ikpt < mf.get_n_kpoints(); ++ikpt) {
            
            
//             for (int ib = 0; ib < mf.get_n_bands(); ++ib) {
//                 double energy = mf.get_eigenvals()[ispin](ikpt, ib);
//                 mu = energy;
//                 total_occupation = calculate_total_occupation(mf, mu, temperature);
//                 // 如果能量比mu小，更新vbm
//                 if (total_occupation < total_electrons + tolerance) {
//                     if (energy > vbm) {
//                         vbm = energy;
//                         std::cout << "ikpt1: " << ikpt << std::endl;
//                         std::cout << "ib1: " << ib << std::endl;
//                     } 
//                 }
//                 // 如果能量比mu大，更新cbm
//                 else{
//                     if (energy < cbm) {
//                         cbm = energy;
//                         std::cout << "ikpt2: " << ikpt << std::endl;
//                         std::cout << "ib2: " << ib << std::endl;
//                     } 
//                 }
//             }
//         }
//     }
//     // 最终费米能级取 vbm 和 cbm 的中间值
//     mu = (vbm + cbm) * 0.5;
//     gap = cbm - vbm ;
//     std::cout << "Final VBM: " << vbm * HA2EV<< ", CBM: " << cbm * HA2EV << ", Final Fermi level: " << mu * HA2EV << std::endl;
//     std::cout << "Hamiltonian_gap: " << gap * HA2EV << " eV, "<< std::endl;
//     return mu;  
// }

//calculate_semiconductor_gap
double calculate_fermi_energy(const MeanField &mf, double temperature, double total_electrons) {
    double tolerance = 1e-5;   
    double mu = 0.0;
    double gap = 0.0;
    double vbm = -10000.0;  // 比mu小的最大值
    double cbm = 10000.0;  // 比mu大的最小值

    for (int ispin = 0; ispin < mf.get_n_spins(); ++ispin) {
        for (int ikpt = 0; ikpt < mf.get_n_kpoints(); ++ikpt) {
            double local_vbm = -10000.0;
            double local_cbm = 10000.0;
            double local_occupation = 0.0;
            for (int ib = 0; ib < mf.get_n_bands(); ++ib) {
                double local_mu = mf.get_eigenvals()[ispin](ikpt, ib);
                local_occupation = calculate_local_occupation(mf, local_mu, temperature, ispin, ikpt);
                
                // update vbm,cbm;
                if (local_occupation <= total_electrons + tolerance ) {
                    local_vbm = local_mu;                    
                }
                if (local_occupation > total_electrons + tolerance ) {
                    if (local_mu < local_cbm) {
                        local_cbm = local_mu;
                    }   
                }
            }
            if (local_vbm > vbm) {
                vbm = local_vbm;
            } 
            if (local_cbm < cbm){
                cbm = local_cbm;
            }
        }
    }
    mu = (vbm + cbm) * 0.5;
    gap = cbm - vbm ;
    std::cout << "Final VBM: " << vbm * HA2EV<< ", CBM: " << cbm * HA2EV << ", Final Fermi level: " << mu * HA2EV << std::endl;
    std::cout << "Hamiltonian_gap: " << gap * HA2EV << " eV, "<< std::endl;
    return mu;  
}


//calculate_local_occupation_for each ispin-ikpoint_0-temperature
static double calculate_eqp_local_occupation(const MeanField &mf, std::map<int, std::map<int, std::map<int, double>>> e_qp_all, double mu, double temperature, int ispin, int ikpt) {
    double local_occupation = 0.0;
    for (int ib = 0; ib < mf.get_n_bands(); ++ib) {
        double energy = e_qp_all[ispin][ikpt][ib];
        double occupation = fermi_dirac(energy, mu, temperature) * 2.0 / mf.get_n_spins();
        local_occupation += occupation;
    }
    return local_occupation;
}

static double calculate_eqp_total_occupation(const MeanField &mf, std::map<int, std::map<int, std::map<int, double>>> e_qp_all, double mu, double temperature) {
    double total_occupation = 0.0;

    for (int ispin = 0; ispin < mf.get_n_spins(); ++ispin) {
        for (int ikpt = 0; ikpt < mf.get_n_kpoints(); ++ikpt) {
            for (int ib = 0; ib < mf.get_n_bands(); ++ib) {
                double energy = e_qp_all[ispin][ikpt][ib];
                double occupation = fermi_dirac(energy, mu, temperature) * 2.0 / (mf.get_n_kpoints() * mf.get_n_spins());
                total_occupation += occupation;
            }
        }
    }

    return total_occupation;
}
// //2
// double calculate_eqp_fermi_energy(const MeanField &mf,
//                                   std::map<int, std::map<int, std::map<int, double>>> e_qp_all, 
//                                   double temperature, 
//                                   double total_electrons) {
//     double tolerance = 1e-4;  
//     double total_occupation = 0.0;
//     double mu = 0.0;
//     double gap = 0.0;
//     double vbm = -10000.0;  // 比mu小的最大值
//     double cbm = 10000.0;  // 比mu大的最小值

//     for (int ispin = 0; ispin < mf.get_n_spins(); ++ispin) {
//         for (int ikpt = 0; ikpt < mf.get_n_kpoints(); ++ikpt) {
            
            
//             for (int ib = 0; ib < mf.get_n_bands(); ++ib) {
//                 double energy = e_qp_all[ispin][ikpt][ib];
//                 mu = energy;
//                 total_occupation = calculate_eqp_total_occupation(mf, e_qp_all, mu, temperature);
//                 // 如果能量比mu小，更新vbm
//                 if (total_occupation < total_electrons + tolerance) {
//                     if (energy > vbm) {
//                         vbm = energy;
//                         std::cout << "ikpt3: " << ikpt << std::endl;
//                         std::cout << "ib3: " << ib << std::endl;
//                     } 
//                 }
//                 // 如果能量比mu大，更新cbm
//                 else{
//                     if (energy < cbm) {
//                         cbm = energy;
//                         std::cout << "ikpt4: " << ikpt << std::endl;
//                         std::cout << "ib4: " << ib << std::endl;
//                     } 
//                 }
//             }
//         }
//     }

//     // 最终费米能级取 vbm 和 cbm 的中间值
//     mu = (vbm + cbm) * 0.5;
//     gap = cbm - vbm ;
//     std::cout << "Final eqp_VBM: " << vbm* HA2EV << ", eqp_CBM: " << cbm* HA2EV << ", Final eqp_Fermi level: " << mu * HA2EV<< std::endl;
//     std::cout << "eqp_gap: " << gap * HA2EV << " eV, "<< std::endl;
//     return gap;  
// }
//22
double calculate_eqp_fermi_energy(const MeanField &mf,
                                  std::map<int, std::map<int, std::map<int, double>>> e_qp_all, 
                                  double temperature, 
                                  double total_electrons) {
    double tolerance = 1e-5;                                 
    double mu = 0.0;
    double gap = 0.0;
    double vbm = -10000.0;  // 比mu小的最大值
    double cbm = 10000.0;  // 比mu大的最小值

    for (int ispin = 0; ispin < mf.get_n_spins(); ++ispin) {
        for (int ikpt = 0; ikpt < mf.get_n_kpoints(); ++ikpt) {
            double local_vbm = -10000.0;
            double local_cbm = 10000.0;
            double local_occupation = 0.0;
            for (int ib = 0; ib < mf.get_n_bands(); ++ib) {
                double local_mu = e_qp_all[ispin][ikpt][ib];
                local_occupation = calculate_eqp_local_occupation(mf,e_qp_all ,local_mu, temperature, ispin, ikpt);
                std::cout << "local_occupation: " << local_occupation << std::endl;
                // update vbm,cbm;
                if (local_occupation <= total_electrons + tolerance ) {
                    local_vbm = local_mu;                    
                }
                if (local_occupation > total_electrons + tolerance) {
                    if (local_mu < local_cbm) {
                        local_cbm = local_mu;
                    }   
                }
            }
            if (local_vbm > vbm) {
                vbm = local_vbm;
            } 
            if (local_cbm < cbm){
                cbm = local_cbm;
            }
        }
    }

    // 最终费米能级取 vbm 和 cbm 的中间值
    mu = (vbm + cbm) * 0.5;
    gap = cbm - vbm ;
    std::cout << "Final eqp_VBM: " << vbm* HA2EV << ", eqp_CBM: " << cbm* HA2EV << ", Final eqp_Fermi level: " << mu * HA2EV<< std::endl;
    std::cout << "eqp_gap: " << gap * HA2EV << " eV, "<< std::endl;
    return gap;  
}





void update_fermi_energy_and_occupations(MeanField &mf, const double temperature, const double efermi)
{
    double total_electrons1 = 0.0;
    // 更新占据数
    for (int ispin = 0; ispin < mf.get_n_spins(); ++ispin)
    {
        for (int ikpt = 0; ikpt < mf.get_n_kpoints(); ++ikpt)
        {
            for (int ib = 0; ib < mf.get_n_bands(); ++ib)
            {
                const double energy = mf.get_eigenvals()[ispin](ikpt, ib);
                mf.get_weight()[ispin](ikpt, ib) = fermi_dirac(energy, efermi, temperature) * 2.0 / (mf.get_n_kpoints() * mf.get_n_spins());
                total_electrons1 += (mf.get_weight()[ispin](ikpt, ib)*mf.get_n_kpoints());  // 计算总占据数
            }
        }
    }
    total_electrons1 = total_electrons1 / mf.get_n_kpoints();
    // 输出 total_electrons
    std::cout << "Total electrons: " << total_electrons1 << std::endl;
    std::cout << "efermi: " << efermi << std::endl;
    mf.get_efermi() = efermi;
}
