#pragma once
#include "meanfield.h"
double calculate_total_occupation(const MeanField &mf, double mu, double temperature);
double calculate_fermi_energy(const MeanField &mf, double temperature, double total_electrons);
double calculate_eqp_fermi_energy(const MeanField &mf,
                                  std::map<int, std::map<int, std::map<int, double>>> e_qp_all, 
                                  double temperature, 
                                  double total_electrons) ;
double fermi_dirac(double energy, double mu, double temperature);
void update_fermi_energy_and_occupations(MeanField &meanfield, const double temperature, const double efermi);
