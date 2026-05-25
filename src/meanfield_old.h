/*!
 @file meanfield.h
 @brief Utilities to handle the mean-field starting point for many-body calculation
 */
#ifndef MEANFIELD_H
#define MEANFIELD_H

#include <map>
#include <vector>

#include "math/complexmatrix.h"
#include "math/matrix.h"
#include "math/vector3_order.h"
#include "mpi/base_mpi.h"

using librpa_int::matrix;
using librpa_int::ComplexMatrix;
using librpa_int::Vector3_Order;
using librpa_int::MpiCommHandler;

//! Object of the meanfield input of Green's function
/*!
  @note Energies are saved in Rydberg unit.
        Currently only n_spins = 1 is implemented.
 */
class MeanField
{
private:
    //! number of spin channels
    int n_spins;
    //! number of atomic orbitals
    int n_aos;
    //! number of bands
    int n_bands;
    //! number of kpoints
    int n_kpoints;
    //! number of up down when SOC, n_soc = 1 for non-SOC, n_soc = 2 for SOC.
    //! all physical quantities related to n_aos should be refactored.
    int n_soc;
    //! eigenvalues, (n_spins, n_kpoints, n_bands)
    std::vector<matrix> eskb;
    //! occupation weight, scaled by n_kpoints. (n_spins, n_kpoints, n_bands)
    std::vector<matrix> wg;
    std::vector<matrix> wg0;
    //! eigenvector, (n_spins, n_soc, n_kpoint, n_bands, n_aos)
    std::vector<std::vector<std::vector<ComplexMatrix>>> wfc;
    std::vector<std::vector<std::vector<ComplexMatrix>>> wfc0;
    //! unit: eV*m , velocity_matrix, (n_spins, n_kpoint, n_alpha, n_bands, n_bands)
    std::vector<std::vector<std::vector<ComplexMatrix>>> velocity;
    //! Fermi energy
    double efermi;
    void resize(int ns, int nk, int nb, int nao);

public:
    MeanField() : n_spins(0), n_aos(0), n_bands(0), n_kpoints(0) {};
    MeanField(int ns, int nk, int nb, int nao);
    ~MeanField() {};
    void set(int ns, int nk, int nb, int nao);
    MeanField(const MeanField&);
    inline int get_n_bands() const { return n_bands; }
    inline int get_n_spins() const { return n_spins; }
    inline int get_n_kpoints() const { return n_kpoints; }
    inline int get_n_aos() const { return n_aos; }
    inline int get_n_soc() const { return n_soc; }
    inline double& get_efermi() { return efermi; }
    inline const double& get_efermi() const { return efermi; }
    std::vector<matrix>& get_eigenvals() { return eskb; }
    const std::vector<matrix>& get_eigenvals() const { return eskb; }
    std::vector<matrix>& get_weight() { return wg; }
    double get_total_weight() const;
    const std::vector<matrix>& get_weight() const { return wg; }
    std::vector<matrix> & get_weight0() { return wg0; }
    const std::vector<matrix> & get_weight0() const { return wg0; }
    //! get the density matrix of a particular spin and kpoint
    ComplexMatrix get_dmat_cplx(int ispin, int isoc1, int isoc2, int ikpt) const;
    ComplexMatrix get_dmat_cplx_R(int ispin, int isoc1, int isoc2,
                                  const std::vector<Vector3_Order<double>>& kfrac_list,
                                  const Vector3_Order<int>& R) const;
    std::vector<std::vector<std::vector<ComplexMatrix>>>& get_eigenvectors() { return wfc; }
    std::vector<std::vector<std::vector<ComplexMatrix>>>& get_eigenvectors0() { return wfc0; }
    const std::vector<std::vector<std::vector<ComplexMatrix>>>& get_eigenvectors() const
    {
        return wfc;
    }
    const std::vector<std::vector<std::vector<ComplexMatrix>>>& get_eigenvectors0() const
    {
        return wfc0;
    }
    std::vector<std::vector<std::vector<ComplexMatrix>>>& get_velocity() { return velocity; }
    double get_E_min_max(double& emin, double& emax) const;
    double get_band_gap() const;
    std::map<double, std::map<Vector3_Order<int>, ComplexMatrix>> get_gf_cplx_imagtimes_Rs(
        int ispin, int isoc1, int isoc2, const std::vector<Vector3_Order<double>>& kfrac_list,
        std::vector<double> imagtimes, const std::vector<Vector3_Order<int>>& Rs) const;
    std::map<double, std::map<Vector3_Order<int>, matrix>> get_gf_real_imagtimes_Rs(
        int ispin, int isoc1, int isoc2, const std::vector<Vector3_Order<double>>& kfrac_list,
        std::vector<double> imagtimes, const std::vector<Vector3_Order<int>>& Rs) const;
    void allredue_wfc_isk();
    void broadcast(const MpiCommHandler& comm_hdl, int root);
};

//! A global MeanField object
extern MeanField meanfield;
extern MeanField pyatb_meanfield;

#endif  // ! MEANFIELD_H
