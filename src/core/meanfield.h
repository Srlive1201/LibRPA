/*!
 @file meanfield.h
 @brief Utilities to handle the mean-field starting point for many-body calculation
 */
#ifndef MEANFIELD_H
#define MEANFIELD_H

#include <vector>
#include <map>
#include "../math/matrix.h"
#include "../math/complexmatrix.h"
#include "../math/vector3_order.h"

namespace librpa_int {

//! Object of the meanfield input
/*!
  @note Energies are saved in Hartree unit.
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

    // Relevant local dimensions within parallel distribution.
    // Only for saving eigenvectors
    int n_aos_local;
    int n_bands_local;
    int n_kpoints_local;
    std::vector<int> iks_local;

    //! eigenvalues, (n_spins, n_kpoints, n_bands)
    std::vector<matrix> eskb;
    //! occupation weight, scaled by n_kpoints. (n_spins, n_kpoints, n_bands)
    std::vector<matrix> wg;
    //! eigenvector, (n_spins, n_kpoint, n_bands, n_aos)
    // TODO: parallelize
    std::vector<std::vector<ComplexMatrix>> wfc;
    //! Fermi energy
    double efermi;
    void resize(int ns, int nk, int nb, int nao, int nk_local, int nb_local, int nao_local);
public:
    MeanField()
        : n_spins(0),
        n_aos(0),
        n_bands(0),
        n_kpoints(0),
        n_aos_local(0),
        n_bands_local(0),
        n_kpoints_local(0),
        iks_local(),
        eskb(),
        wg(),
        wfc(),
        efermi(0) {};
    MeanField(int ns, int nk, int nb, int nao);
    MeanField(int ns, int nk, int nb, int nao, int nk_local, int nb_local, int nao_local);
    ~MeanField() {};
    void set(int ns, int nk, int nb, int nao);
    void set(int ns, int nk, int nb, int nao, int nk_local, int nb_local, int nao_local);
    MeanField(const MeanField&);
    inline int get_n_bands() const { return n_bands; }
    inline int get_n_states() const { return n_bands; } // alias
    inline int get_n_spins() const { return n_spins; }
    inline int get_n_kpoints() const { return n_kpoints; }
    inline int get_n_aos() const { return n_aos; }
    inline double& get_efermi() { return efermi; }
    inline const double& get_efermi() const { return efermi; }
    std::vector<matrix>& get_eigenvals() { return eskb; }
    const std::vector<matrix>& get_eigenvals() const { return eskb; }
    std::vector<matrix>& get_weight() { return wg; }
    const std::vector<matrix>& get_weight() const { return wg; }
    //! get the density matrix of a particular spin and kpoint
    ComplexMatrix get_dmat_cplx(int ispin, int ikpt) const;
    ComplexMatrix get_dmat_cplx_R(int ispin, const std::vector<Vector3_Order<double>>& kfrac_list,
                                const Vector3_Order<int>& R) const;
    std::vector<std::vector<ComplexMatrix>>& get_eigenvectors() { return wfc; }
    const std::vector<std::vector<ComplexMatrix>>& get_eigenvectors() const { return wfc; }
    double get_E_min_max(double& emin, double& emax) const;
    double get_band_gap() const;
    std::map<double, std::map<Vector3_Order<int>, ComplexMatrix>> get_gf_cplx_imagtimes_Rs(
        int ispin, const std::vector<Vector3_Order<double>>& kfrac_list,
        std::vector<double> imagtimes, const std::vector<Vector3_Order<int>>& Rs) const;
    std::map<double, std::map<Vector3_Order<int>, matrix>> get_gf_real_imagtimes_Rs(
        int ispin, const std::vector<Vector3_Order<double>>& kfrac_list,
        std::vector<double> imagtimes, const std::vector<Vector3_Order<int>>& Rs) const;
    // void allredue_wfc_isk();
};

}

#endif // ! MEANFIELD_H
