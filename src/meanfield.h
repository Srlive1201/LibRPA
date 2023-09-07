/*!
 @file meanfield.h
 @brief Utilities to handle the mean-field starting point for many-body calculation
 */
#ifndef MEANFIELD_H
#define MEANFIELD_H

#include <vector>
#include <map>
#include "matrix.h"
#include "complexmatrix.h"
#include "vector3_order.h"
using std::vector;

//! Object of the meanfield input of Green's function
/*!
  \note Energies are saved in Rydberg unit.
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
        //! eigenvalues, (n_spins, n_kpoints, n_bands)
        vector<matrix> eskb;
	//std::map<int, std::map<int, matrix>> d_eskb;
	map<int, map<int, matrix>> d_eskb;
        //! occupation weight, scaled by n_kpoints. (n_spins, n_kpoints, n_bands)
        vector<matrix> wg;
        //! eigenvector, (n_spins, n_kpoint, n_bands, n_aos)
        vector<vector<ComplexMatrix>> wfc;
        //! Fermi energy
        double efermi;
        void resize(int ns, int nk, int nb, int nao);
    public:
        MeanField():
            n_spins(0), n_aos(0), n_bands(0), n_kpoints(0)
            {};
        MeanField(int ns, int nk, int nb, int nao);
        ~MeanField() {};
        void set(int ns, int nk, int nb, int nao);
        MeanField(const MeanField &);
        int get_n_bands() const { return n_bands; }
        int get_n_spins() const { return n_spins; }
        int get_n_kpoints() const { return n_kpoints; }
        int get_n_aos() const { return n_aos; }
        double & get_efermi() { return efermi; }
        const double & get_efermi() const { return efermi; }
        vector<matrix> & get_eigenvals() { return eskb; }
        map<int, map<int, matrix>> & get_d_eigenvals() { return d_eskb; }
        const vector<matrix> & get_eigenvals() const { return eskb; }
        vector<matrix> & get_weight() { return wg; }
        const vector<matrix> & get_weight() const { return wg; }
        //! get the density matrix of a particular spin and kpoint
        ComplexMatrix get_dmat_cplx(int ispin, int ikpt) const;
        ComplexMatrix get_dmat_cplx_R(int ispin, const std::vector<Vector3_Order<double>>& kfrac_list, const Vector3_Order<int>& R) const;
        vector<vector<ComplexMatrix>> & get_eigenvectors() { return wfc; }
        const vector<vector<ComplexMatrix>> & get_eigenvectors() const { return wfc; }
        double get_E_min_max(double &emin, double &emax);
        double get_band_gap();
};

//! A global MeanField object
extern MeanField meanfield;

#endif // ! MEANFIELD_H
