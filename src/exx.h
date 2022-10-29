/*
 * @file exx.h
 * @brief utilities for computing exact exchange energies, including orbital and total energies.
 */
#include "meanfield.h"
#include "ri.h"

// NOTE: a temporary namespace to avoid redefinition with LibRI
namespace LIBRPA
{

class Exx
{
    private:
        //! refenrence to the MeanField object to compute density matrix
        const MeanField& mf_;
        //! reference to the fractional kpoint list on which the MeanField object is computed
        const vector<Vector3_Order<double>>& kfrac_list_;

        ComplexMatrix get_dmat_cplx_R_global(const int& ispin, const Vector3_Order<int>& R);
        ComplexMatrix extract_dmat_cplx_R_IJblock(const ComplexMatrix& dmat_cplx, const atom_t& I, const atom_t& J);
        void build_dmat_R(const Vector3_Order<int>& R);
        void build_dmat_R(const atom_t& I, const atom_t& J, const Vector3_Order<int>& R);
        void warn_dmat_IJR_nonzero_imag(const ComplexMatrix& dmat_cplx,
                                        const int& ispin, const atom_t& I, const atom_t& J,
                                        const Vector3_Order<int> R);
        void build_exx_orbital_energy_LibRI(const atpair_R_mat_t& LRI_Cs,
                                            const vector<Vector3_Order<int>> &Rlist,
                                            const Vector3_Order<int> &R_period,
                                            const atpair_R_mat_t& coul_mat);
    public:
        //! Density matrix in lattice vector space, dimension (nspins, I, J, R, nao_I, nao_J)
        map<int, atpair_R_mat_t> dmat;
        //! exact-exchange Hamiltonian in k space, dimension (nspins, I, J, k, nao_I, nao_J)
        map<int, atpair_k_cplx_mat_t> Hexx;

        //! exact-exchange Hamiltonian in the basis of KS states, dimension (nspins, n_kpoints, n_bands, n_bands)
        map<int, map<int, ComplexMatrix>> Hexx_KS;

        //! exact-exchange energy of each state, dimension (nspins, n_kpoints, n_bands). This is actually the diagonal elements of Heex_KS.
        map<int, map<int, map<int, double>>> Eexx;

        Exx(const MeanField& mf, const vector<Vector3_Order<double>> &kfrac_list): mf_(mf), kfrac_list_(kfrac_list) {};
        //! Build and store the density matrix from the meanfield object
        void build_exx_orbital_energy(const atpair_R_mat_t& LRI_Cs,
                                      const vector<Vector3_Order<int>> &Rlist,
                                      const Vector3_Order<int> &R_period,
                                      const atpair_R_mat_t& coul_mat);
};

} // namespace LIBRPA
