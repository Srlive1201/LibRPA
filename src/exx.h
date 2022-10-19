/*
 * @file exx.h
 * @brief utilities for computing exact exchange energies, including orbital and total energies.
 */
#include "meanfield.h"
#include "ri.h"

class Exx
{
    private:
        double dmat_threshold;
        //! Density matrix in lattice vector space, dimension (nspins, I, J, R, nao_I, nao_J)
        map<int, atpair_R_mat_t> dmat;
#ifdef __USE_LIBRI
#endif
    public:
        Exx();
        //! Build and store the density matrix from the meanfield object
        void build_density_matrix_R(const MeanField& mf,
                                    const vector<Vector3_Order<double>> &klist,
                                    const Vector3_Order<int>& R);
        void build_exx_orbital_energy(const MeanField& mf,
                                      const vector<Vector3_Order<double>> &klist,
                                      const atpair_R_mat_t& LRI_Cs,
                                      const vector<Vector3_Order<int>> &Rlist,
                                      const Vector3_Order<int> &R_period,
                                      const atpair_R_cplx_mat_t& coul_mat);
};
