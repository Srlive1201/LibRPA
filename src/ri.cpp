#include "ri.h"
#include <memory.h>

int n_irk_points;
int natom;
int ncell;
map<Vector3_Order<double>, double> irk_weight;
map<atom_t, size_t> atom_nw;
map<atom_t, size_t> atom_mu;

atpair_R_mat_t Cs;
atpair_k_cplx_mat_t Vq;

int atom_iw_loc2glo(const int &atom_index, const int &iw_lcoal)
{
    int nb = 0;
    for (int ia = 0; ia != atom_index; ia++)
        nb += atom_nw[ia];
    return iw_lcoal + nb;
}

int atom_mu_loc2glo(const int &atom_index, const int &mu_lcoal)
{
    int nb = 0;
    for (int ia = 0; ia != atom_index; ia++)
        nb += atom_mu[ia];
    return mu_lcoal + nb;
}

vector<int> get_part_range()
{
    vector<int> part_range;
    part_range.resize(atom_mu.size());
    part_range[0] = 0;

    int count_range = 0;
    for (int Mu = 0; Mu != atom_mu.size() - 1; Mu++)
    {
        count_range += atom_mu[Mu];
        part_range[Mu + 1] = count_range;
    }
    return part_range;
}

matrix reshape_Cs(size_t n1, size_t n2, size_t n3, const shared_ptr<matrix> &Csmat) //(n1*n2,n3) -> (n2,n1*n3)
{
    const auto length = sizeof(double) * n3;
    const auto n13 = n1 * n3;
    const double *m_ptr = (*Csmat).c;
    matrix m_new(n2, n1 * n3, false);
    for (size_t i1 = 0; i1 != n1; ++i1)
    {
        double *m_new_ptr = m_new.c + i1 * n3;
        for (size_t i2 = 0; i2 != n2; ++i2, m_ptr += n3, m_new_ptr += n13)
            memcpy(m_new_ptr, m_ptr, length);
    }
    return m_new;
}

matrix reshape_mat(const size_t n1, const size_t n2, const size_t n3, const matrix &mat) //(n1,n2*n3) -> (n2,n1*n3)
{
    const auto length = sizeof(double) * n3;
    const auto n13 = n1 * n3;
    const double *m_ptr = mat.c;
    matrix m_new(n2, n1 * n3, false);
    for (size_t i1 = 0; i1 != n1; ++i1)
    {
        double *m_new_ptr = m_new.c + i1 * n3;
        for (size_t i2 = 0; i2 != n2; ++i2, m_ptr += n3, m_new_ptr += n13)
            memcpy(m_new_ptr, m_ptr, length);
    }
    return m_new;
}

matrix reshape_mat_21(const size_t n1, const size_t n2, const size_t n3, const matrix &mat) //(n1,n2*n3) -> (n1*n2,n3)
{
    const auto length = sizeof(double) * n1 * n2 * n3;
    const auto n13 = n1 * n2 * n3;
    const double *m_ptr = mat.c;
    matrix m_new(n1 * n2, n3, false);
    double *m_new_ptr = m_new.c;
    memcpy(m_new_ptr, m_ptr, length);
    return m_new;
}
