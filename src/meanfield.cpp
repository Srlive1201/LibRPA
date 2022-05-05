#include "meanfield.h"
#include <stdexcept>
#include <iostream>

void MeanField::resize(int ns, int nk, int nb, int nao)
{
    if (ns == 0 || nk == 0 || nb == 0 || nao == 0)
        throw invalid_argument("encounter zero dimension");
    if (n_spins != 0)
    {
        eskb.clear();
        wg.clear();
        wfc.clear();
    }

    n_spins = ns;
    n_kpoints = nk;
    n_bands = nb;
    n_aos = nao;

    eskb.resize(n_spins);
    wg.resize(n_spins);
    wfc.resize(n_spins);

    for (int is = 0; is < n_spins; is++)
    {
        eskb[is].create(n_kpoints, n_bands);
        wg[is].create(n_kpoints, n_bands);
        wfc[is].resize(n_kpoints);
        for (int ik = 0; ik < n_kpoints; ik++)
            wfc[is][ik].create(n_bands, n_aos);
    }
}

MeanField::MeanField(int ns, int nk, int nb, int nao)
{
    resize(ns, nk, nb, nao);
}

void MeanField::set(int ns, int nk, int nb, int nao)
{
    if (n_spins != 0 || n_kpoints != 0 || n_bands != 0 || n_aos != 0)
    {
        std::cout << n_spins << n_kpoints << n_bands << n_aos << std::endl;
        throw invalid_argument("MeanField object already set");
    }
    resize(ns, nk, nb, nao);
}

MeanField::MeanField(const MeanField &mf)
{
    resize(mf.n_spins, mf.n_kpoints, mf.n_bands, mf.n_aos);
    // FIXME: copy data, not tested
    eskb = mf.eskb;
    wg = mf.wg;
    wfc = mf.wfc;
    efermi = mf.efermi;
}

double MeanField::get_E_min_max(double &emin, double &emax)
{
    double lb = eskb[0](0, 0);
    double ub = eskb[0](0, n_bands - 1);
    for (int is = 0; is != n_spins; is++)
        for (int ik = 0; ik != n_kpoints; ik++)
        {
            lb = (lb > eskb[is](ik, 0)) ? eskb[is](ik, 0) : lb;
            ub = (ub < eskb[is](ik, n_bands-1)) ? eskb[is](ik, n_bands-1) : ub;
        }
    double gap = get_band_gap();
    emax = 0.5 * (ub - lb);
    emin = gap;
    return emin;
}

double MeanField::get_band_gap()
{
    double gap = 1e6;
    double homo = -1e6, lumo = 1e6;
    // FIXME: should be nspins/nkpoints?
    double midpoint = 1.0 / (n_spins * n_kpoints);
    for (int is = 0; is != n_spins; is++)
        for (int ik = 0; ik != n_kpoints; ik++)
        {
            int homo_level = 0;
            for (int n = 0; n != n_bands; n++)
            {
                if (wg[is](ik, n) >= midpoint)
                {
                    homo_level = n;
                }
            }
            if (homo < eskb[is](ik, homo_level))
                homo = eskb[is](ik, homo_level);

            if (lumo > eskb[is](ik, homo_level + 1))
                lumo = eskb[is](ik, homo_level + 1);
        }
    gap = 0.5 * (lumo - homo);
    return gap;
}

MeanField meanfield = MeanField();