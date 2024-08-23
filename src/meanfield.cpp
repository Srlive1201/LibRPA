#include "meanfield.h"
#include "lapack_connector.h"
#include "constants.h"
#include <stdexcept>
#include <iostream>
#include "envs_mpi.h"

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
    emax = ub - lb;
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
    {
        //print_matrix("mf.eskb: ",this->eskb[is]);
        for (int ik = 0; ik != n_kpoints; ik++)
        {
            int homo_level = -1;
            for (int n = 0; n != n_bands; n++)
            {
                if (wg[is](ik, n) >= midpoint)
                {
                    homo_level = n;
                }
            }
            //cout<<"|is ik: "<<is<<" "<<ik<<"  homo_level: "<<homo_level<<"   eskb0: "<<eskb[is](ik, homo_level)<<"  eskb1: "<<eskb[is](ik, homo_level + 1)<<endl;
            lumo = eskb[is](ik, homo_level + 1) < lumo ?  eskb[is](ik, homo_level + 1) : lumo;
            if(homo_level != -1)
                homo = eskb[is](ik, homo_level) > homo ? eskb[is](ik, homo_level) : homo;
            
            //cout<<"   homo: "<<homo<<"  lumo: "<<lumo<<endl;
        }
    }
    gap = lumo - homo;
    return gap;
}

ComplexMatrix MeanField::get_dmat_cplx(int ispin, int ikpt) const
{
    assert(ispin < this->n_spins);
    assert(ikpt < this->n_kpoints);
    auto scaled_wfc_conj = conj(wfc[ispin][ikpt]);
    for (int ib = 0; ib != this->n_bands; ib++)
        LapackConnector::scal(this->n_aos, this->wg[ispin](ikpt, ib), scaled_wfc_conj.c + n_aos * ib, 1);
    auto dmat_cplx = transpose(this->wfc[ispin][ikpt], false) * scaled_wfc_conj;
    return dmat_cplx;
}

ComplexMatrix MeanField::get_dmat_cplx_R(int ispin, const std::vector<Vector3_Order<double>>&kfrac_list, const Vector3_Order<int>& R) const
{
    ComplexMatrix dmat_cplx(this->n_aos, this->n_aos);
    for (int ik = 0; ik != this->n_kpoints; ik++)
    {
        auto ang = - (kfrac_list[ik] * R) * TWO_PI;
        complex<double> kphase = complex<double>(cos(ang), sin(ang));
        dmat_cplx += kphase * this->get_dmat_cplx(ispin, ik);
    }
    return dmat_cplx;
}

std::map<double, std::map<Vector3_Order<int>, ComplexMatrix>> MeanField::get_gf_cplx_imagtimes_Rs(
    int ispin, const std::vector<Vector3_Order<double>> &kfrac_list, std::vector<double> imagtimes,
    const std::vector<Vector3_Order<int>> &Rs) const
{
    std::map<double, std::map<Vector3_Order<int>, ComplexMatrix>> gf_tau_R;
    // NOTE: occupation must be copied here, not reference
    auto wg_empty = wg[ispin];
    // cout << "In get_gf_cplx_imagtimes_Rs ispin " << ispin << endl << wg_empty << endl;
    for (int i = 0; i != wg_empty.size; i++)
    {
        wg_empty.c[i] = 1.0 / n_kpoints - wg_empty.c[i] * (0.5 * n_spins);
        if (wg_empty.c[i] < 0) wg_empty.c[i] = 0;
        // printf("%d %f\n", i, wg_empty.c[i]);
    }
    // cout << "wg_empty " << wg_empty << endl;
    const auto wg_occ = wg[ispin] * (0.5 * n_spins);
    // cout << "wg_occ " << wg_occ << endl;
    for (const auto &tau : imagtimes)
    {
        // cout << "tau " << tau << endl;
        const auto &prefac_occ = tau > 0 ? wg_empty : wg_occ;
        // cout << "prefac_occ " << prefac_occ << endl;
        const auto scale = -tau * (eskb[ispin] - efermi);
        for (int ie = 0; ie != scale.size; ie++)
        {
            if (scale.c[ie] > 0) scale.c[ie] = 0;
            scale.c[ie] = std::exp(scale.c[ie]) * prefac_occ.c[ie];
        }
        // cout << "tau " << tau << endl << scale << endl;
        for (int ik = 0; ik != n_kpoints; ik++)
        {
            auto scaled_wfc_conj = conj(wfc[ispin][ik]);
            for (int ib = 0; ib != n_bands; ib++)
                LapackConnector::scal(n_aos, scale(ik, ib), scaled_wfc_conj.c + n_aos * ib, 1);
            const auto gf_k = transpose(wfc[ispin][ik], false) * scaled_wfc_conj;
            for (const auto &R : Rs)
            {
                double ang = -kfrac_list[ik] * R * TWO_PI;
                auto kphase = std::complex<double>(cos(ang), sin(ang));
                auto phase = kphase * (tau > 0 ? 1.0 : -1.0);
                if (gf_tau_R.count(tau) == 0 || gf_tau_R.at(tau).count(R) == 0)
                {
                    gf_tau_R[tau][R].create(n_aos, n_aos);
                }
                gf_tau_R[tau][R] += gf_k * phase;
            }
        }
        // zmy debug, print
        // for (const auto &R: Rs)
        // {
        //     cout << "tau " << tau << " R " << R << endl;
        //     print_complex_matrix("", gf_tau_R[tau][R]);
        // }
    }
    return gf_tau_R;
}

std::map<double, std::map<Vector3_Order<int>, matrix>> MeanField::get_gf_real_imagtimes_Rs(
    int ispin, const std::vector<Vector3_Order<double>> &kfrac_list, std::vector<double> imagtimes,
    const std::vector<Vector3_Order<int>> &Rs) const
{
    std::map<double, std::map<Vector3_Order<int>, matrix>> gf_tau_R;
    for (const auto &tau_gf_cplx_R: this->get_gf_cplx_imagtimes_Rs(ispin, kfrac_list, imagtimes, Rs))
    {
        const auto &tau = tau_gf_cplx_R.first;
        for (const auto &R_gf_cplx: tau_gf_cplx_R.second)
        {
            const auto &R = R_gf_cplx.first;
            gf_tau_R[tau][R] = R_gf_cplx.second.real();
            // cout << "tau " << tau << " R " << R << endl;
            // print_matrix("", gf_tau_R[tau][R]);
        }
    }
    return gf_tau_R;
}

void MeanField::allredue_wfc_isk()
{
    using LIBRPA::envs::mpi_comm_global_h;

    for(int is=0;is!=n_spins;is++)
        for(int ik=0;ik!=n_kpoints;ik++)
            {
                ComplexMatrix loc_wfc(n_bands,n_aos);
                ComplexMatrix glo_wfc(n_bands,n_aos);
                // if(mpi_comm_world_h.is_root())
                // {
                //     loc_wfc=wfc[is][ik];
                // }
                // mpi_comm_world_h.allreduce_ComplexMatrix(loc_wfc,glo_wfc);
                // if(!mpi_comm_world_h.is_root())
                // {
                //     wfc[is][ik]=glo_wfc;
                // }
                mpi_comm_global_h.allreduce_ComplexMatrix(wfc[is][ik],glo_wfc);
                wfc[is][ik]=glo_wfc;
            }
}
MeanField meanfield = MeanField();
