#include "meanfield.h"

#include <iostream>
#include <limits>

#include "../utils/constants.h"
#include "../utils/error.h"
#include "../math/lapack_connector.h"
#include "../math/utils_matrix_mpi.h"
#include "../mpi/global_mpi.h"

namespace librpa_int {

void MeanField::resize(int ns, int nk, int nb, int nao, int st_ib, int nb_local, int st_iao, int nao_local)
{
    if (ns == 0 || nk == 0 || nb == 0 || nao == 0)
        throw LIBRPA_RUNTIME_ERROR("encounter zero dimension");
    if (n_spins != 0)
    {
        eskb.clear();
        wg.clear();
        wfc.clear();
    }

    n_spins = ns;
    n_kpoints = nk;
    n_states = nb;
    n_aos = nao;
    i_ao_start = st_iao;
    n_aos_local = nao_local;
    i_state_start = st_ib;
    n_states_local = nb_local;

    eskb.resize(n_spins);
    wg.resize(n_spins);
    wfc.clear();

    for (int is = 0; is < n_spins; is++)
    {
        eskb[is].create(n_kpoints, n_states);
        wg[is].create(n_kpoints, n_states);
    }
}

std::vector<int> MeanField::get_iks_local() const
{
    std::vector<int> iks_local;
    auto it = this->wfc.find(0);  // assume same k-points in all spin channels
    if (it != this->wfc.cend())
    {
        for (const auto &[ik, _]: it->second)
        {
            iks_local.push_back(ik);
        }
    }
    return iks_local;
}

MeanField::MeanField(int ns, int nk, int nb, int nao)
{
    resize(ns, nk, nb, nao, 0, nb, 0, nao);
}

MeanField::MeanField(int ns, int nk, int nb, int nao, int st_ib, int nb_local, int st_iao, int nao_local)
{
    resize(ns, nk, nb, nao, st_ib, nb_local, st_iao, nao_local);
}

void MeanField::set(int ns, int nk, int nb, int nao)
{
    if (n_spins != 0 || n_kpoints != 0 || n_states != 0 || n_aos != 0)
    {
        std::cout << n_spins << n_kpoints << n_states << n_aos << std::endl;
        throw LIBRPA_RUNTIME_ERROR("MeanField object already set");
    }
    resize(ns, nk, nb, nao, 0, nb, 0, nao);
}

void MeanField::set(int ns, int nk, int nb, int nao, int st_ib, int nb_local, int st_iao, int nao_local)
{
    if (n_spins != 0 || n_kpoints != 0 || n_states != 0 || n_aos != 0)
    {
        std::cout << n_spins << n_kpoints << n_states << n_aos << std::endl;
        throw LIBRPA_RUNTIME_ERROR("MeanField object already set");
    }
    resize(ns, nk, nb, nao, st_ib, nb_local, st_iao, nao_local);
}

MeanField::MeanField(const MeanField &mf)
{
    resize(mf.n_spins, mf.n_kpoints, mf.n_states, mf.n_aos,
           mf.i_state_start, mf.n_states_local, mf.i_ao_start, mf.n_aos_local);
    // FIXME: copy data, not tested
    eskb = mf.eskb;
    wg = mf.wg;
    wfc = mf.wfc;
    efermi = mf.efermi;
}

double MeanField::get_E_min_max(double &emin, double &emax) const
{
    // double lb = eskb[0](0, 0);
    // double ub = eskb[0](0, n_states - 1);
    double lb = std::numeric_limits<double>::max();
    double ub = std::numeric_limits<double>::min();
    for (int is = 0; is != n_spins; is++)
        for (int ik = 0; ik != n_kpoints; ik++)
        {
            lb = (lb > eskb[is](ik, 0)) ? eskb[is](ik, 0) : lb;
            ub = (ub < eskb[is](ik, n_states-1)) ? eskb[is](ik, n_states-1) : ub;
        }
    double gap = get_band_gap();
    emax = ub - lb;
    emin = gap;
    return emin;
}

double MeanField::get_band_gap() const
{
    double homo = -1e6, lumo = 1e6;
    double gap = lumo - homo;
    for (int is = 0; is != n_spins; is++)
    {
        // FIXME: should be nspins/nkpoints?
        const double midpoint = 1.0 / (n_spins * n_kpoints);
        //print_matrix("mf.eskb: ",this->eskb[is]);
        for (int ik = 0; ik != n_kpoints; ik++)
        {
            int homo_level = -1;
            for (int n = 0; n != n_states; n++)
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

// get_dmat_cplx can be used for both serial and k-porallel version
ComplexMatrix MeanField::get_dmat_cplx(int ispin, int ikpt) const
{
    assert(ispin < this->n_spins);
    assert(ikpt < this->n_kpoints);
    auto scaled_wfc_conj = conj(wfc.at(ispin).at(ikpt));
    const double occ_thres = 1e-4 / n_kpoints;
    int nocc;
    for (nocc = 0; nocc != this->n_states; nocc++)
    {
        // Renormalize to single spin channel. Need to adpat SOC case (n_spins = 1 but remove 0.5)
        const auto weight = this->wg[ispin](ikpt, nocc) * 0.5 * n_spins;
        if (weight < occ_thres) break;
        LapackConnector::scal(this->n_aos, weight, scaled_wfc_conj.c + n_aos * nocc, 1);
    }
    ComplexMatrix dmat_cplx(n_aos, n_aos);
    const auto &wfc_sk = this->wfc.at(ispin).at(ikpt);
    LapackConnector::gemm('T', 'N', n_aos, n_aos, nocc, 1.0,
                          wfc_sk.c, n_aos, scaled_wfc_conj.c, n_aos,
                          0.0, dmat_cplx.c, n_aos);
    // auto dmat_cplx = transpose(wfc_sk, false) * scaled_wfc_conj;
    return dmat_cplx;
}

ComplexMatrix MeanField::get_dmat_cplx_R(int ispin, const std::vector<Vector3_Order<double>>&kfrac_list, const Vector3_Order<int>& R) const
{
    ComplexMatrix dmat_cplx(this->n_aos, this->n_aos);
    for (int ik = 0; ik != this->n_kpoints; ik++)
    {
        auto ang = - (kfrac_list[ik] * R) * TWO_PI;
        complex<double> kphase = complex<double>(cos(ang), sin(ang));
        global::ofs_myid << "R " << R << " ik " << ik << " kfrac " << kfrac_list[ik] << " phase " << kphase << std::endl;
        dmat_cplx += kphase * this->get_dmat_cplx(ispin, ik);
    }
    return dmat_cplx;
}

std::map<Vector3_Order<int>, ComplexMatrix> MeanField::get_dmat_cplx_Rs(
    int ispin, const std::vector<Vector3_Order<double>> &kfrac_list,
    const std::vector<Vector3_Order<int>> &Rs) const
{
    std::map<Vector3_Order<int>, ComplexMatrix> dmat_cplx_all;
    if (Rs.size() == 0) return dmat_cplx_all;
    for (const auto &R: Rs)
    {
        dmat_cplx_all[R] = ComplexMatrix(this->n_aos, this->n_aos);
    }
    for (int ik = 0; ik != this->n_kpoints; ik++)
    {
        const auto &kmat = this->get_dmat_cplx(ispin, ik);
        for (auto &[R, rmat]: dmat_cplx_all)
        {
            auto ang = - (kfrac_list[ik] * R) * TWO_PI;
            complex<double> kphase = complex<double>(cos(ang), sin(ang));
            rmat += kphase * rmat;
        }
    }
    return dmat_cplx_all;
}

ComplexMatrix MeanField::get_gf_cplx_imagtime(int ispin, int ikpt, double tau) const
{
    assert(ispin < this->n_spins);
    assert(ikpt < this->n_kpoints);

    const double scale_spin = 0.5 * n_spins;  // TODO: adapt SOC

    std::vector<double> wg_sk(wg[ispin].c + n_states * ikpt, wg[ispin].c + n_states * (ikpt + 1));
    std::vector<double> wg_empty_sk(wg_sk);
    for (int ib = 0; ib < n_states; ib++)
    {
        wg_sk[ib] *= scale_spin;
        wg_empty_sk[ib] = 1.0 / n_kpoints - wg_empty_sk[ib] * scale_spin;
        if (wg_empty_sk[ib] < 0.0) wg_empty_sk[ib] = 0.0;
    }
    const auto &prefac_occ = tau > 0 ? wg_empty_sk : wg_sk;

    std::vector<double> scale(eskb[ispin].c + n_states * ikpt, eskb[ispin].c + n_states * (ikpt + 1));
    for (int ib = 0; ib < n_states; ib++)
    {
        scale[ib] = -tau * (scale[ib] - efermi);
        if (scale[ib] > 0) scale[ib] = 0.0;
        scale[ib] = std::exp(scale[ib]) * prefac_occ[ib];
        if (tau <= 0) scale[ib] *= -1.0;
    }
    auto scaled_wfc_conj = conj(wfc.at(ispin).at(ikpt));
    for (int ib = 0; ib < n_states; ib++)
        LapackConnector::scal(n_aos, scale[ib], scaled_wfc_conj.c + n_aos * ib, 1);
    return transpose(wfc.at(ispin).at(ikpt), false) * scaled_wfc_conj;
}

std::map<double, std::map<Vector3_Order<int>, ComplexMatrix>> MeanField::get_gf_cplx_imagtimes_Rs(
    int ispin, const std::vector<Vector3_Order<double>> &kfrac_list, std::vector<double> imagtimes,
    const std::vector<Vector3_Order<int>> &Rs) const
{
    std::map<double, std::map<Vector3_Order<int>, ComplexMatrix>> gf_tau_R;
    const double scale_spin = 0.5 * n_spins;  // TODO: adapt SOC
    // NOTE: occupation must be copied here, not reference
    auto wg_empty = wg[ispin];
    // cout << "In get_gf_cplx_imagtimes_Rs ispin " << ispin << endl << wg_empty << endl;
    for (size_t i = 0; i != wg_empty.size; i++)
    {
        wg_empty.c[i] = 1.0 / n_kpoints - wg_empty.c[i] * scale_spin;
        if (wg_empty.c[i] < 0) wg_empty.c[i] = 0;
        // printf("%d %f\n", i, wg_empty.c[i]);
    }
    // cout << "wg_empty " << wg_empty << endl;
    const auto wg_occ = wg[ispin] * scale_spin;
    // cout << "wg_occ " << wg_occ << endl;
    for (const auto &tau : imagtimes)
    {
        gf_tau_R[tau] = {};
        // cout << "tau " << tau << endl;
        // Empty local R, cycle after initialize the tau container
        if (Rs.size() == 0) continue;
        const auto &prefac_occ = tau > 0 ? wg_empty : wg_occ;
        // cout << "prefac_occ " << prefac_occ << endl;
        const auto scale = -tau * (eskb[ispin] - efermi);
        for (size_t ie = 0; ie != scale.size; ie++)
        {
            if (scale.c[ie] > 0) scale.c[ie] = 0;
            scale.c[ie] = std::exp(scale.c[ie]) * prefac_occ.c[ie];
        }
        // ofs_myid cout << "tau " << tau << endl << scale << endl;
        for (int ik = 0; ik != n_kpoints; ik++)
        {
            auto scaled_wfc_conj = conj(wfc.at(ispin).at(ik));
            for (int ib = 0; ib != n_states; ib++)
                LapackConnector::scal(n_aos, scale(ik, ib), scaled_wfc_conj.c + n_aos * ib, 1);
            const auto gf_k = transpose(wfc.at(ispin).at(ik), false) * scaled_wfc_conj;
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
        gf_tau_R[tau] = {};
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

// void MeanField::allredue_wfc_isk()
// {
//     using librpa_int::global::mpi_comm_global_h;
// 
//     for(int is=0;is!=n_spins;is++)
//         for(int ik=0;ik!=n_kpoints;ik++)
//             {
//                 ComplexMatrix loc_wfc(n_states,n_aos);
//                 ComplexMatrix glo_wfc(n_states,n_aos);
//                 // if(mpi_comm_world_h.is_root())
//                 // {
//                 //     loc_wfc=wfc[is][ik];
//                 // }
//                 // mpi_comm_world_h.allreduce_ComplexMatrix(loc_wfc,glo_wfc);
//                 // if(!mpi_comm_world_h.is_root())
//                 // {
//                 //     wfc[is][ik]=glo_wfc;
//                 // }
//                 librpa_int::allreduce_ComplexMatrix(wfc[is][ik],glo_wfc,mpi_comm_global_h.comm);
//                 wfc[is][ik]=glo_wfc;
//             }
// }

}
