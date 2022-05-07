#include <iostream>
#include "chi0.h"
#include "matrix.h"
#include "complexmatrix.h"
#include "lapack_connector.h"
#include "input.h"
#include "constants.h"

void Chi0::compute(const apair_R_mat_t &LRI_Cs,
                   const vector<Vector3_Order<int>> &Rlist,
                   const Vector3_Order<int> &R_period,
                   TFGrids::GRID_TYPES gt, bool use_space_time)
{
    gf_save = gf_discard = 0; 
    // prepare the time-freq grids
    switch (gt)
    {
        case ( TFGrids::GRID_TYPES::Minimax):
        {
            double emin, emax;
            mf.get_E_min_max(emin, emax);
            tfg.generate_minimax(emin, emax);
            break;
        }
        case ( TFGrids::GRID_TYPES::EvenSpaced_TF ):
        {
            // WARN: only used to quick test, adjust emin and interval manually
            tfg.generate_evenspaced_tf(0.0304601, 0.0, 0.0116073, 0.0);
            break;
        }
        default:
            throw invalid_argument("chi0 on requested grid is not implemented");
    }

    if (use_space_time && !tfg.has_time_grids())
        throw invalid_argument("current grid type does not have time grids");

    if (use_space_time)
    {
        // use space-time method
        for (int iR = 0; iR < Rlist.size(); iR++)
            for (int itau = 0; itau < tfg.size(); itau++)
            {
                build_gf_Rt(iR, Rlist[iR], itau, 'O');
                build_gf_Rt(iR, Rlist[iR], itau, 'V');
            }
        // For debug use
        /* printf("spin 0, first tau GF: %f\n", tfg.get_time_nodes()[0]); */
        /* for (int iR = 0; iR != Rlist.size(); iR++) */
        /* { */
        /*     for (int I = 0; I != natom; I++) */
        /*         for (int J = 0; J != natom; J++) */
        /*         { */
        /*             cout << "Green_atom  IJ R:  " << I << "  " << J << "    " << Rlist[iR] << endl; */
        /*             print_matrix("unocc Green-atom", gf_unocc_Rt[0][I][J][iR][0]); */
        /*             print_matrix("occ Green-atom", gf_occ_Rt[0][I][J][iR][0]); */
        /*         } */
        /* } */
        /* printf("spin 0 GF: I %d J %d iR %d tau %f R %d %d %d\n", 0, 0, 0, tfg.get_time_nodes()[0], Rlist[0].x, Rlist[0].y, Rlist[0].z); */
        /* print_matrix("unocc Green-atom", gf_unocc_Rt[0][0][0][0][0]); */
        /* print_matrix("occ Green-atom", gf_occ_Rt[0][0][0][0][0]); */
        cout << " FINISH GREEN FUN" << endl;
        cout << " Green threshold: " << gf_R_threshold << endl;
        cout << " Green save num: " << gf_save << endl;
        cout << " Green discard num: " << gf_discard << endl;
    }
    else
    {
        // conventional method does not need to build Green's function explicitly
    }
    return;
}

void Chi0::build_gf_Rt(size_t iR, Vector3_Order<int> R, size_t itau, char ov)
{
    double tau;
    auto nkpts = mf.get_n_kpoints();
    auto nspins = mf.get_n_spins();
    auto nbands = mf.get_n_bands();
    auto naos = mf.get_n_aos();

    if ( ov == 'O') // occupied Green's function for negative imaginary time
        tau = - tfg.get_time_nodes()[itau];
    else if ( ov == 'V') // unoccupied Green's function for positive imaginary time
        tau = tfg.get_time_nodes()[itau];
    else
        throw invalid_argument("Invalid occupancy type");
    // alias the real-space Green's function to compute
    auto & gf_Rt = ov == 'O' ? gf_occ_Rt : gf_unocc_Rt;

    // temporary Green's function
    matrix gf_Rt_is_global(naos, naos);

    for (int is = 0; is != nspins; is++)
    {
        gf_Rt_is_global.zero_out();
        auto wg = meanfield.get_weight()[is];
        if ( tau > 0)
            for (int i = 0; i != wg.size; i++)
            {
                wg.c[i] = 1.0 / nkpts * nspins - wg.c[i];
                if (wg.c[i] < 0) wg.c[i] = 0;
            }
        matrix scale(nkpts, nbands);
        // tau-energy phase
        scale = - 0.5 * tau * (mf.get_eigenvals()[is] - mf.get_efermi());
        /* print_matrix("-(e-ef)*tau", scale); */
        for (int ie = 0; ie != scale.size; ie++)
        {
            // NOTE: enforce non-positive phase
            if ( scale.c[ie] > 0) scale.c[ie] = 0;
            scale.c[ie] = std::exp(scale.c[ie]) * wg.c[ie];
        }
        /* print_matrix("exp(-dE*tau)", scale); */
        for (int ik = 0; ik != nkpts; ik++)
        {
            double ang = - klist[ik] * (R * latvec) * TWO_PI;
            complex<double> kphase = complex<double>(cos(ang), sin(ang));
            /* printf("kphase %f %fj\n", kphase.real(), kphase.imag()); */
            auto scaled_wfc_conj = conj(mf.get_eigenvectors()[is][ik]);
            for ( int ib = 0; ib != nbands; ib++)
                LapackConnector::scal(naos, scale(ik, ib), scaled_wfc_conj.c + naos * ib, 1);
            gf_Rt_is_global += (kphase * transpose(mf.get_eigenvectors()[is][ik], false) * scaled_wfc_conj).real();
        }
        if ( tau < 0 ) gf_Rt_is_global *= -1.;
        for ( int I = 0; I != natom; I++)
        {
            const auto I_num = atom_nw[I];
            for (int J = 0; J != natom; J++)
            {
                const size_t J_num = atom_nw[J];
                // Green_atom[is][I][J][time_tau][R].create(I_num,J_num);
                matrix tmp_green(I_num, J_num);
                for (size_t i = 0; i != I_num; i++)
                {
                    size_t i_glo = atom_iw_loc2glo(I, i);
                    for (size_t j = 0; j != J_num; j++)
                    {
                        size_t j_glo = atom_iw_loc2glo(J, j);
                        tmp_green(i, j) = gf_Rt_is_global(i_glo, j_glo);
                    }
                }
                if (tmp_green.absmax() > gf_R_threshold)
                {
                    // cout<<" max_green_ele:  "<<tmp_green.absmax()<<endl;
                    gf_Rt[is][I][J][iR][itau] = std::move(tmp_green);
                    gf_save++;
                }
                else
                {
                    gf_discard++;
                }
            }
        }
    }
}
