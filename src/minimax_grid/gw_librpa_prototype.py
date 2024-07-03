#!/usr/bin/env python
from __future__ import print_function
import sys
import math
import glob
import subprocess as sp
import numpy as np
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from typing import Tuple
from itertools import product
from io import StringIO
from scipy.optimize import least_squares

try:
    from minimax import init_minimax_grid, available_minimax_grid
except ImportError as e:
    raise ImportError("Require generated minimax module") from e


cols = 200
if sys.platform.lower() in ["linux", "darwin"]:
    _, cols = sp.check_output(['stty', 'size']).split()
    cols = int(cols)
np.set_printoptions(precision=5, linewidth=cols)

def _parser():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--ngrid", type=int, default=6, choices=available_minimax_grid,
                   help="# of minimax grid")
    return p


def phi_nu_x(nu, x):
    """1/(x+iy) + 1/(x-iy)"""
    return (2 * x) / (x**2 + nu**2)


def phi_tau_x(tau, x):
    """exp(-|tau*x|)"""
    return np.exp(-np.absolute(tau * x))


def target_t2f_cos(p, x, freq, taus):
    """target error function of cosine transfrom from time to freq"""
    s = 0
    for param, tau in zip(p, taus):
        s += param * phi_tau_x(tau, x)
    return s - phi_nu_x(freq, x)


def target_t2f_sin(p, x, freq, taus):
    """target error function of sine transfrom from time to freq"""
    s = np.zeros(np.shape(x))
    for param, tau in zip(p, taus):
        s += param * phi_tau_x(x, tau)
    return s - phi_nu_x(x, freq)


def target_f2t_cos(p, x, tau, freqs):
    """target error function of cosine transfrom from time to freq"""
    s = np.zeros(np.shape(x))
    for param, freq in zip(p, freqs):
        s += param * phi_nu_x(freq, x)
    return s - phi_tau_x(tau, x)


def target_f2t_sin(p, x, tau, freqs):
    """target error function of cosine transfrom from time to freq"""
    s = np.zeros(np.shape(x))
    for param, freq in zip(p, freqs):
        s += param * phi_nu_x(x, freq)
    return s - phi_tau_x(x, tau)


def generate_x(ngrids, erange):
    num_points_per_magnitude = 200
    num_x_node = (int(math.log10(erange)) + 1) * num_points_per_magnitude
    num_x_node = max(num_x_node, ngrids)

    x = np.power((erange)**(1.0 / (num_x_node - 1.0)),
                 np.arange(0, num_x_node))
    return x


def generate_trans_mat(omega_g, tau_g, erange, trans, kind):
    """generate the cosin/sine transformation matrix

    Args:
        omega_g, tau_g: arrays
        erange (float): energy range
        trans: "cos" or "sin"
        kind: "t2f" or "f2t"
    """
    assert trans in ["cos", "sin"]
    assert kind in ["t2f", "f2t"]

    ngrids = len(tau_g)
    x = generate_x(ngrids, erange)
    num_x_node = len(x)
    # print(x)
    p0 = np.ones(ngrids)

    dict_target = {
        ("t2f", "cos"): target_t2f_cos,
        ("t2f", "sin"): target_t2f_sin,
        ("f2t", "cos"): target_f2t_cos,
        ("f2t", "sin"): target_f2t_sin,
    }
    tfunc = dict_target[(kind, trans)]

    mat = np.zeros((ngrids, ngrids))
    if kind == "t2f":
        outer = omega_g
        inner = tau_g
    else:
        outer = tau_g
        inner = omega_g
    for i, grid in enumerate(outer):
        para_i = least_squares(tfunc, p0, args=(x, grid, inner))
        mat[i, :] = para_i.x

    return mat


def read_all_Cs_data():
    """"""
    def handle_one(fn, d_naos, d_nabfs, d_Cs):
        with open(fn, 'r') as h:
            line = h.readline()
            natoms = int(line.split()[0])
            lines = h.readlines()
            i = 0
            while i < len(lines):
                I, J, R1, R2, R3, N_I = map(int, lines[i].split())
                N_J, N_MU = map(int, lines[i + 1].split())
                I = I - 1
                size = N_I * N_J * N_MU
                mat = np.loadtxt(StringIO("".join(lines[i + 2:i + 2 + size])))
                mat = np.moveaxis(mat.reshape(N_I, N_J, N_MU), 2, 0)
                J = J - 1
                i += 2 + size
                d_Cs[(I, J, (R1, R2, R3))] = mat
                d_naos[I] = N_I
                d_naos[J] = N_J
                d_nabfs[I] = N_MU
                # print("maxval:", np.max(mat), "shape:", mat.shape)

    Cs_files = glob.glob("Cs_data_*.txt")
    Cs = {}
    abfs = {}
    aos = {}
    for fn in Cs_files:
        handle_one(fn, aos, abfs, Cs)
    abfs = [abfs[i] for i in range(len(abfs.items()))]
    aos = [aos[i] for i in range(len(aos.items()))]
    return aos, abfs, Cs


def read_band():
    """"""
    with open('band_out', 'r') as h:
        lines = h.readlines()
    nkpts, nspins, nbands, naos = map(int, lines[0:4])
    efermi = float(lines[4])
    lines = lines[5:]
    occ = np.zeros((nspins, nkpts, nbands))
    eigen = np.zeros((nspins, nkpts, nbands))
    for i in range(nkpts * nspins):
        ikpt, ispin = map(int, lines[i * (1 + nbands)].split())
        ikpt -= 1
        ispin -= 1
        occ_ks, eigen_ks = np.loadtxt(
            StringIO("".join(lines[i * (nbands + 1) + 1:(i + 1) * (nbands + 1)])),
            usecols=[1, 2], unpack=True)
        occ[ispin, ikpt, :] = occ_ks
        eigen[ispin, ikpt, :] = eigen_ks
    return eigen, occ, naos, efermi


def read_struct():
    """"""
    latt_recplat = np.loadtxt('stru_out', max_rows=6)
    latt = latt_recplat[0:3, 0:3]
    recplatt = latt_recplat[3:6, 0:3]
    kgrids = np.loadtxt('stru_out', skiprows=6, max_rows=1, dtype=int)
    assert (len(kgrids) == 3)
    kpts = np.loadtxt('stru_out', skiprows=7, max_rows=np.product(kgrids))
    kpts = kpts.reshape((np.product(kgrids), 3))
    return latt, recplatt, kgrids, kpts


def read_all_eigenvec(nspins, nkpts, nbands, naos):
    """"""
    def handle_one(fn, arr):
        with open(fn, 'r') as h:
            lines = h.readlines()
        block = 1 + nspins * nbands * naos
        nkpts_local = len(lines) // block
        for ik_local in range(nkpts_local):
            ik = int(lines[ik_local * block]) - 1
            eigenvec_k_real, eigenvec_k_imag = np.loadtxt(
                StringIO("".join(lines[1 + block * ik_local:block * (ik_local + 1)])), unpack=True)
            eigenvec_k = eigenvec_k_real + eigenvec_k_imag * 1j
            eigenvec_k = eigenvec_k.reshape(naos, nbands, nspins)
            eigenvec_k = np.swapaxes(eigenvec_k, 0, 2)
            arr[:, ik, :, :] = eigenvec_k

    arr = np.zeros((nspins, nkpts, nbands, naos), dtype='complex64')
    eigenvec_files = glob.glob('KS_eigenvector_*.txt')
    for fn in eigenvec_files:
        handle_one(fn, arr)
    return arr


def read_coulomb_abf(abfs, kpts_int, glob_pattern):
    """"""
    def handle_one(fn, d_couloumb, d_weights):
        with open(fn, 'r') as h:
            lines = h.readlines()
        nq_local = int(lines[0])
        nabfs = int(lines[1].split()[0])
        lines = lines[1:]
        i = 0
        while i < len(lines):
            nabfs, lb_row, ub_row, lb_col, ub_col = map(int, lines[i].split())
            # print(nabfs, lb_row, ub_row, lb_col, ub_col)
            blocksize = (ub_row - lb_row + 1) * (ub_col - lb_col + 1)
            ik_q, weight = lines[i + 1].split()
            ik_q = int(ik_q) - 1
            weight = float(weight)
            vq_real, vq_imag = np.loadtxt(StringIO("".join(lines[i + 2:i + 2 + blocksize])),
                                          unpack=True)
            vq = vq_real + vq_imag * 1j
            vq = vq.reshape((ub_row - lb_row + 1, ub_col - lb_col + 1))
            if ik_q not in d_couloumb:
                d_couloumb[ik_q] = np.zeros((nabfs, nabfs), dtype='complex64')
            d_couloumb[ik_q][lb_row - 1:ub_row, lb_col - 1:ub_col] = vq
            d_weights[ik_q] = weight
            i += 2 + blocksize

    coulomb_abf_files = glob.glob(glob_pattern)
    weights = {}
    coulomb = {}
    for fn in coulomb_abf_files:
        handle_one(fn, coulomb, weights)
    # divide into atom blocks
    natoms = len(abfs)
    coulomb_IJ = {}
    for iq, c in coulomb.items():
        q = tuple(kpts_int[iq])
        if coulomb_IJ.get(q) is None:
            coulomb_IJ[q] = {}
        for I in range(natoms):
            for J in range(I, natoms):
                coulomb_IJ[q][(I, J)] = np.zeros((abfs[I], abfs[J]), dtype='complex64')
                coulomb_IJ[q][(I, J)] = c[sum(abfs[:I]):sum(abfs[:I + 1]),
                                          sum(abfs[:J]):sum(abfs[:J + 1])]

    return coulomb_IJ, weights


def get_density_mat(ispin, eigen, occ, eigenvec, aos, kpts_int, kgrids, Rs):
    """"""
    nspins, nkpts, nbands, naos = eigenvec.shape
    eigen_s = eigen[ispin, :, :]
    occ_s = occ[ispin, :, :]
    eigenvec_s = eigenvec[ispin, :, :, :]
    dm_global = {}
    for ik in range(nkpts):
        scale = occ_s * (0.5 * nspins)
        dm_k = np.matmul(
            (eigenvec_s[ik, :, :] * scale[ik, :][:, np.newaxis]).T,
            eigenvec_s[ik, :, :].conj())
        # print(gf_k)
        for R in Rs:
            dm_temp = dm_k * np.exp(- 1.0j * get_kphase(kpts_int[ik], R, kgrids)) / nkpts
            if tuple(R) not in dm_global:
                dm_global[tuple(R)] = np.zeros((naos, naos))
            dm_global[tuple(R)] += dm_temp.real
    # divide into atom pairs
    dm_IJR = {}
    for I in range(len(aos)):
        for J in range(len(aos)):
            for R, dm in dm_global.items():
                dm_IJR[(I, J, R)] = dm[sum(aos[:I]):sum(aos[:I + 1]),
                                       sum(aos[:J]):sum(aos[:J + 1])]
    return dm_IJR


def get_greens_func(ispin, eigen, occ, eigenvec, aos, kpts_int, kgrids, Rs, efermi, taus):
    """"""
    nspins, nkpts, nbands, naos = eigenvec.shape
    eigen_s = eigen[ispin, :, :]
    occ_s = occ[ispin, :, :]
    eigenvec_s = eigenvec[ispin, :, :, :]
    gf_tIJR = {}
    for tau in taus:
        if tau not in gf_tIJR:
            gf_tIJR[tau] = {}
        if tau > 0:
            scale = np.exp(-(eigen_s - efermi) * tau) * (1 - occ_s * (0.5 * nspins))
        else:
            scale = - np.exp(-(eigen_s - efermi) * tau) * occ_s * (0.5 * nspins)
        gf_global = {}
        for ik in range(nkpts):
            gf_k = np.matmul(
                (eigenvec_s[ik, :, :] * scale[ik, :][:, np.newaxis]).T,
                eigenvec_s[ik, :, :].conj())
            # print(gf_k)
            for R in Rs:
                gf_temp = gf_k * np.exp(- 1.0j * get_kphase(kpts_int[ik], R, kgrids)) / nkpts
                if tuple(R) not in gf_global:
                    gf_global[tuple(R)] = np.zeros((naos, naos))
                gf_global[tuple(R)] += gf_temp.real
        # divide into atom pairs
        for I in range(len(aos)):
            for J in range(len(aos)):
                for R, gf in gf_global.items():
                    gf_tIJR[tau][(I, J, R)] = gf[sum(aos[:I]):sum(aos[:I + 1]),
                                                 sum(aos[:J]):sum(aos[:J + 1])]
    return gf_tIJR


def add_Rs(period: Tuple[int], *Rs: Tuple[int]) -> Tuple[int]:
    """"""
    R = np.mod(np.sum(np.array(Rs), axis=0), period)
    R = np.subtract(R, np.array(period) // 2)
    return tuple(R)


def get_kphase(k_int, R_int, period):
    """"""
    period = np.array(period)
    return np.sum(np.multiply(k_int, R_int) / period / period) * 2.0 * np.pi


# pep8: disable=E501
def compute_chi_t_IJR(Cs, taus, Rs, aos, abfs, period, gf_posi, gf_nega):
    """"""
    chi = {}
    nspins = len(gf_posi.items())
    nz = np.zeros
    for ispin in range(nspins):
        for tau in taus:
            if chi.get(tau) is None:
                chi[tau] = {}
            for Mu in range(len(abfs)):
                for (I, K, R1), C_mu in Cs.items():
                    if Mu != I:
                        continue
                    for R in Rs:
                        RmR1 = add_Rs(period, R, -np.array(R1))
                        for Nu in range(len(abfs)):
                            if chi[tau].get((Mu, Nu, tuple(R))) is None:
                                chi[tau][(Mu, Nu, tuple(R))] = nz((abfs[Mu], abfs[Nu]))
                            temp = nz((abfs[Nu], aos[I], aos[K]), dtype='complex64')
                            gf_ij_R_p = gf_posi[ispin][tau][(I, Nu, tuple(R))]
                            gf_ij_R_p_counterpart = nz((abfs[Nu], aos[I], aos[K]), dtype='complex64')
                            gf_ij_R_n_conj = gf_nega[ispin][-tau][(I, Nu, tuple(R))].conj()
                            gf_ij_R_n_conj_counterpart = nz((abfs[Nu], aos[I], aos[K]), dtype='complex64')
                            gf_kj_RmR1_p = gf_posi[ispin][tau][(K, Nu, RmR1)]
                            gf_kj_RmR1_p_counterpart = nz((abfs[Nu], aos[Nu], aos[I]), dtype='complex64')
                            gf_kj_RmR1_n_conj = gf_nega[ispin][-tau][(K, Nu, RmR1)].conj()
                            gf_kj_RmR1_n_conj_counterpart = nz((abfs[Nu], aos[Nu], aos[I]), dtype='complex64')
                            for (J, L, R2), C_nu in Cs.items():
                                if J != Nu:
                                    continue
                                mRpR1mR2 = add_Rs(period, -np.array(R2), -np.array(RmR1))
                                mRmR2 = add_Rs(period, -np.array(R2), -np.array(R))
                                gf_ij_R_p_counterpart += np.einsum('njl,lk->njk', C_nu, gf_nega[ispin][-tau][(L, K, mRpR1mR2)])
                                gf_ij_R_n_conj_counterpart += np.einsum('njl,lk->njk', C_nu, gf_posi[ispin][tau][(L, K, mRpR1mR2)].conj())
                                gf_kj_RmR1_n_conj_counterpart += np.einsum('njl,li->nji', C_nu, gf_posi[ispin][tau][(L, I, mRmR2)].conj())
                                gf_kj_RmR1_p_counterpart += np.einsum('njl,li->nji', C_nu, gf_nega[ispin][-tau][(L, I, mRmR2)])
                            temp += np.einsum('ij,njk->nik', gf_ij_R_p, gf_ij_R_p_counterpart)
                            temp += np.einsum('ij,njk->nik', gf_ij_R_n_conj, gf_ij_R_n_conj_counterpart)
                            temp += np.einsum('kj,nji->nik', gf_kj_RmR1_p, gf_kj_RmR1_p_counterpart)
                            temp += np.einsum('kj,nji->nik', gf_kj_RmR1_n_conj, gf_kj_RmR1_n_conj_counterpart)
                            chi[tau][(Mu, Nu, tuple(R))] += 2.0 * np.einsum('mik,nik->mn', C_mu, temp.real) / nspins
    return chi


def power_hemat(mat, power, thres):
    """"""
    ev, evec = np.linalg.eigh(mat, 'U')
    ev[ev < thres] = 0.0
    ev = np.power(ev, power)
    return np.matmul(evec, np.multiply(ev[:, np.newaxis], evec.conj().T))


def compute_Wc_fq(abfs, chi_fq_IJ, coulomb, coulomb_cut):
    """"""
    Wc_fq = {}
    nz = np.zeros
    for freq, q_chi_IJ in chi_fq_IJ.items():
        Wc_fq[freq] = {}
        for q, chi_IJ in q_chi_IJ.items():
            Wc_fq[freq][q] = {}
            coul_q_IJ = coulomb[q]
            coul_cut_q_IJ = coulomb_cut[q]
            chi_all = nz((sum(abfs), sum(abfs)), dtype='complex64')
            coul_q_all = nz((sum(abfs), sum(abfs)), dtype='complex64')
            for (I, J), chi in chi_IJ.items():
                chi_all[sum(abfs[:I]):sum(abfs[:I + 1]),
                        sum(abfs[:J]):sum(abfs[:J + 1])] = chi[:, :]
            for (I, J), c in coul_q_IJ.items():
                coul_q_all[sum(abfs[:I]):sum(abfs[:I + 1]),
                           sum(abfs[:J]):sum(abfs[:J + 1])] = c[:, :]
                # if I != J:
                #     coul_q_all[sum(abfs[:J]):sum(abfs[:J + 1]),
                #                sum(abfs[:I]):sum(abfs[:I + 1])] = c[:, :].conj().T
            coul_q_all = power_hemat(coul_q_all, 0.5, 1e-4)
            eps = np.identity(sum(abfs)) - np.matmul(coul_q_all, np.matmul(chi_all, coul_q_all))
            coul_q_all[:, :] = 0.0
            for (I, J), c in coul_cut_q_IJ.items():
                coul_q_all[sum(abfs[:I]):sum(abfs[:I + 1]), sum(abfs[:J]):sum(abfs[:J + 1])] = c[:, :]
                # if I != J:
                #     coul_q_all[sum(abfs[:J]):sum(abfs[:J + 1]), sum(abfs[:I]):sum(abfs[:I + 1])] = c[:, :].conj().T
            coul_q_all = power_hemat(coul_q_all, 0.5, 1e-4)
            Wc = np.matmul(coul_q_all, np.matmul(np.linalg.inv(eps) - np.identity(sum(abfs)), coul_q_all))
            for I in range(len(abfs)):
                for J in range(len(abfs)):
                    Wc_fq[freq][q][(I, J)] = nz((abfs[I], abfs[J]), dtype='complex64')
                    Wc_fq[freq][q][(I, J)][:, :] = Wc[sum(abfs[:I]):sum(abfs[:I + 1]), sum(abfs[:J]):sum(abfs[:J + 1])]
    return Wc_fq


def contraction_ACBC(aos, Rs, period, Cs, A_IJR, B_MuNuR):
    """"""
    # use equation 29-32 of selfenergy.pdf
    contracted = {}
    for R in Rs:
        for Rp in Rs:
            for (I, K, R1mR), C_Mu in Cs.items():
                for (J, L, R2mRp), C_Nu in Cs.items():
                    temp = np.zeros((aos[I], aos[K], aos[J], aos[L]))
                    for (Mu, Nu, RmRp), B in B_MuNuR.items():
                        if (I, J) != (Mu, Nu) or RmRp != add_Rs(period, R, -np.array(Rp)):
                            continue
                        temp += np.einsum('mik,mn,njl->ikjl', C_Mu, B, C_Nu)
                    R2mR1 = add_Rs(period, R2mRp, Rp, R1mR, R)
                    R2mR = add_Rs(period, R2mR1, R1mR)
                    RpmR1 = add_Rs(period, Rp, R, R1mR)
                    RpmR = add_Rs(period, RpmR1, R1mR)
                    if (I, J, RpmR) not in contracted:
                        contracted[(I, J, RpmR)] = np.zeros((aos[I], aos[J]))
                    contracted[(I, J, RpmR)] += np.einsum('ikjl,kl->ij', temp, A_IJR[(K, L, R2mR1)])
                    if (K, J, RpmR1) not in contracted:
                        contracted[(K, J, RpmR1)] = np.zeros((aos[K], aos[J]))
                    contracted[(K, J, RpmR1)] += np.einsum('ikjl,il->kj', temp, A_IJR[(I, L, R2mR)])
                    if (I, L, R2mR) not in contracted:
                        contracted[(I, L, R2mR)] = np.zeros((aos[I], aos[L]))
                    contracted[(I, L, R2mR)] += np.einsum('ikjl,kj->il', temp, A_IJR[(K, J, RpmR1)])
                    if (K, L, R2mR1) not in contracted:
                        contracted[(K, L, R2mR1)] = np.zeros((aos[K], aos[L]))
                    contracted[(K, L, R2mR1)] += np.einsum('ikjl,ij->kl', temp, A_IJR[(I, J, RpmR)])
    return contracted


def rotate_IJ_to_KS_basis(A_IJ, aos, eigenvec, ispin, ik):
    """"""
    nspins, nkpts, nbands, naos = eigenvec.shape
    A_full = np.zeros((sum(aos), sum(aos)), dtype='complex64')
    for I in range(len(aos)):
        for J in range(len(aos)):
            A_full[sum(aos[:I]):sum(aos[:I + 1]),
                   sum(aos[:J]):sum(aos[:J + 1])] = A_IJ[(I, J)][:, :]
    eigenvec_sk = eigenvec[ispin, ik, :, :]
    return np.matmul(eigenvec_sk.conj(), np.matmul(A_full, eigenvec_sk.T))


if __name__ == "__main__":
    args = _parser().parse_args()
    aos, abfs, Cs = read_all_Cs_data()
    nabfs = sum(abfs)
    # print(aos, abfs)
    eigen, occ, naos, efermi = read_band()
    nspins, nkpts, nbands = eigen.shape
    # print(nspins, nkpts, nbands, naos)
    latt, recplatt, kgrids, kpts = read_struct()
    # print(kgrids)
    Rs = np.array(list(product(*list(map(range, kgrids)))), dtype='int64') - np.array(kgrids) // 2
    kpts_frac = np.matmul(kpts, latt.T) / 2.0 / np.pi
    kpts_frac[np.abs(kpts_frac) < 1e-10] = 0.0
    kpts_int = np.asarray(np.rint(kpts_frac * kgrids), dtype='int64')
    # print(kpts_frac)
    # print(kpts_int)
    eigenvec = read_all_eigenvec(nspins, nkpts, nbands, naos)
    # print(eigenvec)
    coul_all, coul_weights = read_coulomb_abf(abfs, kpts_int, 'coulomb_mat_*.txt')
    coul_cut, coul_weights = read_coulomb_abf(abfs, kpts_int, 'coulomb_cut_*.txt')
    # print(coul_all)
    # print(coul_cut)
    # print(coul_weights)

    emin = np.min(eigen[occ < 0.001 * 2 / nspins]) - np.max(eigen[occ > 0.999 * 2 / nspins])
    emax = np.max(eigen) - np.min(eigen)
    erange = emax / emin
    print("emin/emax/erange:", emin, emax, erange)
    ngrids = args.ngrid
    freq_grids, freq_weights, time_grids, time_weights = init_minimax_grid(ngrids, emin, emax)
    cos_t2f = generate_trans_mat(freq_grids / emin, time_grids * emin, erange, "cos", "t2f") / emin
    sin_t2f = generate_trans_mat(freq_grids / emin, time_grids * emin, erange, "sin", "t2f") / emin
    cos_f2t = generate_trans_mat(freq_grids / emin, time_grids * emin, erange, "cos", "f2t") * emin
    # cos_f2t = np.linalg.inv(cos_t2f)  # force identity
    sin_f2t = generate_trans_mat(freq_grids / emin, time_grids * emin, erange, "sin", "f2t") * emin
    # sin_f2t = np.linalg.inv(sin_t2f)  # force identity
    # print(cos_f2t)
    # print(np.matmul(cos_t2f, cos_f2t))
    # print(np.matmul(sin_t2f, sin_f2t))

    coul_cut_IJR = {}
    for R in Rs:
        for q, IJcoul in coul_cut.items():
            for (I, J), c in IJcoul.items():
                if (I, J, tuple(R)) not in coul_cut_IJR.keys():
                    coul_cut_IJR[(I, J, tuple(R))] = np.zeros((abfs[I], abfs[J]), dtype='complex64')
                scale = np.exp(- 1.0j * get_kphase(q, R, kgrids)) / nkpts
                coul_cut_IJR[(I, J, tuple(R))] += c * scale
    for IJR in list(coul_cut_IJR.keys()):
        coul_cut_IJR[IJR] = coul_cut_IJR[IJR].real

    gf_posi = {}
    gf_nega = {}
    for ispin in range(nspins):
        gf_posi[ispin] = get_greens_func(ispin, eigen, occ, eigenvec, aos, kpts_int,
                                         kgrids, Rs, efermi, time_grids)
        gf_nega[ispin] = get_greens_func(ispin, eigen, occ, eigenvec, aos, kpts_int,
                                         kgrids, Rs, efermi, -time_grids)
    # for t, Rgf in gf_posi.items():
    #     for IJR, gf in Rgf.items():
    #         print(t, IJR)
    #         print(gf)
    #         print(-t, IJR)
    #         print(gf_nega[-t][IJR])

    print("Start compute_chi_t_IJR")
    chi_t_IJR = compute_chi_t_IJR(Cs, time_grids, Rs, aos, abfs, kgrids, gf_posi, gf_nega)
    print("Done compute_chi_t_IJR")
    print("time-domain real-space chi")
    for t, IJRchi in chi_t_IJR.items():
        for R, chi in IJRchi.items():
            print(t, R)
            print(chi)

    chi_wq_IJ = {}
    for iw, w in enumerate(freq_grids):
        chi_wq_IJ[w] = {}
        for iq in coul_weights.keys():
            q_int = tuple(kpts_int[iq])
            chi_wq_IJ[w][q_int] = {}
            for tau, IJRchi in chi_t_IJR.items():
                it = tuple(time_grids).index(tau)
                for (I, J, R), chi in IJRchi.items():
                    if chi_wq_IJ[w][q_int].get((I, J)) is None:
                        chi_wq_IJ[w][q_int][(I, J)] = np.zeros((abfs[I], abfs[J]), dtype='complex64')
                    scale = np.exp(1.0j * get_kphase(q_int, R, kgrids)) * cos_t2f[iw, it]
                    print(iw, it, chi[0, 0], scale, chi[0, 0] * scale)
                    chi_wq_IJ[w][q_int][(I, J)] += chi * scale

    print("", "==========", sep="\n")
    print("frequency-domain reciprocal-space chi")
    # for w, q_IJchi in chi_wq_IJ.items():
    #     for q, IJchi in q_IJchi.items():
    #         for IJ, chi in IJchi.items():
    #             print(w, q, IJ)
    #             print(chi)

    print("", "==========", sep="\n")
    print("compute frequency-domain reciprocal-space screened interaction")
    Wc_fq_IJ = compute_Wc_fq(abfs, chi_wq_IJ, coul_all, coul_cut)
    # for w, q_IJWc in Wc_fq_IJ.items():
    #     for q, IJWc in q_IJWc.items():
    #         for IJ, Wc in IJWc.items():
    #             print(w, q, IJ)
    #             print(Wc)

    print("", "==========", sep="\n")
    print("transfer to time-domain real-space screened interaction")
    Wc_t_IJR = {}
    for it, t in enumerate(time_grids):
        Wc_t_IJR[t] = {}
        for R in Rs:
            for w, q_IJWc in Wc_fq_IJ.items():
                iw = tuple(freq_grids).index(w)
                for q, IJWc in q_IJWc.items():
                    for (I, J), Wc in IJWc.items():
                        if (I, J, tuple(R)) not in Wc_t_IJR[t].keys():
                            Wc_t_IJR[t][(I, J, tuple(R))] = np.zeros((abfs[I], abfs[J]), dtype='complex64')
                        scale = np.exp(- 1.0j * get_kphase(q, R, kgrids)) * cos_f2t[it, iw] / nkpts
                        Wc_t_IJR[t][(I, J, tuple(R))] += Wc * scale
        for IJR in list(Wc_t_IJR[t].keys()):
            Wc_t_IJR[t][IJR] = Wc_t_IJR[t][IJR].real

    for t, IJRWc in Wc_t_IJR.items():
        for IJR, Wc in IJRWc.items():
            print(t, IJR)
            print(Wc)

    print("", "==========", sep="\n")
    print("Perform contration")
    sigc_posi_is_t_IJR = {}
    sigc_nega_is_t_IJR = {}
    for ispin in range(nspins):
        sigc_posi_is_t_IJR[ispin] = {}
        sigc_nega_is_t_IJR[ispin] = {}
        for t in time_grids:
            sigc_posi_is_t_IJR[ispin][t] = contraction_ACBC(aos, Rs, kgrids, Cs, gf_posi[ispin][t], Wc_t_IJR[t])
            sigc_nega_is_t_IJR[ispin][-t] = contraction_ACBC(aos, Rs, kgrids, Cs, gf_nega[ispin][-t], Wc_t_IJR[t])

    for ispin in range(nspins):
        print("Spin", ispin + 1)
        for t in time_grids:
            print("Time", t)
            print("Positive Sigc c")
            for IJR, sigc in sigc_posi_is_t_IJR[ispin][t].items():
                print(IJR)
                print(sigc)
            print("Negative Sigc c")
            for IJR, sigc in sigc_nega_is_t_IJR[ispin][-t].items():
                print(IJR)
                print(sigc)

    print("", "==========", sep="\n")
    print("Transform Sigc to frequency-domain reciprocal-space")

    sigc_ispin_fq_IJ = {}
    for ispin in range(nspins):
        sigc_ispin_fq_IJ[ispin] = {}
        for iw, w in enumerate(freq_grids):
            sigc_ispin_fq_IJ[ispin][w] = {}
            for iq in coul_weights.keys():
                q = tuple(kpts_int[iq])
                sigc_ispin_fq_IJ[ispin][w][q] = {}
                for it, t in enumerate(time_grids):
                    for I, J, R in sigc_posi_is_t_IJR[ispin][t].keys():
                        if sigc_ispin_fq_IJ[ispin][w][q].get((I, J)) is None:
                            sigc_ispin_fq_IJ[ispin][w][q][(I, J)] = np.zeros((aos[I], aos[J]), dtype='complex64')
                        # print(sigc_posi_is_t_IJR[ispin][t][(I, J, R)][0, 0])
                        # print(sigc_nega_is_t_IJR[ispin][-t][(I, J, R)][0, 0])
                        sigc_cos = sigc_posi_is_t_IJR[ispin][t][(I, J, R)] + sigc_nega_is_t_IJR[ispin][-t][(I, J, R)]
                        sigc_sin = sigc_posi_is_t_IJR[ispin][t][(I, J, R)] - sigc_nega_is_t_IJR[ispin][-t][(I, J, R)]
                        # NOTE: half here to get almost same result as FHI-aims, but I don't know where this comes from
                        temp = 0.5 * (sigc_cos * cos_t2f[iw, it] + sigc_sin * sin_t2f[iw, it] * 1.0j)
                        # temp = sigc_cos * sin_t2f[iw, it] + sigc_sin * cos_t2f[iw, it] * 1.0j
                        sigc_ispin_fq_IJ[ispin][w][q][(I, J)] += temp * np.exp(1.0j * get_kphase(q, R, kgrids))

    print("", "==========", sep="\n")
    print("Rotate to KS band basis")
    for ispin in range(nspins):
        for iw, w in enumerate(freq_grids):
            for iq in coul_weights.keys():
                q = tuple(kpts_int[iq])
                sigc_IJ = sigc_ispin_fq_IJ[ispin][w][q]
                print("Sigc for ispin", ispin + 1, "freq", w, "q", q)
                # print("I J block")
                # for (I, J), sigc in sigc_IJ.items():
                #     print(I, J)
                #     print(sigc)
                sigc_KS = rotate_IJ_to_KS_basis(sigc_IJ, aos, eigenvec, ispin, iq)
                print(sigc_KS)
                print("Diagonal terms")
                print(sigc_KS.diagonal())
