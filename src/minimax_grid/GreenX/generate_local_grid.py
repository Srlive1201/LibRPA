"""write the minimax grids and cosine/sine transformation matrix"""
from __future__ import print_function
import math
import sys
import os
import numpy as np
from scipy.optimize import least_squares
# from scipy.optimize import minimize,Bounds,least_squares
# import matplotlib.pyplot as plt

fnerror = IOError
if sys.version_info.major >= 3:
    fnerror = FileNotFoundError


def phi_nu_x(nu, x):
    """1/(x+iy) + 1/(x-iy)"""
    return (2*x)/(x**2+nu**2)


def phi_tau_x(tau, x):
    """exp(-|tau*x|)"""
    return np.exp(-np.absolute(tau*x))


def target_t2f_cos(p, x, freq, taus):
    """target error function of cosine transfrom from time to freq"""
    s = 0
    for param, tau in zip(p, taus):
        s += param*phi_tau_x(tau, x)
    return s - phi_nu_x(freq, x)


def target_t2f_sin(p, x, freq, taus):
    """target error function of sine transfrom from time to freq"""
    s = np.zeros(np.shape(x))
    for param, tau in zip(p, taus):
        s += param*phi_tau_x(x, tau)
    return s - phi_nu_x(x, freq)


def target_f2t_cos(p, x, tau, freqs):
    """target error function of cosine transfrom from time to freq"""
    s = np.zeros(np.shape(x))
    for param, freq in zip(p, freqs):
        s += param*phi_nu_x(freq, x)
    return s - phi_tau_x(tau, x)


def target_f2t_sin(p, x, tau, freqs):
    """target error function of cosine transfrom from time to freq"""
    s = np.zeros(np.shape(x))
    for param, freq in zip(p, freqs):
        s += param*phi_nu_x(x, freq)
    return s - phi_tau_x(x, tau)


def read_grid(ngrids, erange):
    """read the tabulated freq and time grids and write to local files

    Args:
        ngrids (float): number of grids
        erange (float): energy range

    Returns:
        freq and time grids.
        (2*ngrids,) array, first ngrids as points and left are weights.
    """
    path = os.path.dirname(__file__) + "/"

    grids = np.zeros((2, 2*ngrids))

    fn_freq = str(ngrids)+"_freq_points.dat"
    fn_time = str(ngrids)+"_time_points.dat"
    for ifn, fn in enumerate([fn_freq, fn_time]):
        try:
            with open(path+fn, 'r') as ff:
                while True:
                    lines = ff.readline()
                    if lines[:6] == 'Erange':
                        _, lb, ub = lines.split()
                        if eval(lb) <= erange <= eval(ub):
                            lines = next(ff)
                            i = 0
                            while i < 2*ngrids:
                                lines = next(ff).strip()
                                grids[ifn, i] = eval(lines)
                                i += 1
                            break
        except fnerror:
            raise fnerror("%d-point freq/time grids is not found" % ngrids)
    for ifn, fn in enumerate([fn_freq, fn_time]):
        with open('local_'+fn, 'w') as fw:
            print("\n".join(str(x) for x in grids[ifn, :].tolist()), file=fw)

    return grids[0, :], grids[1, :]


def generate_x(ngrids, erange):
    num_points_per_magnitude = 200
    num_x_node = (int(math.log10(erange))+1)*num_points_per_magnitude
    num_x_node = max(num_x_node, ngrids)

    x = np.power((erange)**(1.0/(num_x_node-1.0)),
                 np.arange(0, num_x_node))
    return x


def generate_trans_mat(omegas, taus, erange, trans, kind):
    """generate and write the cosin/sine transformation matrix

    Args:
        omegas, taus: array generated from ``read_grid``
        erange (float): energy range
        trans: "cos" or "sin"
        kind: "t2f" or "f2t"
    """
    assert trans in ["cos", "sin"]
    assert kind in ["t2f", "f2t"]

    ngrids = len(taus) // 2
    x = generate_x(ngrids, erange)
    num_x_node = len(x)
    omega_g = omegas[:ngrids]
    tau_g = taus[:ngrids]
    p0 = np.ones(ngrids)

    trans_str = "Cosine"
    if trans == "sin":
        trans_str = "Sine"
    kind_str = "time2freq"
    if kind == "f2t":
        kind_str = "freq2time"
    dict_target = {
            ("t2f", "cos"): target_t2f_cos,
            ("t2f", "sin"): target_t2f_sin,
            ("f2t", "cos"): target_f2t_cos,
            ("f2t", "sin"): target_f2t_sin,
            }
    tfunc = dict_target[(kind, trans)]

    mat = np.zeros((ngrids, ngrids))
    with open("INFO.txt", "at") as f1:
        print(kind_str + "_transform   N = ", ngrids,
              "   Range_x: ", np.min(x), np.max(x),
              "n(x) = ", num_x_node, file=f1)
        print(trans_str + " transform:", file=f1)
    fn = str(ngrids) + "_" + kind_str + "_grid_" + trans + ".txt"
    with open(fn, "wt") as f0:
        if kind == "t2f":
            for i_tau in omega_g:
                print("Freq: ", i_tau, file=f0)
                para_i = least_squares(tfunc, p0, args=(x, i_tau, tau_g))
                for gamma in para_i.x:
                    print(gamma, file=f0)
                with open("INFO.txt", "at") as f1:
                    print("Freq_point: ", i_tau, "         transform_err: ",
                          np.max(abs(para_i.fun)), file=f1)
        else:
            for i_tau in tau_g:
                print("Tau: ", i_tau, file=f0)
                para_i = least_squares(tfunc, p0, args=(x, i_tau, omega_g))
                for gamma in para_i.x:
                    print(gamma, file=f0)
                with open("INFO.txt", "at") as f1:
                    print("Tau_point: ", i_tau, "         transform_err: ",
                          np.max(abs(para_i.fun)), file=f1)

    with open("INFO.txt", "at") as f1:
        print("Finish " + kind_str + " " + trans_str + " transformation!", file=f1)
    return mat


def main():
    """main stream"""
    ngrids = int(sys.argv[1])
    erange = float(sys.argv[2])
    print("N: ", ngrids)
    print("R: ", erange)
    omegas, taus = read_grid(ngrids, erange)
    # TODO: internal check for Delta matrix?
    generate_trans_mat(omegas, taus, erange, "cos", "t2f")
    generate_trans_mat(omegas, taus, erange, "sin", "t2f")
    generate_trans_mat(omegas, taus, erange, "cos", "f2t")
    generate_trans_mat(omegas, taus, erange, "sin", "f2t")


if __name__ == '__main__':
    main()
