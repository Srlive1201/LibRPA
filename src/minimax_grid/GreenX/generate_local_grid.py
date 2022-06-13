"""write the minimax grids and cosine/sine transformation matrix"""
from __future__ import print_function
from argparse import ArgumentParser
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

    with open("INFO.txt", "at") as f1:
        print(kind_str + "_transform   N = ", ngrids,
              "   Range_x: ", np.min(x), np.max(x),
              "n(x) = ", num_x_node, file=f1)
        print(trans_str + " transform:", file=f1)
    fn = str(ngrids) + "_" + kind_str + "_grid_" + trans + ".txt"
    mat = np.zeros((ngrids, ngrids))
    with open(fn, "wt") as f0:
        if kind == "t2f":
            prefix = "Freq"
            outer = omega_g
            inner = tau_g
        else:
            prefix = "Tau"
            outer = tau_g
            inner = omega_g
        for i, grid in enumerate(outer):
            print(prefix + ": ", grid, file=f0)
            para_i = least_squares(tfunc, p0, args=(x, grid, inner))
            for gamma in para_i.x:
                print(gamma, file=f0)
            with open("INFO.txt", "at") as f1:
                print(prefix+"_point: ", grid, "         transform_err: ",
                      np.max(abs(para_i.fun)), file=f1)
            mat[i, :] = para_i.x

    with open("INFO.txt", "at") as f1:
        print("Finish " + kind_str + " " + trans_str + " transformation!", file=f1)
    return mat


def parser():
    """the parser"""
    p = ArgumentParser()
    p.add_argument("ngrids", type=int, help="Number of grid points",
                   choices=[6, 8, 10, 12, 14, 16, 18, 20, 26, 28, 30, 32, 34])
    p.add_argument("erange", type=float, help="energy range")
    p.add_argument("--check-delta", action="store_true",
                   help="check the Delta matrix of cosine and sine transform")
    p.add_argument("--show", action="store_true")
    p.add_argument("--plot", action="store_true")
    p.add_argument("--dpi", type=int, default=300, help="DPI")
    return p

def main():
    """main stream"""
    args = parser().parse_args()
    ngrids = int(args.ngrids)
    erange = float(args.erange)
    print("N: ", ngrids)
    print("R: ", erange)
    omegas, taus = read_grid(ngrids, erange)
    cos_t2f = generate_trans_mat(omegas, taus, erange, "cos", "t2f")
    sin_t2f = generate_trans_mat(omegas, taus, erange, "sin", "t2f")
    cos_f2t = generate_trans_mat(omegas, taus, erange, "cos", "f2t")
    sin_f2t = generate_trans_mat(omegas, taus, erange, "sin", "f2t")

    # internal check for Delta matrix
    if args.check_delta:
        np.set_printoptions(precision=3, linewidth=200)
        import matplotlib.pyplot as plt
        fig, axs = plt.subplots(1, 2, figsize=(12, 6))
        fig.suptitle("t2f.f2t, ngrid = %d, erange = %f" % (ngrids, erange))
        for title, ax, mat in zip(["Cos", "Sin"],
                                  axs, [np.matmul(cos_t2f, cos_f2t),
                                  np.matmul(sin_t2f, sin_f2t)]):
            delta = np.absolute(np.identity(ngrids) - mat)
            # print(delta)
            norm = np.linalg.norm(delta)
            vmax = np.max(delta)
            print("max diff element: %f, 2-norm of %s trans: %f" % (vmax, title, norm))
            ax.set_title(title + " Delta, 2-norm: %.4f" % norm)
            c = ax.matshow(delta, cmap="Blues", vmin=0, vmax=vmax)
            fig.colorbar(c, ax=ax, fraction=0.046, pad=0.04)
        if args.plot:
            plt.savefig("Delta_n%d_e%.5f.png" % (ngrids, erange), dpi=args.dpi)
        if args.show:
            plt.show()


if __name__ == '__main__':
    main()
