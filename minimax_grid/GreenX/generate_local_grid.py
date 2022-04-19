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


def _fun(x, y):
    """1/(x+iy) + 1/(x-iy)"""
    return (2*x)/(x**2+y**2)


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
            with open(path+fn,'r') as ff:
                while True:
                    lines=ff.readline()
                    if lines[:6]=='Erange':
                        _, lb, ub = lines.split()
                        if eval(lb) <= erange <= eval(ub):
                            lines=next(ff)
                            i=0
                            while i<2*ngrids:
                                lines=next(ff).strip()
                                grids[ifn, i] = eval(lines)
                                i+=1
                            break
        except fnerror:
            raise fnerror("%d-point freq/time grids is not found" % ngrids)
    for ifn, fn in enumerate([fn_freq, fn_time]):
        with open('local_'+fn, 'w') as fw:
            print("\n".join(str(x) for x in grids[ifn, :].tolist()),file=fw)

    return grids[0, :], grids[1, :]


def time2freq_grid(omegas, taus, erange):
    """generate and write the cosin/sine transformation matrix

    Args:
        omegas, taus: array generated from ``read_grid``
        erange (float): energy range
    """
    def tau_fun(p,x,tau):
        s = 0
        i = 0
        while i< len(p):
            s += p[i]*np.exp(-x*abs(tau[i]))
            i += 1
        return s

    def error(p, x, y, tau):                    # 拟合残差
        return tau_fun(p, x, tau) - y

    ngrids = len(taus) // 2
    num_points_per_magnitude = 200
    num_x_node = (int(math.log10(erange))+1)*num_points_per_magnitude
    num_x_node = max(num_x_node, ngrids)

    x = np.power((erange)**(1.0/(num_x_node-1.0)),
                 np.arange(0, num_x_node))

    p0 = np.ones(ngrids)
    with open("INFO.txt", "at") as f1:
        print("Time_to_freq_transform   N = ", ngrids,
              "   Range_x: ", np.min(x), np.max(x),
              "n(x) = ", num_x_node, file=f1)
        print("Cosine transform:", file=f1)
    with open(str(ngrids)+"_time2freq_grid_cos.txt", "wt") as f0:
        omega_g=omegas[:ngrids]
        for i_freq in omega_g:
            print("Freq: ", i_freq, file=f0)
            y=_fun(x, i_freq)
            para_i=least_squares(error, p0, args=(x, y, taus))
            for gamma in para_i.x:
                print(gamma, file=f0)
            with open("INFO.txt", "at") as f1:
                print("Freq_point: ", i_freq, "         transform_err: ",
                      np.max(abs(para_i.fun)), file=f1)

    with open("INFO.txt", "at") as f1:
        print("", file=f1)
        print("Sine transform:", file=f1)
    with open(str(ngrids)+"_time2freq_grid_sin.txt", "wt") as f0:
        omega_g=omegas[:ngrids]
        for i_freq in omega_g:
            print("Freq: ", i_freq, file=f0)
            y=_fun(i_freq, x)
            para_i=least_squares(error, p0, args=(x, y, taus))
            for lambd in para_i.x:
                print(lambd, file=f0)
            with open("INFO.txt", "at") as f1:
                print("Freq_point: ", i_freq, "         transform_err: ",
                      np.max(abs(para_i.fun)), file=f1)

    with open("INFO.txt", "at") as f1:
        print("Finish!", file=f1)


def main():
    """main stream"""
    ngrids=int(sys.argv[1])
    erange=float(sys.argv[2])
    print("N: ", ngrids)
    print("R: ", erange)
    omegas, taus = read_grid(ngrids, erange)
    time2freq_grid(omegas, taus, erange)


if __name__=='__main__':
    main()

