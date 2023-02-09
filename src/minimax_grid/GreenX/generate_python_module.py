#!/usr/bin/env python
from __future__ import print_function
import warnings
from argparse import ArgumentParser
import numpy as np

FUNCTION_NAME_TEMPLATE = "_init_minimax_%s_grid_%d"
WRAPPER_FUNCTION_NAME = "init_minimax_grid"


def get_grids_weights_from_datafile(datafile):
    """read the Minimax data file"""
    with open(datafile, 'r') as h:
        lines = h.readlines()
    i = 0
    try:
        if lines[0].startswith("###"):
            ngrids = int(lines[3].split()[2])
            lines = lines[18:]
        else:
            ngrids = int(lines[1].split()[2])
            lines = lines[15:]
        data = []
        while True:
            try:
                l = lines[i]
            except IndexError:
                break
            if l.strip() == '':
                i += 1
                continue
            if l.strip().startswith("Erange"):
                eratio = tuple(map(float, l.split()[1:]))
                lines_eratio = lines[i + 2:i + 2 + 2 * ngrids]
                grids = tuple(map(float, [x.strip() for x in lines_eratio[:ngrids]]))
                grids = np.array(grids, dtype="float128")
                weights = tuple(
                    map(float, [
                        x.strip() for x in lines_eratio[ngrids:2 * ngrids]
                    ]))
                weights = np.array(weights, dtype="float128")
                data.append({"eratio": eratio, "grids": grids, "weights": weights})
                i += 2 * ngrids + 2
            else:
                i += 1
    except ValueError:
        msg = "Invalid value around line %d in datafile %s" % (i + 1, datafile)
        warnings.warn(msg)
        return None, None

    # check whether eratio is valid
    # check eratio[1] of dataset is eratio[0] of the next set
    for i in range(len(data) - 1):
        if data[i]["eratio"][0] > data[i]["eratio"][1]:
            raise ValueError("eratio[0] < eratio[1] in set %d, datafile %s" % (i, datafile))
        if abs(data[i]["eratio"][1] - data[i + 1]["eratio"][0]) > 1:
            raise ValueError(
                "eratio[1] of set %d != eratio[0] of set %d, %f != %f, datafile %s" %
                (i, i + 1, data[i]["eratio"][1], data[i + 1]["eratio"][0], datafile))
    return ngrids, data


def format_tf_module(ngrids, data, keyword, indent=4, doubletype="float64"):
    """base function to format the fortran subroutine"""
    header = [
        "def " + FUNCTION_NAME_TEMPLATE % (keyword, ngrids) + "(eratio, scale):",
    ]
    variables = []

    assignments = []
    for i in range(len(data)):
        eratio = data[i]["eratio"]
        grids = data[i]["grids"]
        weights = data[i]["weights"]

        if i == 0:
            assignments.append("if eratio <= %.4f:" % eratio[1])
        elif i == len(data) - 1:
            assignments.append("else:")
        else:
            assignments.append("elif eratio < %.4f:" % eratio[1])

        assignments.append(" " * indent + "%s_grids = [" % keyword)

        for grid in grids:
            assignments.append(" " * 2 * indent + "%.16e," % grid)
        assignments.append(" " * indent + "]")

        assignments.append(" " * indent + "%s_weights = [" % keyword)
        for grid in weights:
            assignments.append(" " * 2 * indent + "%.16e," % grid)
        assignments.append(" " * indent + "]")

    postprocess = [
        "%s_grids = np.array(%s_grids, dtype='%s') * scale" % (keyword, keyword, doubletype),
        "%s_weights = np.array(%s_weights, dtype='%s') * scale" % (keyword, keyword, doubletype),
        "return %s_grids, %s_weights" % (keyword, keyword),
    ]

    body = []
    for x in variables + assignments + postprocess:
        if x.strip() == '':
            body.append(x)
        else:
            body.append(" " * indent + x)
    footer = []
    return header + body + footer


def format_freq_function(ngrids, freqdata, indent=4, doubletype="float64"):
    """format the fortran subroutine for a frequency ngrids"""
    return format_tf_module(ngrids, freqdata, "freq", indent, doubletype)


def format_time_function(ngrids, timedata, indent=4, doubletype="float64"):
    """format the fortran subroutine for a frequency ngrids"""
    return format_tf_module(ngrids, timedata, "time", indent, doubletype)


def format_wrapper_function(ngrids_list, indent=4, doubletype="float64"):
    """format the wrapper subroutine of a list of ngrids"""
    header = [
        "def " + WRAPPER_FUNCTION_NAME + "(ngrids, emin, emax, no_scale=False):"
    ]
    variables = []
    assignments = [
        "",
        "if emin < 0.0:",
        " " * indent + "raise ValueError('Metal system not implemented')",
        "",
        "eratio = emax / emin",
        "scale_freq = emin",
        "if no_scale:",
        " " * indent + "scale_freq = 1.0",
        "scale_time = 1.0 / scale_freq",
        ""
    ]

    selection = [
        "grid_func_dict = {",
    ]
    for ngrids in ngrids_list:
        selection.append(" " * indent + "%d: " % ngrids + "(" + FUNCTION_NAME_TEMPLATE %
                         ("freq", ngrids) + ", " + FUNCTION_NAME_TEMPLATE %
                         ("time", ngrids) + "),")
    selection.append("}")
    selection.extend([
        "if ngrids not in grid_func_dict:",
        " " * indent + "raise ValueError('grid size not supported')",
        "freq_func, time_func = grid_func_dict[ngrids]",
        "freq_grids, freq_weights = freq_func(eratio, scale_freq)",
        "time_grids, time_weights = time_func(eratio, scale_time)",
    ])

    ret = ["return freq_grids, freq_weights, time_grids, time_weights", ]

    body = []
    for x in variables + assignments + selection + ret:
        if x.strip() == '':
            body.append(x)
        else:
            body.append(" " * indent + x)
    footer = []
    return header + body + footer


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--freqdata", nargs="+", type=str, help="filenames of frequency grids")
    parser.add_argument("--timedata", nargs="+", type=str, help="filenames of time grids")
    parser.add_argument("--indent", default=4, type=int, help="indent inside subroutines")
    parser.add_argument("--dt", default="float64",
                        type=str, choices=["float64", "float128"],
                        help="indent inside subroutines")
    parser.add_argument("--disable-preprocess-freq-weight", action="store_true", help="")

    args = parser.parse_args()
    freq_functions = []
    time_functions = []
    freq_ngrids_list = []
    time_ngrids_list = []
    preproc_weight_scale = 0.25
    if args.disable_preprocess_freq_weight:
        preproc_weight_scale = None
    preproc_ngrids_list = list(range(1, 25))

    if args.freqdata is not None:
        for freqfile in args.freqdata:
            ngrids, data = get_grids_weights_from_datafile(freqfile)
            if ngrids is None:
                continue
            if ngrids in preproc_ngrids_list and preproc_weight_scale is not None:
                for d in data:
                    d["weights"] *= preproc_weight_scale
            freq_ngrids_list.append(ngrids)
            freq_functions.append(
                format_freq_function(ngrids,
                                     data,
                                     args.indent,
                                     doubletype=args.dt))

    if args.timedata is not None:
        for timefile in args.timedata:
            ngrids, data = get_grids_weights_from_datafile(timefile)
            if ngrids is None:
                continue
            time_ngrids_list.append(ngrids)
            time_functions.append(
                format_time_function(ngrids,
                                     data,
                                     args.indent,
                                     doubletype=args.dt))

    functions = []
    ngrids_list = []
    for ngrids in freq_ngrids_list:
        if ngrids in time_ngrids_list:
            ngrids_list.append(ngrids)
            functions.append(freq_functions[freq_ngrids_list.index(ngrids)])
            functions.append(time_functions[time_ngrids_list.index(ngrids)])

    module_imports = ["import numpy as np", ""]
    module_headers = ["__all__ = ['%s', 'available_minimax_grid ']" % WRAPPER_FUNCTION_NAME, "",
                      ""]
    module_body = ["available_minimax_grid = (" + ", ".join(str(x) for x in ngrids_list) + ")", "", ""]

    for sub in [
            format_wrapper_function(
                ngrids_list, args.indent, doubletype=args.dt),
    ] + functions:
        for s in sub:
            if s.strip() == '':
                module_body.append("")
            else:
                module_body.append(s)
        module_body.append("")
        module_body.append("")

    print("\n".join(module_imports + module_headers + module_body))
