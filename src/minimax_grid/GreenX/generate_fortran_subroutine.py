#!/usr/bin/env python
from __future__ import print_function
import warnings
from argparse import ArgumentParser
import numpy as np

SUBROUTINE_NAME_TEMPLATE = "init_minimax_%s_grid_%d"
WRAPPER_SUBROUTINE_NAME = "init_minimax_grid"


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


def format_tf_subroutine(ngrids, data, keyword, indent=4, doubletype="real*8"):
    """base function to format the fortran subroutine"""
    header = [
        "subroutine " + SUBROUTINE_NAME_TEMPLATE % (keyword, ngrids)
        + "(%s_grids, %s_weights, eratio, scale)" % (keyword, keyword),
    ]
    variables = [
        "",
        "%s,intent(out) :: %s_grids(%d)" % (doubletype, keyword, ngrids),
        "%s,intent(out) :: %s_weights(%d)" % (doubletype, keyword, ngrids),
        "%s,intent(in)  :: eratio" % doubletype,
        "%s,intent(in)  :: scale" % doubletype,
        "",
    ]

    assignments = []
    for i in range(len(data)):
        eratio = data[i]["eratio"]
        grids = data[i]["grids"]
        weights = data[i]["weights"]

        if i == 0:
            assignments.append("if (eratio .le. %.4f) then" % eratio[1])
        elif i == len(data) - 1:
            assignments.append("else")
        else:
            assignments.append("else if (eratio .lt. %.4f) then" % eratio[1])

        assignments.append(" " * indent + "%s_grids(:) = (/ &" % keyword)

        for grid in grids[:-1]:
            assignments.append(" " * 2 * indent + "%.16e, &" % grid)
        assignments.append(" " * 2 * indent + "%.16e &" % grids[-1])
        assignments.append(" " * indent + "/)")

        assignments.append(" " * indent + "%s_weights(:) = (/ &" % keyword)
        for grid in weights[:-1]:
            assignments.append(" " * 2 * indent + "%.16e, &" % grid)
        assignments.append(" " * 2 * indent + "%.16e &" % weights[-1])
        assignments.append(" " * indent + "/)")

        if i == len(data) - 1:
            assignments.append("end if")

    postprocess = [
        "%s_grids(:) = %s_grids(:) * scale" % (keyword, keyword),
        "%s_weights(:) = %s_weights(:) * scale" % (keyword, keyword)
    ]

    body = []
    for x in variables + assignments + postprocess:
        if x.strip() == '':
            body.append(x)
        else:
            body.append(" " * indent + x)
    footer = ["end subroutine"]
    return header + body + footer


def format_freq_subroutine(ngrids, freqdata, indent=4, doubletype="real*8"):
    """format the fortran subroutine for a frequency ngrids"""
    return format_tf_subroutine(ngrids, freqdata, "freq", indent, doubletype)


def format_time_subroutine(ngrids, timedata, indent=4, doubletype="real*8"):
    """format the fortran subroutine for a frequency ngrids"""
    return format_tf_subroutine(ngrids, timedata, "time", indent, doubletype)


def format_wrapper_subroutine(ngrids_list, indent=4, doubletype="real*8"):
    """format the wrapper subroutine of a list of ngrids"""
    header = [
        "subroutine " + WRAPPER_SUBROUTINE_NAME +
        "(ngrids, emin, emax, freq_grids, freq_weights, time_grids, time_weights, ierr)"
    ]
    variables = [
        "",
        "integer,intent(in)  :: ngrids",
        "integer,intent(out) :: ierr",
        "%s,intent(in)        :: emin" % doubletype,
        "%s,intent(in)        :: emax" % doubletype,
        "%s,intent(out)       :: freq_grids(ngrids)" % doubletype,
        "%s,intent(out)       :: freq_weights(ngrids)" % doubletype,
        "%s,intent(out)       :: time_grids(ngrids)" % doubletype,
        "%s,intent(out)       :: time_weights(ngrids)" % doubletype,
        "",
        "%s :: eratio" % doubletype,
        "%s :: scale_freq" % doubletype,
        "%s :: scale_time" % doubletype,
        "",
    ]
    assignments = [
        "",
        "if (emin .le. 0.0) then",
        " " * indent + "ierr = 2",
        " " * indent + "return",
        "end if",
        "",
        "eratio = emax / emin",
        "scale_freq = emin",
        "scale_time = 1.0 / scale_freq",
        "ierr = 0",
        ""
    ]

    selection = [
        "select case (ngrids)",
    ]
    for ngrids in ngrids_list:
        selection.extend([
            "case (%d)" % ngrids,
            " " * indent + "call " + SUBROUTINE_NAME_TEMPLATE % ("freq", ngrids) + "(freq_grids, freq_weights, eratio, scale_freq)",
            " " * indent + "call " + SUBROUTINE_NAME_TEMPLATE % ("time", ngrids) + "(time_grids, time_weights, eratio, scale_time)",
        ])
    selection.extend(["case default", " " * indent + "ierr = 1", "end select"])

    body = []
    for x in variables + assignments + selection:
        if x.strip() == '':
            body.append(x)
        else:
            body.append(" " * indent + x)
    footer = ["end subroutine"]
    return header + body + footer


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--freqdata", nargs="+", type=str, help="filenames of frequency grids")
    parser.add_argument("--timedata", nargs="+", type=str, help="filenames of time grids")
    parser.add_argument("--indent", default=4, type=int, help="indent inside subroutines")
    parser.add_argument("--dt", default="double precision",
                        type=str, choices=["double precision", "real*8"],
                        help="indent inside subroutines")
    parser.add_argument("--disable-preprocess-freq-weight", action="store_true", help="")

    args = parser.parse_args()
    freq_subroutines = []
    time_subroutines = []
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
            freq_subroutines.append(
                format_freq_subroutine(ngrids,
                                       data,
                                       args.indent,
                                       doubletype=args.dt))

    if args.timedata is not None:
        for timefile in args.timedata:
            ngrids, data = get_grids_weights_from_datafile(timefile)
            if ngrids is None:
                continue
            time_ngrids_list.append(ngrids)
            time_subroutines.append(
                format_time_subroutine(ngrids,
                                       data,
                                       args.indent,
                                       doubletype=args.dt))

    subroutines = []
    ngrids_list = []
    for ngrids in freq_ngrids_list:
        if ngrids in time_ngrids_list:
            ngrids_list.append(ngrids)
            subroutines.append(freq_subroutines[freq_ngrids_list.index(ngrids)])
            subroutines.append(time_subroutines[time_ngrids_list.index(ngrids)])

    module_header = ["module mod_minimax_grids", ]
    module_declare = [
        "", "private",
        "public :: %s" % WRAPPER_SUBROUTINE_NAME,
        "", "contains", ""
    ]
    module_footer = [
        "end module",
    ]

    module_body = []
    for s in module_declare:
        if s.strip() == '':
            module_body.append("")
        else:
            module_body.append(" " * args.indent + s)

    for sub in [
            format_wrapper_subroutine(
                ngrids_list, args.indent, doubletype=args.dt),
    ] + subroutines:
        for s in sub:
            if s.strip() == '':
                module_body.append("")
            else:
                module_body.append(" " * args.indent + s)
        module_body.append("")

    print("\n".join(module_header + module_body + module_footer))
