import math

def abs_diff(tolerance, precision=3):
    """
    returns true, if the two objects can be converted to floats and are their
    absolute difference is smaller than the given tolerance. The optional
    precision argument determines the number of digits behind the decimal dot.
    """
    tolerance = float(tolerance)
    prec = int(precision)

    def inner_float_absdiff(fnobj1, fnobj2):
        """inner function of closure"""
        msg = "tol. "
        msg = r"max abs diff = {:.%iE} (tol = {:.%iE})" % (prec, prec)
        diff = 0.0
        fns = set([*fnobj1.keys(), *fnobj2.keys()])
        for fn in fns:
            try:
                objs1 = fnobj1[fn]
                objs2 = fnobj1[fn]
                for o1, o2 in zip(objs1, objs2):
                    diff = max(abs(float(o1) - float(o2)), diff)
            except KeyError:
                diff = max(math.nan, diff)
        return diff <= tolerance, msg.format(diff, tolerance)
    return inner_float_absdiff
