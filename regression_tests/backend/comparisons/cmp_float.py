"""Compare sequences of extracted scalar values, grouped by filename."""

import math


__all__ = ["abs_diff"]


def abs_diff(tolerance, precision=3):
    """Create an absolute-difference comparator for scalar sequences.

    Each input maps filenames to sequences of values convertible to float.
    A single string is treated as a one-value sequence; numeric strings may
    use Fortran D exponents. Corresponding values are compared in order across
    matching files, and the largest absolute difference determines the result.
    Missing files, unequal sequence lengths, and empty inputs fail. Matching
    NaN pairs are accepted, while a NaN on only one side fails.

    Args:
        tolerance: Maximum allowed absolute difference, inclusive.
        precision: Digits after the decimal point in scientific-notation
            diagnostics. Defaults to 3; does not affect comparison accuracy.

    Returns:
        A callable accepting test and reference dictionaries and returning
        (passed, message). The message reports the maximum difference and its
        filename and one-based value index, or explains a structural mismatch.

    Examples:
        >>> compare = abs_diff(1e-4)
        >>> compare({}, {})
        (False, 'no files found')
        >>> compare({"librpa.out": ["1.00001"]},
        ...         {"librpa.out": ["1.00002"]})[0]
        True
        >>> compare({"librpa.out": ["nan"]},
        ...         {"librpa.out": ["NaN"]})[0]
        True
        >>> compare({"librpa.out": ["nan"]}, {"librpa.out": ["1.0"]})
        (False, 'nan mismatch in librpa.out value 1: nan != 1.0')
        >>> compare({"librpa.out": ["1.0"]}, {"librpa.out": ["nan"]})
        (False, 'nan mismatch in librpa.out value 1: 1.0 != nan')
        >>> compare({"librpa.out": []}, {"librpa.out": ["-1.234"]})
        (False, 'value count mismatch in librpa.out: 0 != 1')
        >>> compare({"librpa.out": []}, {"librpa.out": []})
        (False, 'no scalar values found in librpa.out')
    """
    tolerance = float(tolerance)
    prec = int(precision)

    def inner_float_absdiff(fnobj1, fnobj2):
        """inner function of closure"""
        msg = r"max abs diff = {:.%iE} (tol = {:.%iE})" % (prec, prec)
        diff = 0.0
        diff_loc = None
        fns = set([*fnobj1.keys(), *fnobj2.keys()])
        if not fns:
            return False, "no files found"

        for fn in fns:
            try:
                objs1 = _as_objects(fnobj1[fn])
                objs2 = _as_objects(fnobj2[fn])
            except KeyError:
                return False, "missing file {}".format(fn)

            if len(objs1) == 0 and len(objs2) == 0:
                return False, "no scalar values found in {}".format(fn)
            if len(objs1) != len(objs2):
                return False, "value count mismatch in {}: {} != {}".format(
                    fn, len(objs1), len(objs2)
                )

            for iobj, (o1, o2) in enumerate(zip(objs1, objs2), 1):
                value1 = _parse_float(o1)
                value2 = _parse_float(o2)
                if math.isnan(value1) or math.isnan(value2):
                    if math.isnan(value1) and math.isnan(value2):
                        continue
                    return False, (
                        "nan mismatch in {} value {}: {} != {}"
                        .format(fn, iobj, o1, o2)
                    )
                d = abs(value1 - value2)
                if d > diff:
                    diff = d
                    diff_loc = (fn, iobj)

        msg = msg.format(diff, tolerance)
        if diff_loc is not None:
            fn, iobj = diff_loc
            msg += ", max at {} value {}".format(fn, iobj)
        return diff <= tolerance, msg
    return inner_float_absdiff


def _as_objects(raw):
    if isinstance(raw, str):
        return [raw]
    return list(raw)


def _parse_float(value):
    if isinstance(value, str):
        value = value.replace("D", "E").replace("d", "E")
    return float(value)
