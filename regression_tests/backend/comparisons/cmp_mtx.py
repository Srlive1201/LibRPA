import math


__all__ = ["abs_diff"]


def abs_diff(tolerance, precision=3, comments=False):
    tolerance = float(tolerance)
    if not math.isfinite(tolerance) or tolerance < 0:
        raise ValueError("tolerance must be finite and nonnegative")
    precision = int(precision)
    comments = _as_bool(comments)

    def inner_mtx_absdiff(fnobj1, fnobj2):
        msg = r"max abs diff = {:.%iE} (tol = {:.%iE}) over {} matrix entries" % (
            precision, precision
        )
        diff = 0.0
        diff_loc = None
        nentries = 0

        fns = set([*fnobj1.keys(), *fnobj2.keys()])
        if not fns:
            return False, "no files found"

        for fn in fns:
            try:
                matrices1 = _as_matrices(fnobj1[fn])
                matrices2 = _as_matrices(fnobj2[fn])
            except KeyError:
                return False, "missing file {}".format(fn)

            if len(matrices1) == 0 and len(matrices2) == 0:
                return False, "no matrices found in {}".format(fn)
            if len(matrices1) != len(matrices2):
                return False, "matrix count mismatch in {}: {} != {}".format(
                    fn, len(matrices1), len(matrices2)
                )

            for imatrix, (matrix1, matrix2) in enumerate(zip(matrices1, matrices2), 1):
                meta1, comments1, entries1 = matrix1
                meta2, comments2, entries2 = matrix2
                if meta1 != meta2:
                    return False, "matrix header mismatch in {} matrix {}: {} != {}".format(
                        fn, imatrix, meta1, meta2
                    )
                if comments and comments1 != comments2:
                    return False, "matrix comments mismatch in {} matrix {}".format(
                        fn, imatrix
                    )
                if len(entries1) != len(entries2):
                    return False, "entry count mismatch in {} matrix {}: {} != {}".format(
                        fn, imatrix, len(entries1), len(entries2)
                    )

                for ientry, (entry1, entry2) in enumerate(zip(entries1, entries2), 1):
                    row1, col1, value1 = entry1
                    row2, col2, value2 = entry2
                    if (row1, col1) != (row2, col2):
                        return False, (
                            "coordinate mismatch in {} matrix {} entry {}: ({}, {}) != ({}, {})"
                            .format(fn, imatrix, ientry, row1, col1, row2, col2)
                        )
                    ok, d = _value_diff(value1, value2)
                    if not ok:
                        return False, (
                            "nonfinite value (nan mismatch) in {} matrix {} row {} column {}"
                            .format(fn, imatrix, row1, col1)
                        )
                    if d > diff:
                        diff = d
                        diff_loc = (fn, imatrix, row1, col1)
                    nentries += 1

        msg = msg.format(diff, tolerance, nentries)
        if diff_loc is not None:
            fn, imatrix, row, col = diff_loc
            msg += ", max at {} matrix {} row {} column {}".format(
                fn, imatrix, row, col
            )
        return diff <= tolerance, msg

    return inner_mtx_absdiff


def _as_matrices(raw):
    if isinstance(raw, str):
        raw = [raw]
    return [_read_coordinate_matrix_market(item) for item in raw]


def _read_coordinate_matrix_market(text):
    lines = iter(text.splitlines())
    header = _next_line(lines)
    parts = header.split()
    if len(parts) != 5 or parts[0] != "%%MatrixMarket":
        raise ValueError("invalid Matrix Market header")
    _, obj, fmt, field, symmetry = [parts[0], *[p.lower() for p in parts[1:]]]
    if obj != "matrix" or fmt != "coordinate":
        raise ValueError("only Matrix Market coordinate matrices are supported")
    if field not in ("real", "integer", "complex", "pattern"):
        raise ValueError("unsupported Matrix Market field {}".format(field))

    size_line, comments = _next_data_line(lines)
    size = size_line.split()
    if len(size) != 3:
        raise ValueError("invalid Matrix Market coordinate size line")
    nrows, ncols, nnz = [int(x) for x in size]
    entries = []
    for line in lines:
        line = line.strip()
        if not line or line.startswith("%"):
            continue
        items = line.split()
        if len(items) != {"complex": 4, "pattern": 2}.get(field, 3):
            raise ValueError("invalid Matrix Market {} entry".format(field))
        row, col = int(items[0]), int(items[1])
        if field == "complex":
            value = (_parse_float(items[2]), _parse_float(items[3]))
        elif field == "pattern":
            value = (1.0,)
        else:
            value = (_parse_float(items[2]),)
        entries.append((row, col, value))

    if len(entries) != nnz:
        raise ValueError("Matrix Market entry count {} != {}".format(len(entries), nnz))
    return (nrows, ncols, field, symmetry), comments, sorted(entries)


def _next_line(lines):
    for line in lines:
        line = line.strip()
        if line:
            return line
    raise ValueError("unexpected end of Matrix Market file")


def _next_data_line(lines):
    comments = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith("%"):
            comments.append(line)
            continue
        return line, tuple(comments)
    raise ValueError("unexpected end of Matrix Market file")


def _as_bool(value):
    if isinstance(value, bool):
        return value
    return str(value).strip("\"'").lower() in ("1", "true", "yes", "on")


def _parse_float(value):
    return float(value.replace("D", "E").replace("d", "E"))


def _value_diff(value1, value2):
    diffs = []
    for x, y in zip(value1, value2):
        if not math.isfinite(x) or not math.isfinite(y):
            return False, 0.0
        diffs.append(x - y)
    if len(diffs) == 1:
        return True, abs(diffs[0])
    return True, math.hypot(*diffs)
