import math
import re


__all__ = ["abs_diff"]


def abs_diff(tolerance, precision=3, columns="all", allow_matching_nan=False):
    """
    Compare numeric text tables with an absolute tolerance.

    Tables are expected to be extracted by Validate, typically with regex,
    headers, and rows selecting the table body after a matched header.
    Selected cells must be finite by default. A documented legacy fixture may
    explicitly allow matching NaNs; one-sided NaNs and infinities still fail.
    """
    tolerance = float(tolerance)
    if not math.isfinite(tolerance) or tolerance < 0:
        raise ValueError("tolerance must be finite and nonnegative")
    precision = int(precision)
    columns = _process_columns(columns)
    allow_matching_nan = str(allow_matching_nan).strip("\"'").lower() == "true"

    def inner_table_absdiff(fnobj1, fnobj2):
        msg = r"max abs diff = {:.%iE} (tol = {:.%iE}) over {} table rows" % (
            precision, precision
        )
        diff = 0.0
        diff_loc = None
        nrows = 0
        matching_nan = 0

        fns = set([*fnobj1.keys(), *fnobj2.keys()])
        if not fns:
            return False, "no files found"

        for fn in fns:
            try:
                tables1 = _as_tables(fnobj1[fn])
                tables2 = _as_tables(fnobj2[fn])
            except KeyError:
                return False, "missing file {}".format(fn)

            if len(tables1) == 0 and len(tables2) == 0:
                return False, "no tables found in {}".format(fn)
            if len(tables1) != len(tables2):
                return False, "table count mismatch in {}: {} != {}".format(
                    fn, len(tables1), len(tables2)
                )

            for itable, (table1, table2) in enumerate(zip(tables1, tables2), 1):
                if len(table1) != len(table2):
                    return False, "row count mismatch in {} table {}: {} != {}".format(
                        fn, itable, len(table1), len(table2)
                    )
                if len(table1) == 0:
                    return False, "no table rows found in {} table {}".format(
                        fn, itable
                    )

                ncells = 0
                for irow, (row1, row2) in enumerate(zip(table1, table2), 1):
                    if len(row1) != len(row2):
                        return False, (
                            "column count mismatch in {} table {} row {}: {} != {}"
                            .format(fn, itable, irow, len(row1), len(row2))
                        )
                    cols = _columns_for_row(columns, len(row1))
                    for icol in cols:
                        value1 = _parse_float(row1[icol])
                        value2 = _parse_float(row2[icol])
                        if math.isnan(value1) or math.isnan(value2):
                            if allow_matching_nan and math.isnan(value1) and math.isnan(value2):
                                matching_nan += 1
                                continue
                            return False, (
                                "nonfinite value (nan mismatch) in {} table {} row {} column {}: {} != {}"
                                .format(fn, itable, irow, icol + 1, row1[icol], row2[icol])
                            )
                        if not math.isfinite(value1) or not math.isfinite(value2):
                            return False, (
                                "nonfinite value in {} table {} row {} column {}: {} != {}"
                                .format(fn, itable, irow, icol + 1, row1[icol], row2[icol])
                            )
                        d = abs(value1 - value2)
                        if d > diff:
                            diff = d
                            diff_loc = (fn, itable, irow, icol + 1)
                        ncells += 1
                    nrows += 1
                if ncells == 0:
                    return False, "no table cells compared in {} table {}".format(
                        fn, itable
                    )

        msg = msg.format(diff, tolerance, nrows)
        if matching_nan:
            msg += "; legacy exception: {} matching NaN cells skipped".format(matching_nan)
        if diff_loc is not None:
            fn, itable, irow, icol = diff_loc
            msg += ", max at {} table {} row {} column {}".format(
                fn, itable, irow, icol
            )
        return diff <= tolerance, msg

    return inner_table_absdiff


def _as_tables(raw):
    if isinstance(raw, str):
        raw = [raw]
    tables = []
    for item in raw:
        rows = []
        for line in item.splitlines():
            if line.strip():
                rows.append(line.split())
        tables.append(rows)
    return tables


def _parse_float(value):
    return float(value.replace("D", "E").replace("d", "E"))


def _process_columns(columns):
    if columns is None:
        return None
    columns = columns.strip()
    if columns.lower() in ("all", "*", "table"):
        return None

    parsed = []
    for item in re.split(r"[\s;]+", columns):
        if not item:
            continue
        if ":" in item:
            parsed.append(_parse_range(item))
        else:
            parsed.append(int(item) - 1)
    return tuple(parsed)


def _parse_range(item):
    start, stop = item.split(":", 1)
    start = int(start) if start else 1
    stop = int(stop) if stop else 0
    return start - 1, stop


def _columns_for_row(columns, ncols):
    if columns is None:
        return range(ncols)
    expanded = []
    for spec in columns:
        if isinstance(spec, tuple):
            start, stop = spec
            if stop == 0:
                stop = ncols
            expanded.extend(range(start, stop))
        else:
            expanded.append(spec)
    for col in expanded:
        if col < 0 or col >= ncols:
            raise ValueError("column {} is outside table width {}".format(col + 1, ncols))
    return expanded
