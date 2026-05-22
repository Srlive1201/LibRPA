import re


__all__ = ["abs_diff"]


def abs_diff(tolerance, precision=3, columns="all"):
    """
    Compare numeric text tables with an absolute tolerance.

    Tables are expected to be extracted by Validate, typically with regex,
    headers, and rows selecting the table body after a matched header.
    """
    tolerance = float(tolerance)
    precision = int(precision)
    columns = _process_columns(columns)

    def inner_table_absdiff(fnobj1, fnobj2):
        msg = r"max abs diff = {:.%iE} (tol = {:.%iE}) over {} table rows" % (
            precision, precision
        )
        diff = 0.0
        diff_loc = None
        nrows = 0

        fns = set([*fnobj1.keys(), *fnobj2.keys()])
        for fn in fns:
            try:
                tables1 = _as_tables(fnobj1[fn])
                tables2 = _as_tables(fnobj2[fn])
            except KeyError:
                return False, "missing file {}".format(fn)

            if len(tables1) != len(tables2):
                return False, "table count mismatch in {}: {} != {}".format(
                    fn, len(tables1), len(tables2)
                )

            for itable, (table1, table2) in enumerate(zip(tables1, tables2), 1):
                if len(table1) != len(table2):
                    return False, "row count mismatch in {} table {}: {} != {}".format(
                        fn, itable, len(table1), len(table2)
                    )

                for irow, (row1, row2) in enumerate(zip(table1, table2), 1):
                    if len(row1) != len(row2):
                        return False, (
                            "column count mismatch in {} table {} row {}: {} != {}"
                            .format(fn, itable, irow, len(row1), len(row2))
                        )
                    cols = _columns_for_row(columns, len(row1))
                    for icol in cols:
                        d = abs(_parse_float(row1[icol]) - _parse_float(row2[icol]))
                        if d > diff:
                            diff = d
                            diff_loc = (fn, itable, irow, icol + 1)
                    nrows += 1

        msg = msg.format(diff, tolerance, nrows)
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
