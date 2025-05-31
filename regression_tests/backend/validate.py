import pathlib
import re
import importlib
from typing import Tuple


def _extract_plain(directory: str, file: str,
                   regex: re.Pattern, headers: int, rows: int, occurences: Tuple[int, int, int]):
    """Extract plain text data"""
    raw = {}
    d = pathlib.Path(directory)
    for f in d.rglob(file):
        rela = f.relative_to(d)
        with open(f, 'r') as h:
            lines = h.readlines()
        # When regex is None, treat the whole file as necessary data
        if regex is None:
            if rows is None:
                raw[rela] = "".join(lines[headers:])
            else:
                raw[rela] = "".join(lines[headers:rows])
            return

        st, ed, intv = occurences
        matched = []
        for i, l in enumerate(lines):
            m = regex.search(l)
            if m is not None:
                # the matched line have the data
                if rows is None:
                    if len(m.groups()) == 1:
                        matched.append(m.group(1))
                    else:
                        matched.append(" ".join(m.groups()))
                else:
                    matched.append("".join(lines[i + headers: i + rows]))
        if ed is None:
            ed = len(matched) - 1
        raw[rela] = [matched[i] for i in range(st, ed + 1, intv)]
    return raw


def _process_regex(regex: str):
    if regex is not None:
        regex = re.compile(regex)
    return regex


def _import_comparison(comparison: str):
    """
    import a function from a given module and construct the closure from
    the additional arguments, if any
    """
    REGEX_IMPORTCHECK = re.compile(r"(?P<mod>[A-z0-9_]+).(?P<func>[A-z0-9_]+)"
                                   r"[ ]*(?:\((?P<args>.*)\)|)$")
    try:
        match = REGEX_IMPORTCHECK.match(comparison)
        modname = "backend.comparisons.%s" % (match.group("mod"))
        module = importlib.import_module(modname)
        if match.group("args") is not None:
            paramlist = match.group("args").split(",")
            args = tuple(x for x in paramlist if "=" not in x)
            keywords = tuple(x.split("=") for x in paramlist if "=" in x)
            keywords = {x[0].strip(): x[1] for x in keywords}
            return getattr(module, match.group("func"))(*args, **keywords)
        else:
            return getattr(module, match.group("func"))
    except ImportError:
        raise ImportError("Failed to import comparision method: " + comparison)


def _import_binary_extract(binary_extract: str):
    return None


def _process_rows(rows: str):
    if rows is not None:
        rows = int(rows)
    return rows


def _process_headers(headers: str):
    if headers is not None:
        headers = int(headers)
    # No header lines
    return 0


def _process_occurences(occurences: str):
    st, ed, intv = 0, None, 1
    if occurences is not None:
        words = occurences.split(":")
        if len(words) == 3:
            if words[0]:
                st = int(words[0])
            if words[1]:
                ed = int(words[1])
            if words[2]:
                intv = int(words[2])
        elif len(words) == 2:
            if words[0] == "":
                st = 0
                ed = int(words[1])
            elif words[1] == "":
                st = words[0]
                ed = None
            else:
                st, ed = map(int, words)
        elif len(words) == 1:
            # single occurence
            st = int(words[0])
            ed = int(words[0])
        else:
            raise ValueError("invalid occurences")
    return st, ed, intv


class Validate():

    def __init__(self, name: str, file: str, comparison: str, headers: str, rows: str,
                 regex: str, occurences: str, binary_extract: str):
        self._name = name
        self._file = file
        self._comparison = _import_comparison(comparison)
        self._headers = _process_headers(headers)
        self._rows = _process_rows(rows)
        self._regex = _process_regex(regex)
        self._occurences = _process_occurences(occurences)
        self._binary_extract = _import_binary_extract(binary_extract)

    def evaluate(self, dir_test, dir_refr):
        if self._binary_extract is None:
            test = _extract_plain(dir_test, self._file, self._regex, self._headers, self._rows, self._occurences)
            refr = _extract_plain(dir_refr, self._file, self._regex, self._headers, self._rows, self._occurences)
        else:
            test = self._binary_extract(dir_test, self._file)
            refr = self._binary_extract(dir_refr, self._file)

        # print(dir_test, dir_refr)
        return self._comparison(test, refr)
