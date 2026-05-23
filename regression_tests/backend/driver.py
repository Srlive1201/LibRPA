import os
import pathlib
import shutil
import tarfile
from collections import OrderedDict
from typing import Tuple

from .utils import run_librpa
from .validate import Validate


__all__ = ["TestDriver"]

PASS_FAIL = {True: "PASS", False: "FAIL"}


def _disable_message(tc: dict):
    labels = tc.get("labels", {})
    disable = labels.get("disable", False)
    if disable is True:
        return ""
    if isinstance(disable, str):
        msg = disable.strip()
        if msg:
            return msg
    return None


def _disable_testcase(tc: dict, message):
    if _disable_message(tc) is not None:
        return
    tc.setdefault("labels", {})["disable"] = message


def _validate_scope_filter(testcases, option: str, values):
    selected = set(values or [])
    known = {tc["directory"] for tc in testcases}
    unknown = sorted(selected - known)
    if unknown:
        raise ValueError("Unknown test case directory in {:s}: {:s}".format(
            option, ", ".join(unknown)))
    return selected


def _check_build_run_filter(tc: dict,
                            ntasks: int, nthreads: int,
                            use_libri: bool) -> Tuple[bool, bool]:
    # Filter build options
    build = tc["build"]
    if build["require_libri"] and not use_libri:
        return True, False

    # Filter runtime options
    run = tc["run"]
    if ntasks in run["ntasks_disable"]:
        return False, True
    if nthreads in run["nthreads_disable"]:
        return False, True
    nt = run["ntasks_enable"]
    if len(nt) > 0 and ntasks not in nt:
        return False, True
    nt = run["nthreads_enable"]
    if len(nt) > 0 and nthreads not in nt:
        return False, True

    return False, False


class TestDriver:

    def __init__(self, dir_input: str, dir_ref: str, workspace: str, groups: dict):
        self._dir_input = pathlib.Path(dir_input)
        self._dir_ref = pathlib.Path(dir_ref)
        self._workspace = pathlib.Path(workspace)
        self._dir_testcase = self._workspace / "testcases"

        # Check
        if not self._dir_input.exists():
            raise FileNotFoundError("Input directory does not exist")
        if not self._dir_ref.exists():
            raise FileNotFoundError("Reference directory does not exist")

        self._groups = groups
        self._ntasks = None
        self._nthreads = None
        # Testcases that are qualified after initialized with build and runtime conditions
        self._testcases_filtered = None

    def initialize(self, ntasks: int, nthreads: int, use_libri: bool,
                   only=None, exclude=None):
        self._testcases_filtered = []
        self._ntasks = ntasks
        self._nthreads = nthreads

        if only and exclude:
            raise ValueError("--only and --exclude cannot be used together")

        all_testcases = [tc for gtcs in self._groups.values() for tc in gtcs]
        only = _validate_scope_filter(all_testcases, "--only", only)
        exclude = _validate_scope_filter(all_testcases, "--exclude", exclude)
        if only:
            for tc in all_testcases:
                if tc["directory"] not in only:
                    _disable_testcase(tc, "not selected by --only")
        elif exclude:
            for tc in all_testcases:
                if tc["directory"] in exclude:
                    _disable_testcase(tc, "excluded by --exclude")

        print("Initializing workspace targeting:")
        print("- MPI tasks :", ntasks)
        print("- threads   :", nthreads)
        print("- LibRI     ?", use_libri)

        # Filter test cases
        skip_due_to_disable = []
        skip_due_to_build = []
        skip_due_to_run = []

        for g, gtcs in self._groups.items():
            for tc in gtcs:
                if _disable_message(tc) is not None:
                    skip_due_to_disable.append(tc)
                    continue
                filter_build, filter_run = _check_build_run_filter(
                    tc, ntasks, nthreads, use_libri)
                if filter_build:
                    skip_due_to_build.append(tc)
                    continue
                if filter_run:
                    skip_due_to_run.append(tc)
                    continue
                self._testcases_filtered.append(tc)

        if skip_due_to_disable:
            print()
            print("Disabled following test cases")
            for tc in skip_due_to_disable:
                msg = _disable_message(tc)
                reminder = " [{:s}]".format(msg) if msg else ""
                print("- {:s}: {:s}{:s}".format(tc["directory"], tc["name"], reminder))

        if skip_due_to_build:
            print()
            print("Skipped following test cases due to executable build condition")
            for tc in skip_due_to_build:
                print("- {:s}: {:s}".format(tc["directory"], tc["name"]))

        if skip_due_to_run:
            print()
            print("Skipped following test cases due to runtime condition")
            for tc in skip_due_to_run:
                print("- {:s}: {:s}".format(tc["directory"], tc["name"]))

    def run(self, exec: str, mpiexec: str, force):
        if self._workspace.exists() and not force:
            raise FileExistsError("Workspace directory exists, please remove")
        self._workspace.mkdir(parents=True, exist_ok=True)

        if self._testcases_filtered is None:
            raise ValueError("initialize() needs to be called before running")

        # Resolve the absolute path of executable to test
        exec = pathlib.Path(exec).resolve()

        os.environ["OMP_NUM_THREADS"] = str(self._nthreads)
        args = []
        if mpiexec == "srun":
            args = [mpiexec, exec]
        elif mpiexec in ["mpirun", "mpiexec"]:
            args = [mpiexec, "-np", str(self._ntasks), exec]

        print()
        if not self._testcases_filtered:
            print("No test case to run")
            return
        for tc in self._testcases_filtered:
            dname = tc["directory"]
            src = self._dir_input / dname
            src_librpa = src / "librpa"
            dst = self._dir_testcase / dname
            dst_librpa = dst / "librpa"
            # prepare inputs
            dst_librpa.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src_librpa / "librpa.in", dst_librpa)
            with tarfile.open(src / "dataset.tar.gz", "r:gz") as tar:
                tar.extractall(path=dst)
            # prepare inputs
            print("Running {} [{}]".format(tc["name"], dname))
            run_librpa(args, dst_librpa)
        print("Finished test calculations")
        print()

    def analyze(self):
        status = 0
        good_all = []
        for tc in self._testcases_filtered:
            # name = tc["name"]
            dname = tc["directory"]
            test = self._dir_testcase / dname
            refr = self._dir_ref / dname
            results = []
            for v in tc["validates"]:
                entry = Validate(v["name"],
                                 v["file"],
                                 v["comparison"],
                                 v["headers"],
                                 v["rows"],
                                 v["regex"],
                                 v["occurences"],
                                 v["binary_extract"],
                                 )
                good, msg = entry.evaluate(refr, test)
                good_all.append(good)
                results.append([good, msg])
            tc["results"] = results
        if not all(good_all):
            return 1
        return status

    # TODO: make output work
    def print(self, output):
        for g, gtcs in self._groups.items():
            gtcs_active = [tc for tc in gtcs if _disable_message(tc) is None]
            if not gtcs_active:
                continue
            print()
            print("Test group: {}".format(g))
            for tc in gtcs_active:
                results = tc.get("results", None)
                print()
                if results:
                    s = "Validate results for {} [directory: {}]: {}"
                    good_all = PASS_FAIL.get(all(x[0] for x in results))
                    print(s.format(tc["name"], tc["directory"], good_all))
                    for v, e in zip(tc["validates"], results):
                        print("- {:4s}: {:s}, {:s}".format(PASS_FAIL[e[0]], v["name"].strip(), e[1]))
                else:
                    s = "No results to validate for {} [directory: {}]"
                    print(s.format(tc["name"], tc["directory"]))
