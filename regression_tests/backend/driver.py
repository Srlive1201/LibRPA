import os
import pathlib
import shutil
import tarfile
from collections import OrderedDict
from typing import Tuple

from .utils import run_librpa


__all__ = ["Driver"]


def _check_build_run_filter(tc: dict, use_libri: bool, use_greenx_api: bool) -> Tuple[bool, bool]:
    # Filter build options
    build = tc["build"]
    if build["require_libri"] and not use_libri:
        return True, False
    if build["require_greenx_api"] and not use_greenx_api:
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


class Driver:

    def __init__(self, dir_input: str, dir_ref: str, workspace: str, groups: dict):
        self._dir_input = pathlib.Path(dir_input)
        self._dir_ref = pathlib.Path(dir_ref)
        self._workspace = pathlib.Path(workspace)

        # Check
        if not self._dir_input.exists():
            raise FileNotFoundError("Input directory does not exist")
        if not self._dir_ref.exists():
            raise FileNotFoundError("Reference directory does not exist")

        self._groups = groups
        self._ntasks = None
        self._nthreads = None
        # Testcases that are qualified after initialized with build and runtime conditions
        self._testcases_run = None

    def initialize(self, ntasks: int, nthreads: int, use_libri: bool,
                   use_greenx_api: bool):

        self._testcases_run = []
        self._ntasks = ntasks
        self._nthreads = nthreads

        print("Initializing workspace targeting:")
        print("- MPI tasks :", ntasks)
        print("- threads   :", nthreads)
        print("- LibRI     ?", use_libri)
        print("- GreenX API?", use_greenx_api)

        # Filter test cases
        skip_due_to_build = []
        skip_due_to_run = []

        for g, gtcs in self._groups.items():
            for tc in gtcs:
                filter_build, filter_run = _check_build_run_filter(tc, use_libri, use_greenx_api)
                if filter_build:
                    skip_due_to_build.append(tc)
                    continue
                if filter_run:
                    skip_due_to_run.append(tc)
                    continue
                self._testcases_run.append(tc)

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

        if self._testcases_run is None:
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
        if not self._testcases_run:
            print("No test case to run")
            return
        for tc in self._testcases_run:
            dname = tc["directory"]
            # prepare inputs
            src = self._dir_input / dname
            dst = self._workspace / dname
            dst.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src / "librpa.in", dst)
            with tarfile.open(src / "input_librpa.tar.gz", "r:gz") as tar:
                tar.extractall(path=dst)
            d = self._workspace / dname
            print("Running {} [{}]".format(tc["name"], dname))
            run_librpa(args, d)

    def get_testcases_run(self):
        return self._testcases_run
