import pathlib
import shutil
import tarfile
from collections import OrderedDict


class Driver:

    def __init__(self, dir_input: str, dir_ref: str, workspace: str, groups: dict, force: bool = False):
        self._dir_input = pathlib.Path(dir_input)
        self._dir_ref = pathlib.Path(dir_ref)
        self._workspace = pathlib.Path(workspace)

        # Check
        if not self._dir_input.exists():
            raise FileNotFoundError("Input directory does not exist")
        if not self._dir_ref.exists():
            raise FileNotFoundError("Reference directory does not exist")
        if self._workspace.exists() and not force:
            raise FileExistsError("Workspace directory exists, please remove")
        else:
            self._workspace.mkdir(parents=True, exist_ok=True)

        self._groups = groups
        # Testcases that are qualified after initialized with build and runtime conditions
        self._testcases_run = None

    def initialize(self, ntasks: int, nthreads: int, use_libri: bool, use_greenx_api: bool):
        self._testcases_run = []

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
                name = tc["name"]
                directory = tc["directory"]
                build = tc["build"]

                if build["require_libri"] and not use_libri:
                    skip_due_to_build.append(tc)
                    continue
                if build["require_greenx_api"] and not use_greenx_api:
                    skip_due_to_build.append(tc)
                    continue

                run = tc["run"]
                if ntasks in run["ntasks_disable"]:
                    skip_due_to_run.append(tc)
                    continue

                if ntasks in run["nthreads_disable"]:
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

        # Create workspace directories and copy input data
        for tc in self._testcases_run:
            dname = tc["directory"]
            src = self._dir_input / dname
            dst = self._workspace / dname
            dst.mkdir()
            shutil.copy2(src / "librpa.in", dst)
            with tarfile.open(src / "input_librpa.tar.gz", "r:gz") as tar:
                tar.extractall(path=dst)

    def run(self, exec: str, mpiexec: str, ):
        if self._testcases_run is None:
            raise ValueError("initialize() needs to be called before running")
