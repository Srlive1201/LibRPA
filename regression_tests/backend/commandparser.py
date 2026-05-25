from argparse import ArgumentParser, RawDescriptionHelpFormatter


def get_parser():
    p = ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)

    subp = p.add_subparsers(dest="mode")
    p_full = subp.add_parser("full", help="Run the regression tests and perform diff analysise")
    p_run = subp.add_parser("run", help="Run the regression tests")
    p_list = subp.add_parser("list", help="List test cases")
    p_analyze = subp.add_parser("analyze", help="Perform diff analysis")

    _add_testsuite_options(p_full)
    _add_testsuite_options(p_run)
    _add_testsuite_options(p_list)
    _add_testsuite_options(p_analyze)

    _add_execute_options(p_full)
    _add_execute_options(p_run)

    _add_analyze_options(p_full)
    _add_analyze_options(p_analyze)

    return p


def _add_testsuite_options(p: ArgumentParser):
    p.add_argument("--xml", type=str, default="testsuite.xml",
                   help="XML file containing test case configurations, default: testsuite.xml")
    p.add_argument("-d", dest="workspace", type=str, default="workspace",
                   help="Working directory to run regression test, default: workspace/")
    p.add_argument("--use-libri", action="store_true",
                   help="Flag that the test executable is built with LibRI")
    p.add_argument("-n", "--ntasks", dest="ntasks", type=int, default=1,
                   help="Number of MPI tasks to run the tests, default: 1")
    p.add_argument("--nthreads", type=int, default=1,
                   help="Number of OpenMP threads to run the tests, default: 1")
    p.add_argument("--dir-input", type=str, default="testcases/", help="Directory that contain input data for testcases")
    p.add_argument("--dir-ref", type=str, default="refs/", help="Directory that contain reference results")
    scope = p.add_mutually_exclusive_group()
    scope.add_argument("--only", nargs="+", metavar="TESTCASE",
                       help="Only enable test cases with these directory names")
    scope.add_argument("--exclude", nargs="+", metavar="TESTCASE",
                       help="Disable test cases with these directory names")


def _add_execute_options(p: ArgumentParser):
    p.add_argument("librpa_exec", help="Path to LibRPA executable for test")
    p.add_argument("-f", "--force", dest="force", action="store_true", help="Force running")
    p.add_argument("--mpiexec", type=str, default="mpirun", help="MPI command")


def _add_analyze_options(p: ArgumentParser):
    p.add_argument("-o", "--output", default="regression.log", help="output file of regression test")
