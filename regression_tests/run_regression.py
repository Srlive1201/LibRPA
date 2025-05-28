#!/usr/bin/env python3
"""Driver for running regression test for LibRPA
"""
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from backend.xmlparser import XMLParser
from backend.driver import Driver


def _parser():
    p = ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
    p.add_argument("librpa_exec", help="Path to LibRPA executable for test")
    p.add_argument("--xml", type=str, default="testsuite.xml",
                   help="XML file containing test case configurations, default: testsuite.xml")
    p.add_argument("-d", dest="workspace", type=str, default="workspace",
                   help="Working directory to run regression test, default: workspace/")
    p.add_argument("-n", "--ntasks", dest="ntasks", type=int, default=1,
                   help="Number of MPI tasks to run the tests, default: 1")
    p.add_argument("--nthreads", type=int, default=1,
                   help="Number of OpenMP threads to run the tests, default: 1")
    p.add_argument("--mpiexec", type=str, default="mpirun", help="MPI command")
    p.add_argument("-f", "--force", dest="force", action="store_true", help="Force running")
    p.add_argument("--use-libri", action="store_true",
                   help="Flag that the test executable is built with LibRI")
    p.add_argument("--use-greenx-api", action="store_true",
                   help="Flag that the test executable is built with GreenX API")
    p.add_argument("--dir-input", type=str, default="testcases/", help="Directory that contain input data for testcases")
    p.add_argument("--dir-ref", type=str, default="refs/", help="Directory that contain reference results")
    return p


if __name__ == '__main__':
    args = _parser().parse_args()

    # Parse the test suite
    suite = XMLParser(args.xml)

    # Create the test driver
    driver = Driver(args.dir_input, args.dir_ref, args.workspace, suite.groups, args.force)

    # Initialize workspace and run the tests
    driver.initialize(args.ntasks, args.nthreads, args.use_libri, args.use_greenx_api)
    driver.run(args.librpa_exec, args.mpiexec)
