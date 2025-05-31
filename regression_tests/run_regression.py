#!/usr/bin/env python3
"""Driver for running regression test for LibRPA
"""
import sys
from backend.commandparser import get_parser
from backend.xmlparser import XMLParser
from backend.driver import Driver


if __name__ == '__main__':
    args = get_parser().parse_args()

    # Parse the test suite
    suite = XMLParser(args.xml)

    # Create the test driver
    driver = Driver(args.dir_input, args.dir_ref, args.workspace, suite.groups)

    # Initialize workspace and run the tests
    driver.initialize(args.ntasks, args.nthreads, args.use_libri, args.use_greenx_api)

    if args.mode in ["run", "full"]:
        driver.run(args.librpa_exec, args.mpiexec, args.force)

    status = 0
    if args.mode in ["analyze", "full"]:
        status = driver.analyze()
        driver.print(args.output)
    sys.exit(status)
