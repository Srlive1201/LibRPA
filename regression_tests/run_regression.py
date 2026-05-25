#!/usr/bin/env python3
"""Driver for running regression test for LibRPA
"""
import sys
from backend.commandparser import get_parser
from backend.xmlparser import XMLParser
from backend.driver import TestDriver


if __name__ == '__main__':
    parser = get_parser()
    args = parser.parse_args()

    if args.mode is None:
        parser.print_help()
        sys.exit(1)

    # Parse the test suite
    suite = XMLParser(args.xml)

    # Create the test driver
    driver = TestDriver(args.dir_input, args.dir_ref,
                        args.workspace, suite.groups)

    # Initialize workspace and run the tests
    try:
        driver.initialize(args.ntasks, args.nthreads, args.use_libri,
                          args.only, args.exclude)
    except ValueError as exc:
        parser.error(str(exc))

    status = 0

    if args.mode == "list":
        driver.list()
        sys.exit(status)

    if args.mode in ["run", "full"]:
        driver.run(args.librpa_exec, args.mpiexec, args.force)

    if args.mode in ["analyze", "full"]:
        status = driver.analyze()
        driver.print(args.output)
        print()
        if status == 0:
            print("All selected tests PASSED :)")
        else:
            print("Some tests FAILED, please check above details")

    sys.exit(status)
