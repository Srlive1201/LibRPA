#!/usr/bin/env python3
"""Driver for running regression test for LibRPA
"""
from argparse import ArgumentParser, RawDescriptionHelpFormatter
from backend.xmlparser import XMLParser


def _parser():
    p = ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
    p.add_argument("--xml", type=str, default="testsuite.xml",
                   help="XML file containing test case configurations, default: testsuite.xml")
    p.add_argument("-d", dest="workspace", type=str, default="workspace",
                   help="Working directory to run regression test, default: workspace/")
    p.add_argument("-n", "--ntasks", dest="ntasks", type=int, default=1,
                   help="Number of MPI tasks to run the tests, default: 1")
    p.add_argument("--nthreads", type=int, default=1,
                   help="Number of OpenMP threads to run the tests, default: 1")
    return p


if __name__ == '__main__':
    args = _parser().parse_args()
    xml = XMLParser(args.xml)
