#!/usr/bin/env python3

import argparse


def parse_command_line():
    parser = argparse.ArgumentParser(
        description='Copy geometry from a file (inpfile) and into a qe input file (qefile)'
    )
    parser.add_argument(
        'inpfile', metavar='inpfile', type=str, help='input geometry file to read in and then convert'
    )
    parser.add_argument(
        'qefile', metavar='qefile', type=str, help='qe input file to be write to'
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_command_line()
    import pw2py as pw
    geo = pw.atomgeo.from_file(args.inpfile)
    inp = pw.qeinp.from_file(args.qefile)
    print(inp.load_geo(geo))
