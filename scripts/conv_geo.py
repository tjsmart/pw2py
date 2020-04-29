#!/usr/bin/env python3

import argparse


def parse_command_line():
    parser = argparse.ArgumentParser(description='Convert a geometry file (filepath) to a different format')
    parser.add_argument(
        'ftype', metavar='format', type=str, help='format type of output, can be: [qeinp|vasp|xsf|xyz|jdftx]'
    )
    parser.add_argument(
        'inpfile', metavar='filepath', type=str, help='input geometry file to read in and then convert'
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_command_line()
    import pw2py as pw
    if args.ftype == 'qeinp':
        print(pw.atomgeo.from_file(args.inpfile))
    else:
        print(pw.atomgeo.from_file(args.inpfile).__str__(ftype=args.ftype))
    # print(pw.atomgeo.from_file(args.inpfile).__str__(ftype=args.ftype))
