#!/usr/bin/env python

import argparse


def parse_command_line():
    ''' parse command line options '''
    parser = argparse.ArgumentParser(
        description="Lint QE input file"
    )
    parser.add_argument(
        'file_in', metavar='file_in', type=str, help="file to read in and linted"
    )
    parser.add_argument(
        '-i', '--inplace', action='store_true', help="if set then qe input file will be overwritten"
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_command_line()
    import pw2py as pw
    inp = pw.qeinp.from_file(args.file_in)
    if not args.inplace:
        print(inp, end='')
    else:
        inp.write_file(args.file_in)

