#!/usr/bin/env python3

import argparse


def parse_command_line():
    ''' parse command line options '''
    parser = argparse.ArgumentParser(
        description='Copy geometry from a file (inpfile) and into a qe input file (qefile)'
    )
    parser.add_argument(
        'inpfile', metavar='inpfile', type=str, help='input geometry file to read in and then convert'
    )
    parser.add_argument(
        'qefile', metavar='qefile', type=str, help='qe input file to be write to'
    )
    parser.add_argument(
        '-o', '--only', metavar='arg', type=str, help='only copy (c)ell parameters | (a)tomic positions'
    )

    return parser.parse_args()


def parse_only_option(only):
    ''' parse argument passed to -o, --only (return argument to pass to load_geo) '''
    if (args.only == 'c') or (args.only == 'cell'):
        option = 'cell'
    elif (args.only == 'a') or (args.only == 'atoms'):
        option = 'atoms'
    elif args.only is None:
        option = 'default'
    else:
        raise ValueError('Unrecognized arg for -o, --only option: {}'.format(only))

    return option


if __name__ == "__main__":
    args = parse_command_line()
    only = parse_only_option(args.only)
    import pw2py as pw
    geo = pw.atomgeo.from_file(args.inpfile)
    inp = pw.qeinp.from_file(args.qefile)
    inp.load_geo(geo, load=only)
    print(inp)
