#!/usr/bin/env python

import argparse
import pw2py as pw
import numpy as np
import warnings


def parse_command_line() -> argparse.ArgumentParser:
    ''' parse command line options '''
    parser = argparse.ArgumentParser(
        description="Read in two files and reorder atoms so that they match (according to distance)"
    )
    parser.add_argument(
        'file1', metavar='file1', type=str, help="first qe input or atomic geometry file to read in"
    )
    parser.add_argument(
        'file2', metavar='file2', type=str, help="second qe input or atomic geometry file to read in"
    )
    parser.add_argument(
        '-i', '--inplace', action='store_true', help="overwrite the contents of file2"
    )
    parser.add_argument(
        '-f', '--ftype', metavar='arg', type=str, help='format type of output, can be: [qeinp|vasp|xsf|xyz|jdftx]',
        default=None
    )

    return parser.parse_args()


def find_nearest(pos: np.ndarray, geo: pw.atomgeo) -> int:
    '''
    calculated distance of pos w.r.t. geo.pos, return closest
    '''
    gdistance = geo.calc_distance(pos)
    return np.argmin(gdistance)


if __name__ == '__main__':
    args = parse_command_line()

    # check for ftype
    ftype = pw._determine_ftype(args.file2) if args.ftype is None \
        else args.ftype.lower()

    # read in geo
    geo1 = pw.atomgeo.from_file(args.file1)
    geo2 = pw.atomgeo.from_file(args.file2)

    # read in copy of file2 to overwrite
    fix2 = pw.qeinp.from_file(args.file2) if ftype == 'qeinp' \
        else pw.atomgeo.from_file(args.file2)

    # convert atomic positions to angstrom (so distances are true)
    for obj in [geo1, geo2, fix2]:
        obj.pos_units = 'angstrom'

    # swap positions based on nearest neighbor index
    for i, pos in enumerate(geo1.pos):
        corrected_i = find_nearest(pos, geo2)
        fix2.pos[i] = geo2.pos[corrected_i]

    # convert back to crystal
    fix2.pos_units = 'crystal'

    if args.inplace:
        fix2.write_file(args.file2)
    else:
        if ftype == 'qeinp':
            print(fix2)
        else:
            print(fix2.to_string(ftype=ftype))
