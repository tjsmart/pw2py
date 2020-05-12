#!/usr/bin/env python3

import argparse
import numpy as np
import os
from copy import deepcopy


default_ratios = np.linspace(0, 1, 11)
default_nonrad_ratios = np.array([0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.5, 0.75, 0.8, 0.85, 0.9, 0.95, 1])


def parse_command_line():
    ''' parse command line options '''
    parser = argparse.ArgumentParser(
        description="Linear mix geometry between two files. Creates subfolders ratio-* and stores files within."
    )
    parser.add_argument(
        'file1', metavar='file1', type=str, help="first qe input or atomic geometry file to read in"
    )
    parser.add_argument(
        'file2', metavar='file2', type=str, help="second qe input or atomic geometry file to read in"
    )
    parser.add_argument(
        '-r', '--ratios', metavar='arg', type=str, help="list of ratio for mixing (pass as space separated str)",
        default=default_ratios
    )
    parser.add_argument(
        '-g', '--geo', action='store_true', help="treat input file1 and file2 as atomgeo objects (default is qeinp)"
    )
    parser.add_argument(
        '-n', '--nonrad', action='store_true', help="use default nonrad list of ratios"
    )

    return parser.parse_args()


def mkdir_ignore(path):
    '''
    call os.mkdir but ignore FileExistsError
    '''
    try:
        os.mkdir(path)
    except FileExistsError:
        pass


if __name__ == "__main__":
    args = parse_command_line()
    # parse ratios
    if args.nonrad:
        ratios = default_nonrad_ratios
    if not np.array_equiv(args.ratios, default_ratios):
        ratios = np.fromstring(args.ratio, sep=' ', dtype=np.float64)
    else:
        ratios = default_ratios
    # comments to user
    print('Generating mix structures from: {} and {}'.format(args.file1, args.file2))
    print('Mixing ratios: ' + ' '.join('{:6.4f}'.format(r) for r in ratios))

    import pw2py as pw

    # determine basename of files to generate
    if not args.geo:
        basename = 'scf.in'
    else:
        ftype = pw._common.resource._determine_ftype(args.file1)
        if ftype == 'vasp':
            basename = 'POSCAR'
        elif ftype in ['xyz', 'xsf']:
            basename = 'mix.' + ftype
        elif ftype == 'jdftx':
            basename = 'mix.pos'
        else:
            basename = 'scf.in'

    # make directories for putting ratio-* folders
    lin_dir1 = 'lin1'
    lin_dir2 = 'lin2'
    mkdir_ignore(lin_dir1)
    mkdir_ignore(lin_dir2)

    if args.geo:
        obj1 = pw.atomgeo.from_file(args.file1)
        obj2 = pw.atomgeo.from_file(args.file2)
    else:
        obj1 = pw.qeinp.from_file(args.file1)
        obj2 = pw.qeinp.from_file(args.file2)

    mix1 = deepcopy(obj1)
    mix2 = deepcopy(obj2)

    for ratio in ratios:
        print('.', end='')
        # make ratio dir
        dir1 = os.path.join(lin_dir1, 'ratio-{:6.4f}'.format(ratio))
        dir2 = os.path.join(lin_dir2, 'ratio-{:6.4f}'.format(ratio))
        mkdir_ignore(dir1)
        mkdir_ignore(dir2)
        # mix positions
        mix1.pos = (1 - ratio) * obj1.pos + ratio * obj2.pos
        mix2.pos = (1 - ratio) * obj1.pos + ratio * obj2.pos
        # write files
        path1 = os.path.join(dir1, basename)
        path2 = os.path.join(dir2, basename)
        mix1.write_file(path1)
        mix2.write_file(path2)

    print('\nDone! :)')
