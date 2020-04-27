from collections import Counter
import f90nml
from numpy import ones, array

from .. import atomgeo
from .._common.resource import _determine_ftype, _read_atomic_positions, _resolve_continuation_lines


@classmethod
def from_file(cls, filename, ftype='auto'):
    '''
    generate qegeo object from file
    '''
    if ftype == 'auto':
        ftype = _determine_ftype(filename)
    elif ftype not in ['vasp', 'xyz', 'xsf', 'jdftx', 'qeinp', 'qeout']:
        raise ValueError('Value of ftype not recognized')

    # load atomgeo object from file
    geo = atomgeo.from_file(filename, ftype=ftype)

    # determine other parameters (i.e. if_pos, ibrav, A, B, C, celldm)
    # defaults
    nml = f90nml.Namelist({'ibrav': 0})
    if_pos = ones((geo.nat, 3), dtype=int)

    if ftype == 'jdftx':
        # reread lines where all continuation lines have been resolved (i.e. lines that end in '\')
        lines = iter(_resolve_continuation_lines(filename))
        if_pos = []
        for line in lines:
            if line.startswith('ion '):
                line = line.split('#')[0].strip()
                if_pos.append(array(line.split()[5] * 3, dtype=int))

    elif ftype == 'qeout':
        with open(filename) as f:
            lines = iter(f.readlines())
        # then file is output
        for line in lines:
            if "ATOMIC_POSITIONS" in line:
                _, _, _, if_pos = _read_atomic_positions(lines, geo.nat, no_if_pos=False)

    elif ftype == 'qeinp':
        with open(filename) as f:
            lines = iter(f.readlines())

        nml = f90nml.read(filename)
        # filter namelist
        nml = nml['system']
        keys = list(nml.keys())
        for key in keys:
            if key not in ['ibrav', 'a', 'b', 'c', 'celldm', 'cosab', 'cosac', 'cosbc',
                           'nat', 'ntyp']:
                nml.pop(key)

        # ibrav = nml['system']['ibrav']
        # if ibrav != 0:
        #     if 'celldm' in nml['system']:
        #         celldm = nml['system']['celldm']
        #     elif 'A' in nml['system']:
        #         A = nml['system']['A']
        #         if 'B' in nml['system']:
        #             B = nml['system']['B']
        #         elif 'C' in nml['system']:
        #             C = nml['system']['C']
        #         elif 'cosAB' in nml['system']:
        #             cosAB = nml['system']['cosAB']
        #         elif 'cosAC' in nml['system']:
        #             cosAC = nml['system']['cosAC']
        #         elif 'cosBC' in nml['system']:
        #             cosBC = nml['system']['cosBC']

        #     else:
        #         raise ValueError('ibrav !=0 but celldm or A was not set in input file')

        for line in lines:
            if 'ATOMIC_POSITIONS' in line:
                line = line.split('!')[0].split('#')[0]   # trim away comments
                _, _, _, if_pos = _read_atomic_positions(lines, geo.nat, no_if_pos=False)

    if 'nat' not in nml:
        nml['nat'] = geo.nat
    if 'ntyp' not in nml:
        nml['ntyp'] = len(Counter(geo.ion))

    # return cls(geo=geo, if_pos=if_pos, ibrav=ibrav, A=A, B=B, C=C, cosAB=cosAB, cosAC=cosAC, \
    # cosBC=cosBC, celldm=celldm)
    return cls(geo=geo, if_pos=if_pos, nml=nml)


def write_file(self, filename):
    '''
    write qegeo to file (QE format)
    '''
    with open(filename, 'w') as f:
        f.write(str(self))
