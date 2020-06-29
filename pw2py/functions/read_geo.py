import f90nml
import numpy as np
import os

from .. import element
from .._common.constants import bohr_to_angstrom
from .._common.resource import _ibrav_to_par, _read_qe_card_option, _determine_ftype


def _read_qe_atomic_positions(f, nat, read_if_pos=False):
    ''' read ion, pos, if_pos from qe atomic_positions '''
    ion, pos, if_pos = [], [], []
    for _ in range(nat):
        nl = f.readline().split('!')[0].split('#')[0].split()
        ion.append(nl[0])
        pos.append(nl[1:4])
        if read_if_pos:
            try:
                # cannot use [4:7], because IndexError will not be thrown
                if_pos.append([nl[4], nl[5], nl[6]])
            except IndexError:
                if_pos.append([1, 1, 1])  # default values
    ion = np.array(ion)
    pos = np.array(pos, dtype=np.float64)
    if read_if_pos:
        if_pos = np.array(if_pos, dtype=np.int32)

    return ion, pos, if_pos


def read_vasp(filename):
    '''
    read geometry from vasp file

    returns ion(array), par(array), par_units(str), pos(array), pos_units(str)
    '''
    with open(filename) as f:
        # skip first line
        f.readline()
        # read cell par
        alat = np.float64(f.readline())
        par = np.array([f.readline().split()[0:3] for _ in range(3)], dtype=np.float64) * alat
        par_units = "angstrom"
        # read ions
        species = f.readline().split()
        count = np.fromstring(f.readline(), sep=' ', dtype=int)
        ion = []
        for s, c in zip(species, count):
            ion += [s] * c
        ion = np.array(ion)
        nat = len(ion)
        # read pos type
        pos_units = "crystal" if f.readline().strip().lower() == "direct" else "angstrom"
        pos = np.array([f.readline().split()[0:3] for _ in range(nat)], dtype=np.float64)

    return ion, par, par_units, pos, pos_units


def read_xyz(filename, a=20.0, par_units='angstrom'):
    '''
    read geometry from xyz file

    returns ion(array), par(array), par_units(str), pos(array), pos_units(str)
    '''
    par_units = par_units.lower()
    if par_units not in ['angstrom', 'bohr']:
        raise ValueError('Invalid value for par_units: {}'.format(par_units))
    with open(filename) as f:
        # then file is xyz format
        par_units = 'angstrom'
        # set par to a box with dimension of 'a' (units set by par_units)
        par = a * np.eye(3)
        # first line is number of atoms
        nat = int(f.readline())
        # next line is a comment
        f.readline()
        pos_units = 'angstrom'
        # next nat lines are positions
        ion, pos = [], []
        for _ in range(nat):
            nl = f.readline().split()[0:4]
            ion.append(nl[0])
            pos.append(nl[1:4])
    ion = np.array(ion)
    pos = np.array(pos, dtype=np.float64)

    return ion, par, par_units, pos, pos_units


def read_xsf(filename):
    '''
    read geometry from xsf file

    returns ion(array), par(array), par_units(str), pos(array), pos_units(str)
    '''
    par_units = 'angstrom'
    pos_units = 'angstrom'
    with open(filename) as f:
        for line in f:
            if line.startswith('PRIMVEC') or ' PRIMVEC' in line:
                # if statement above avoids reading 'RECIP-PRIMVEC'
                # read cell parameters
                par = np.array([f.readline().split()[0:3] for _ in range(3)], dtype=np.float64)
            elif 'PRIMCOORD' in line:
                # read ions and pos
                nat = int(f.readline().split()[0])
                ion, pos = [], []
                for _ in range(nat):
                    nl = f.readline().split()
                    try:
                        ion.append(
                            element.request(int(nl[0]), 'symbol')
                        )
                    except ValueError:
                        ion.append(nl[0])
                    pos.append(nl[1:4])
                ion = np.array(ion)
                pos = np.array(pos, dtype=np.float64)
                break

    return ion, par, par_units, pos, pos_units


def read_jdftx(filename):
    '''
    read geometry from jdftx file

    returns ion(array), par(array), par_units(str), pos(array), pos_units(str)
    '''
    def resolve_continuation_lines(line):
        while line.endswith('\\'):
            line = line[:-1] + f.readline().rstrip('\n')
        return line

    with open(filename) as f:
        par_units = 'bohr'
        pos_units = 'crystal'
        ion, pos = [], []
        for line in f:
            line = line.strip()
            if line.startswith('lattice '):
                # read par
                lines = resolve_continuation_lines(line)
                lines = lines.strip('lattice').split('#')[0]
                par = np.fromstring(lines, dtype=np.float64, sep=' ').reshape((3, 3)).T
            elif line.startswith('ion '):
                # read ion, pos
                ion.append(line.split()[1])
                pos.append(line.split()[2:5])
            elif line.startswith('coords-type '):
                # read pos_units
                pos_units = 'crystal' if line.split()[1] == 'lattice' else 'bohr'
    ion = np.array(ion)
    pos = np.array(pos, dtype=np.float64)

    return ion, par, par_units, pos, pos_units


def read_qeout(filename, par_units='angstrom', read_if_pos=False):
    '''
    read geometry from qeout file

    filename can be a prefix (i.e. read 'filename.out' and 'filename.in')

    par_units may be overwritten if par_units are read from qeinp
    code tries to read par from qeinp since precision of par qeout is low if not vc-relax

    returns ion(array), par(array), par_units(str), pos(array), pos_units(str)
    '''
    if not os.path.exists(filename) and os.path.exists(filename + '.out'):
        filename += '.out'

    par_units = par_units.lower()
    if par_units not in ['angstrom', 'bohr']:
        raise ValueError('Invalid value for par_units: {}'.format(par_units))

    # try to find input file to read cell parameters
    par = None
    prefix, extension = os.path.splitext(filename)
    if extension == '.out' and os.path.exists(prefix + '.in'):
        nml = f90nml.read(prefix + '.in')
        if ('calculation' not in nml['control']) or (nml['control']['calculation'] in ['scf', 'nscf', 'bands']):
            # no relaxations preformed in output file just read qeinp file
            return read_qeinp(prefix + '.in', read_if_pos=read_if_pos)
        elif ('calculation' in nml['control']) and (nml['control']['calculation'] in ['vc-relax', 'vc-md']):
            # don't read input file, want cell parameters from output
            pass
        elif nml['system']['ibrav'] == 0:
            with open(prefix + '.in') as f:
                for line in f:
                    if 'CELL_PARAMETERS' in line:
                        par_units = _read_qe_card_option(line, 'CELL_PARAMETERS')
                        par = np.array([f.readline().split()[0:3] for _ in range(3)], dtype=np.float64)
                        break
        else:
            # need to convert ibrav to par
            par = _ibrav_to_par(nml['system'], units='angstrom')
            par_units = 'angstrom'

    # now read actual provided qeout file
    with open(filename) as f:
        for line in f:
            if 'celldm(1)=' in line:
                if par is not None:
                    continue
                alat = np.float64(line.split()[1])
                if par_units == 'angstrom':
                    alat *= bohr_to_angstrom
            elif 'number of atoms/cell' in line:
                nat = int(line.split()[-1])
            elif 'crystal axes:' in line:
                if par is not None:
                    continue
                par = np.array([f.readline().split()[3:6] for _ in range(3)], dtype=np.float64) * alat
            elif 'CELL_PARAMETERS' in line:
                alat = np.float64(line.strip().split('alat=')[1].rstrip(')'))
                if par_units == 'angstrom':
                    alat *= bohr_to_angstrom
                par = np.array([f.readline().split()[:3] for _ in range(3)], dtype=np.float64) * alat
            elif 'ATOMIC_POSITIONS' in line:
                pos_units = _read_qe_card_option(line, 'ATOMIC_POSITIONS')
                ion, pos, if_pos = _read_qe_atomic_positions(f, nat, read_if_pos=read_if_pos)

    if read_if_pos:
        return ion, par, par_units, pos, pos_units, if_pos
    else:
        return ion, par, par_units, pos, pos_units


def read_qeinp(filename, read_if_pos=False):
    '''
    read geometry from qeinp file

    returns ion(array), par(array), par_units(str), pos(array), pos_units(str)
    '''
    alat = None
    # read files namelist
    nml = f90nml.read(filename)
    if 'a' in nml['system']:
        alat = nml['system']['a']
    elif 'celldm' in nml['system']:
        alat = nml['system']['celldm'][0] * bohr_to_angstrom
    nat = int(nml['system']['nat'])
    # now iterate through file
    with open(filename) as f:
        for line in f:
            if 'CELL_PARAMETERS' in line:
                par_units = _read_qe_card_option(line, 'CELL_PARAMETERS')
                par = np.array([f.readline().split()[0:3] for _ in range(3)], dtype=np.float64)
            elif 'ATOMIC_POSITIONS' in line:
                pos_units = _read_qe_card_option(line, 'ATOMIC_POSITIONS')
                ion, pos, if_pos = _read_qe_atomic_positions(f, nat, read_if_pos=read_if_pos)
        if int(nml['system']['ibrav']) != 0:
            # if lattice is specified by ibrav then build par from ibrav
            # note this routine handles alat on its own
            par = _ibrav_to_par(nml['system'], units="angstrom")
            par_units = "angstrom"
        elif alat is not None:
            par *= alat

    if read_if_pos:
        return ion, par, par_units, pos, pos_units, if_pos
    else:
        return ion, par, par_units, pos, pos_units


def read_geo(filename, ftype='auto', read_if_pos=False):
    '''
    read geometry of an arbitrary geometry file (filename)

    input
    ---
        filename (str)
            - path of file to be read
        ftype (str) (optional)
            - type of file to be read, can be: ['vasp', 'xyz', 'xsf', 'jdftx', 'qeinp', 'qeout']
            - default is 'auto', wherein ftype is determined based on extension of file or other
        read_if_pos (bool) (optional)
            - only used if ftype='qeinp', in which case if_pos will be read and returned

    returns
    ---
        tuple
            - ion: array of ion names, shape = (nat,)
            - par: array of cell parameters, shape = (3, 3)
            - par_units: str, can be 'angstrom', or 'bohr'
            - pos: array of atomic positons, ndim = 2, shape (nat, 3)
            - pos_units: str, can be 'angstrom', 'crystal', or 'bohr'
    '''
    # check provided file type (ftype)
    valid_ftypes = ['auto', 'vasp', 'qeinp', 'qeout', 'jdftx', 'xsf', 'xyz']
    if ftype.lower() not in valid_ftypes:
        raise ValueError("ftype value not recognized or supported '{}'".format(ftype))
    elif ftype.lower() == 'auto':
        # determine file type based on extension or other
        ftype = _determine_ftype(filename)

    # call appropriate function
    if ftype == 'qeinp':
        return read_qeinp(filename, read_if_pos=read_if_pos)
    elif ftype == 'qeout':
        return read_qeout(filename, read_if_pos=read_if_pos)
    elif ftype == 'jdftx':
        return read_jdftx(filename)
    elif ftype == 'xyz':
        return read_xyz(filename)
    elif ftype == 'xsf':
        return read_xsf(filename)
    elif ftype == 'vasp':
        return read_vasp(filename)
