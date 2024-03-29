import numpy as np
import os
from math import pi

from .constants import bohr_to_angstrom


def _read_qe_card_option(line, card_name):
    ''' strip away everything to read option from qe card '''
    option = line.split('!')[0].split('#')[0].split(card_name)[1]
    option = ''.join(filter(str.isalpha, option)).lower()
    return option


def _resolve_continuation_lines(filename):
    '''
    return lines of file without continuation lines (used for jdftx)
    '''
    lines = []
    with open(filename) as f:
        for line in f:
            line = line.rstrip('\n')
            while line.endswith('\\'):
                line = line[:-1] + next(f).rstrip('\n')
            lines.append(line)
    return lines


def _read_atomic_positions(lines, nat, no_if_pos=True):
    '''
    iterate through lines of atomic positons returning ion, pos, if_pos (QE format)
    '''
    pos, ion, if_pos = [], [], []
    for _ in range(nat):
        nl = next(lines).split('!')[0].split()
        ion.append(nl[0])
        pos.append(np.array(nl[1:4], dtype=np.float64))
        if not no_if_pos:
            try:
                # explicit so indexError will be thrown
                if_pos.append(np.array([nl[4], nl[5], nl[6]], dtype=np.int))
            except IndexError:
                if_pos.append(np.array([1, 1, 1], dtype=np.int))
    if no_if_pos:
        return lines, np.array(pos), ion
    else:
        return lines, np.array(pos), ion, np.array(if_pos)


def _ibrav_to_par(system_nml, units="angstrom"):
    '''
    calculate cell parameters from bravais lattice using qe system namelist as input

    note the default cell parameters are in angstrom
    '''
    # check value of ibrav
    if int(system_nml['ibrav']) == 0:
        return None
    elif int(system_nml['ibrav']) == 1:
        # v1 = a(1,0,0),  v2 = a(0,1,0),  v3 = a(0,0,1)
        par = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64)
    elif int(system_nml['ibrav']) == 2:
        # v1 = (a/2)(-1,0,1),  v2 = (a/2)(0,1,1), v3 = (a/2)(-1,1,0)
        par = np.array([[-1, 0, 1], [0, 1, 1], [-1, 1, 0]],
                       dtype=np.float64) * 0.5
    elif int(system_nml['ibrav']) == 3:
        # v1 = (a/2)(1,1,1),  v2 = (a/2)(-1,1,1),  v3 = (a/2)(-1,-1,1)
        par = np.array([[1, 1, 1], [-1, 1, 1], [-1, -1, 1]],
                       dtype=np.float64) * 0.5
    elif int(system_nml['ibrav']) == 4:
        # v1 = a(1,0,0),  v2 = a(-1/2,sqrt(3)/2,0),  v3 = a(0,0,c/a)
        par = np.array([[1, 0, 0], [-1/2, np.sqrt(3)/2, 0],
                        [0, 0, 1]], dtype=np.float64)
    elif int(system_nml['ibrav']) == 6:
        # v1 = a(1,0,0),  v2 = a(0,1,0),  v3 = a(0,0,c/a)
        par = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64)
    elif int(system_nml['ibrav']) == 8:
        # v1 = (a,0,0),  v2 = (0,b,0), v3 = (0,0,c)
        par = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=np.float64)
    elif int(system_nml['ibrav']) == -12:
        # v1 = (a,0,0),  v2 = (0,b,0), v3 = (c*cos(beta),0,c*sin(beta))
        par = np.array([[1, 0, 0], [0, 1, 0], [1, 0, 1]], dtype=np.float64)
    else:
        raise ValueError(
            "Unsupported value for ibrav '{}'".format(system_nml['ibrav']))

    # rescale lattice parameters
    if 'a' in system_nml.keys():
        par *= np.float64(system_nml['a'])
        if int(system_nml['ibrav']) in [8, -12]:
            par[1, 1] = np.float64(system_nml['b'])
        if int(system_nml['ibrav']) in [4, 6, 8]:
            par[2, 2] = np.float64(system_nml['c'])
        if int(system_nml['ibrav']) in [-12]:
            sinac = 1 - (system_nml['cosac'])**2
            par[2] = np.float64(system_nml['c']) * np.array([system_nml['cosac'], 0, sinac])
    elif 'celldm' in system_nml.keys():
        par *= np.float64(system_nml['celldm'][0]) * bohr_to_angstrom
        if int(system_nml['ibrav']) in [8, -12]:
            par[1, 1] *= np.float64(system_nml['celldm'][1])
        if int(system_nml['ibrav']) in [4, 6, 8]:
            par[2, 2] *= np.float64(system_nml['celldm'][2])
        if int(system_nml['ibrav']) in [-12]:
            par[2] *= np.float64(system_nml['celldm'][2])
            par[2, 0] *= system_nml['celldm'][4]
            sinac = 1 - (system_nml['celldm'][4])**2
            par[2, 2] *= sinac

    # check units
    if units == "bohr":
        par /= bohr_to_angstrom

    return par


def _convert_par(par, in_units, out_units="angstrom", alat=None, alat_units="angstrom"):
    '''
    Convert cell parameters from 'in_units' to 'out_units'
        ['alat', 'angstrom', 'bohr']
    '''
    # convert alat
    if alat_units == "bohr":
        alat *= alat * bohr_to_angstrom
    elif alat_units != "angstrom":
        raise ValueError("Invalid value of alat_units: '{}".format(alat_units))
    # define unit dictionary for unit conversion
    conversion = {'alat': alat, 'angstrom': 1, 'bohr': bohr_to_angstrom}
    # check input
    if in_units == out_units:
        return par
    elif in_units not in conversion:
        raise ValueError("Invalid value of in_units: '{}'".format(in_units))
    elif out_units not in conversion:
        raise ValueError("Invalid value of out_units: '{}'".format(out_units))
    elif (in_units == "alat" or out_units == "alat") and alat is None:
        raise ValueError(
            "Requested alat conversion but value for alat was not provided.")
    # calculate prefactor
    return par / np.float64(conversion[out_units.lower()]) * np.float64(conversion[in_units.lower()])


def _convert_pos(pos, in_units, out_units="angstrom", alat=None, alat_units="angstrom", par=None, par_units=None):
    '''
    Convert atomic positions from 'in_units' to 'out_units'
        ['alat', 'bohr', 'angstrom', 'crystal']
    '''
    if in_units == out_units:
        return pos
    elif in_units != 'crystal' and out_units != 'crystal':
        return _convert_par(pos, in_units=in_units, out_units=out_units, alat=alat, alat_units=alat_units)
    elif par is None or par_units is None:
        raise ValueError(
            "Requested crystal conversion but cell parameters (par) were not provided.")

    # cases with crystal conversion
    prefactor = 1.0
    if (par_units == 'angstrom' and out_units == 'bohr') or (par_units == 'bohr' and in_units == 'angstrom'):
        prefactor /= bohr_to_angstrom
    elif (par_units == 'bohr' and out_units == 'angstrom') or (par_units == 'angstrom' and in_units == 'bohr'):
        prefactor *= bohr_to_angstrom

    if in_units == 'crystal':
        return pos.dot(par) * prefactor
    elif out_units == 'crystal':
        return pos.dot(np.linalg.inv(par)) * prefactor


def _determine_ftype(filename):
    '''
    determine the type of file based on the filename extension
    (also handles a few other cases)

    input
    ---
        filename (str)
            - path to file

    returns
    ---
        ftype (str)
            - determined type of file, can be: ['vasp', 'xyz', 'xsf', 'jdftx', 'qeinp', 'qeout']
    '''
    base, ext = os.path.splitext(filename.lower())

    # cases without extensions
    if not ext:
        filename = base.split('/')[-1]
        if filename in ['poscar', 'contcar']:
            return 'vasp'
        elif os.path.exists(base + '.in') and os.path.exists(base + '.out'):
            # filename is a prefix for 'prefix.in' and 'prefix.out'
            return 'qeout'

    # recognized extensions
    ftype_extensions = {
        'vasp': ['.poscar', '.contcar', '.vasp'],
        'xyz': ['.xyz'],
        'xsf': ['.xsf'],
        'jdftx': ['.pos', '.jdftx'],
        'qeinp': ['.in'],
        'qeout': ['.out']
    }
    for ftype in ftype_extensions:
        # if ext is recognized return corresponding ftype
        if ext in ftype_extensions[ftype]:
            return ftype

    # raise ValueError when determining ftype fails
    raise ValueError(
        'Unrecognizable extension or unable to determine the type of file for: {}'.format(filename))


def savetxt_with_maxcol(f, X, columns=6, fmt='%13.6e'):
    '''
    calls np.savetxt, fits X.reshape((X.size/columns, columns)) and saves to f

    if X.size%columns != 0, the remaining elements are handled appropriately appended to the end
    '''
    X = np.array(X).flatten()
    rows, remainder = divmod(X.size, columns)
    np.savetxt(f, X[:rows*columns].reshape((rows, columns)), fmt=fmt)
    if remainder != 0:
        np.savetxt(f, X[-remainder:].reshape((1, remainder)), fmt=fmt)


def _calc_rec(par: np.ndarray) -> np.ndarray:
    '''
    Calculate reciprocal lattice vectors from par
    (units of 2pi / alat)
    '''
    alat = np.linalg.norm(par[0])
    V = par[0].dot(np.cross(par[1], par[2]))
    return np.cross(par[[1, 2, 0]], par[[2, 0, 1]]) / V * alat
