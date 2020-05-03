import f90nml
from mendeleev import element
from numpy import array, fromstring, float64, eye
from warnings import warn

from .._common.constants import bohr_to_angstrom, nonalpha
from .._common.resource import _ibrav_to_par, _read_atomic_positions, \
    _resolve_continuation_lines


# Warning reading from xyz not supported!
_valid_ftypes = ['qe', 'jdftx', 'xyz', 'xsf', 'vasp']


@classmethod
def from_file(cls, filename, ftype='auto', xyz_par=20):
    '''
    generate atomgeo object from file

    input
    ----
    filename (str)
        - name/path to file to be read

    ftype (str)
        - specify filetype otherwise ftype='auto' will use filename to detect file type
        - acceptable values: ['auto', 'vasp', 'qeinp', 'qeout', 'jdftx', 'xsf', 'xyz']

    xyz_par (float)
        - optional, only used when reading from xyz, par is set to a box of with dimension xyz_par

    output
    ----
    atomgeo (atomgeo object)
    '''
    if ftype.lower() not in ['auto', 'vasp', 'qeinp', 'qeout', 'jdftx', 'xsf', 'xyz']:
        raise ValueError("ftype value not recognized '{}'".format(ftype))

    # store lines of file to iterator
    with open(filename) as f:
        lines = iter(f.readlines())

    if any([filename.lower().endswith(vasp.lower()) for vasp in ["OUTCAR", "CONTCAR", "vasp"]]) \
            or ftype.lower() == "vasp":
        # then file is vasp
        next(lines)  # skip first line
        # read cell par
        alat = float64(next(lines))
        par = array([next(lines).split()[0:3] for _ in range(3)], dtype=float64) * alat
        par_units = "angstrom"
        # read ions
        _ion = next(lines).split()
        _count = fromstring(next(lines), sep=' ', dtype=int)
        ion = []
        for i, c in zip(_ion, _count):
            ion += [i] * c
        nat = len(ion)
        # read pos type
        pos_units = "crystal" if next(lines).strip().lower() == "direct" else "angstrom"
        pos = array([next(lines).split()[0:3] for _ in range(nat)], dtype=float64)

    elif filename.lower().endswith("xyz") or ftype.lower() == "xyz":
        # then file is xyz format
        par_units = 'angstrom'
        # set par to a box with dimension of xyz_par
        par = xyz_par * eye(3)
        # first line is number of atoms
        nat = int(next(lines))
        # next line is a comment
        next(lines)
        pos_units = 'angstrom'
        # next nat lines are positions
        ion = []
        pos = []
        for _ in range(nat):
            nl = next(lines).split()[0:4]
            ion.append(nl[0])
            pos.append(nl[1:4])

    elif filename.lower().endswith("xsf") or ftype.lower() == "xsf":
        for line in lines:
            if line.startswith('PRIMVEC') or ' PRIMVEC' in line:
                # if statement above avoids reading 'RECIP-PRIMVEC'
                # read cell parameters
                par_units = 'angstrom'
                par = array([next(lines).split()[0:3] for _ in range(3)], dtype=float64)
            elif 'PRIMCOORD' in line:
                # read ions and pos
                pos_units = 'angstrom'
                _nat = int(next(lines).split()[0])
                ion, pos = [], []
                for _ in range(_nat):
                    nl = next(lines).split()
                    ion.append(element(int(nl[0])).symbol)
                    pos.append(nl[1:4])

    elif filename.lower().endswith("jdftx") or filename.lower().endswith("pos") \
            or ftype.lower() == "jdftx":
        # reread lines where all continuation lines have been resolved (i.e. lines that end in '\')
        lines = iter(_resolve_continuation_lines(filename))
        ion = []
        pos = []
        for line in lines:
            line = line.split('#')[0].strip()
            if line.startswith('lattice '):
                par = fromstring(line.strip('lattice'), dtype=float, sep=' ').reshape((3, 3)).T
            elif line.startswith('ion '):
                ion.append(line.split()[1])
                pos.append(array(line.split()[2:5], dtype=float))

        pos = array(pos)
        nat = len(ion)
        par_units = "bohr"
        warn("pos_units not set, defaulting to crystal")
        pos_units = "crystal"

    else:
        # then file is QE
        nml = f90nml.read(filename)
        if len(nml) == 0 or ftype.lower() == "qeout":
            # then file is output
            for line in lines:
                if "lattice parameter (alat)" in line:
                    alat = float64(line.split()[-2]) * bohr_to_angstrom
                elif "number of atoms/cell" in line:
                    nat = int(line.split()[-1])
                elif "crystal axes:" in line:
                    # cell parameters in angstrom
                    par_units = "angstrom"
                    par = array([next(lines).split()[3:6] for _ in range(3)], dtype=float64) * alat
                elif "ATOMIC_POSITIONS" in line:
                    pos_units = nonalpha.sub('', line.split('ATOMIC_POSITIONS')[1]).lower()
                    lines, pos, ion = _read_atomic_positions(lines, nat, no_if_pos=True)
        else:
            # ftype.lower() == "qeinp", no need to actually check
            # then file is input
            alat = None
            for line in lines:
                if 'a' in nml['system']:
                    alat = nml['system']['a']
                elif 'celldm' in nml['system']:
                    alat = nml['system']['celldm'][0] * bohr_to_angstrom
                nat = int(nml['system']['nat'])
                line = line.split('!')[0]   # trim away comments
                if 'CELL_PARAMETERS' in line:
                    par_units = nonalpha.sub('', line.split('CELL_PARAMETERS')[1]).lower()
                    par = array([fromstring(next(lines).split('!')[0], sep=' ')
                                 for _ in range(3)], dtype=float64)
                    # par = _convert_par(par, in_units=par_units, alat=alat)
                elif 'ATOMIC_POSITIONS' in line:
                    pos_units = nonalpha.sub('', line.split('ATOMIC_POSITIONS')[1]).lower()
                    lines, pos, ion = _read_atomic_positions(lines, nat, no_if_pos=True)

            if int(nml['system']['ibrav']) != 0:
                # if lattice is specified by ibrav then build par from ibrav
                # note this routine handles alat on its own
                par = _ibrav_to_par(nml['system'], units="angstrom")
                par_units = "angstrom"

            elif alat is not None:
                par *= alat

    return cls(ion=ion, par=par, pos=pos, pos_units=pos_units, par_units=par_units)  # pylint: disable=E1102


def write_file(self, filename, ftype='auto'):
    '''
    write atomgeo to file (QE format)
    '''
    if ftype == 'auto':
        if any([filename.lower().endswith(ext.lower()) for ext in ["OUTCAR", "CONTCAR", "vasp"]]):
            ftype = 'vasp'
        elif filename.lower().endswith("xyz"):
            ftype = 'xyz'
        elif filename.lower().endswith("xsf"):
            ftype = 'xsf'
        elif any([filename.lower().endswith(ext.lower()) for ext in ["pos", "jdftx"]]):
            ftype = 'jdftx'
        else:
            ftype = 'qe'
    elif ftype not in ['vasp', 'xyz', 'xsf', 'jdftx', 'qe']:
        raise ValueError('Value of ftype not recognized')

    with open(filename, 'w') as f:
        f.write(self.__str__(ftype=ftype))

    return None
