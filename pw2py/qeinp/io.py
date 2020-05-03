import f90nml
from numpy import fromstring, zeros
from warnings import warn

from .. import qegeo
from .._common.constants import nonalpha


@classmethod
def from_file(cls, filename):
    '''
    create qeinp object from file
    '''
    # read namelist
    nml = f90nml.read(filename)

    # read qegeo
    geo = qegeo.from_file(filename, ftype='qeinp')

    # acquire other things, i.e. card and kpoints
    # intialize card dictionary
    card = {}

    # intialize kpt (default is None)
    kpt = None

    # store lines of file to iterator
    with open(filename) as f:
        lines = iter(f.readlines())

    # loop through lines for ['CELL_PARAMETERS', 'ATOMIC_SPECIES', 'ATOMIC_POSITIONS', 'K_POINTS']
    for line in lines:
        line = line.split('!')[0].split('#')[0]   # trim away comments
        if 'CELL_PARAMETERS' in line:
            card['CELL_PARAMETERS'] = nonalpha.sub('', line.split('CELL_PARAMETERS')[1]).lower()
            # # read cell parameters
            # par = array([fromstring(next(lines).split('!')[0].split('#')[0], sep=' ') \
            # for _ in range(3)], dtype=float64)

        elif 'ATOMIC_SPECIES' in line:
            card['ATOMIC_SPECIES'] = {}
            for _ in range(int(nml['system']['ntyp'])):
                nl = next(lines)
                try:
                    symbol, mass, pp = nl.split('!')[0].split('#')[0].split()
                except ValueError:
                    raise ValueError('Unable to unpack this line in ATOMIC_SPECIES: {}'.format(nl))
                card['ATOMIC_SPECIES'][symbol] = [float(mass), str(pp)]

        elif 'ATOMIC_POSITIONS' in line:
            card['ATOMIC_POSITIONS'] = nonalpha.sub('', line.split('ATOMIC_POSITIONS')[1]).lower()

        elif 'K_POINTS' in line:
            card['K_POINTS'] = nonalpha.sub('', line.split('K_POINTS')[1]).lower()
            # read k-points
            if card['K_POINTS'] == 'automatic':
                kpt = fromstring(next(lines).split('!')[0].split('#')[0], sep=' ', dtype=int).reshape((2, 3))
            elif card['K_POINTS'] == 'gamma':
                kpt = None
            else:
                # TODO add support for crystal_b etc.
                kpt = None
                warn("K_POINTS option '{}' not supported, please use 'automatic' or 'gamma'")

        elif 'OCCUPATIONS' in line:
            try:
                nbnd = nml['system']['nbnd']
            except AttributeError:
                raise ValueError('Unable to read OCCUPATIONS when nbnd is not set.')
            try:
                nspin = nml['system']['nspin']
            except AttributeError:
                nspin = 1
            occ = zeros((nspin, nbnd))
            for ispin in range(nspin):
                bands_read = 0
                while bands_read < nml['system']['nbnd']:
                    arr = fromstring(next(lines), sep=' ')
                    occ[ispin, bands_read:bands_read + len(arr)] = arr
                    bands_read += len(arr)

            card['OCCUPATIONS'] = occ

        elif 'CONSTRAINTS' in line:
            warn('CONSTRAINTS not implemented and could not be read.')

        elif 'ATOMIC_FORCES' in line:
            warn('ATOMIC_FORCES not implemented and could not be read.')

    return cls(nml, geo, card, kpt)  # pylint: disable=E1102


def write_file(self, filename):
    '''
    create qeinput object from file
    '''
    with open(filename, 'w') as f:
        f.write(str(self))

    return None
