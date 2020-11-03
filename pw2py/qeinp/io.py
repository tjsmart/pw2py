import f90nml
from numpy import fromstring, zeros
from warnings import warn

from ..functions import read_geo
from .._common.resource import _read_qe_card_option


@classmethod
def from_file(cls, filename, is_prefix=None):
    '''
    create qeinp object from file
    '''
    if is_prefix is None:
        # determine whether filename is prefix (i.e. read prefix.in and prefix.out)
        is_prefix = False if filename.endswith('.in') else True

    if is_prefix:
        if filename.endswith('.'):
            inpfile = filename + 'in'
            outfile = filename + 'out'
        elif filename.endswith('.out'):
            inpfile = filename[:-3] + 'in'
            outfile = filename
        else:
            inpfile = filename + '.in'
            outfile = filename + '.out'
    else:
        inpfile = filename

    # read namelist
    nml = f90nml.read(inpfile)
    assert len(nml) != 0, "Input file does not contain a valid namelist!"

    # read atomgeo
    if is_prefix:
        ion, par, par_units, pos, pos_units, if_pos = read_geo(outfile, ftype='qeout', read_if_pos=True)
    else:
        ion, par, par_units, pos, pos_units, if_pos = read_geo(inpfile, ftype='qeinp', read_if_pos=True)

    # acquire other things, i.e. card and kpoints
    # intialize card dictionary
    card = {}

    # intialize kpt (default is None)
    kpt = None

    # parse additional data from inpfile
    with open(inpfile) as f:
        # need to get: ['CELL_PARAMETERS', 'ATOMIC_SPECIES', 'ATOMIC_POSITIONS', 'K_POINTS']
        # also get if_pos and kpt
        for line in f:
            if 'CELL_PARAMETERS' in line:
                # TODO why not just use par_units from above?
                card['CELL_PARAMETERS'] = _read_qe_card_option(line, 'CELL_PARAMETERS')
            elif 'ATOMIC_SPECIES' in line:
                card['ATOMIC_SPECIES'] = {}
                for _ in range(int(nml['system']['ntyp'])):
                    nl = f.readline()
                    try:
                        symbol, mass, pp = nl.split('!')[0].split('#')[0].split()
                    except ValueError:
                        raise ValueError('Unable to unpack this line in ATOMIC_SPECIES: {}'.format(nl))
                    card['ATOMIC_SPECIES'][symbol] = [float(mass), str(pp)]
            elif 'ATOMIC_POSITIONS' in line:
                card['ATOMIC_POSITIONS'] = _read_qe_card_option(line, 'ATOMIC_POSITIONS')
            elif 'K_POINTS' in line:
                card['K_POINTS'] = _read_qe_card_option(line, 'K_POINTS')
                # read k-points
                if card['K_POINTS'] == 'automatic':
                    kpt = fromstring(f.readline().split('!')[0].split('#')[0], sep=' ', dtype=int).reshape((2, 3))
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
                        arr = fromstring(f.readline(), sep=' ')
                        occ[ispin, bands_read:bands_read + len(arr)] = arr
                        bands_read += len(arr)

                card['OCCUPATIONS'] = occ

            elif 'CONSTRAINTS' in line:
                warn('CONSTRAINTS not implemented and could not be read.')

            elif 'ATOMIC_FORCES' in line:
                warn('ATOMIC_FORCES not implemented and could not be read.')

    return cls(nml, card, ion, par, par_units, pos, pos_units, if_pos, kpt)  # pylint: disable=E1102


def write_file(self, filename, **kwargs):
    '''
    create qeinput object from file
    '''
    with open(filename, 'w') as f:
        f.write(str(self))

    return None
