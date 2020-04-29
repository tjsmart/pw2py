import numpy as np

from .methods import _readEigs
from .._common.constants import bohr_to_angstrom, nonalpha


@classmethod
def from_file(cls, filename):
    '''
    create qeout object from file
    '''
    # store lines of file to iterator
    with open(filename) as f:
        lines = iter(f.readlines())

    is_mag = False
    is_exx = False
    conv = {}
    for key in ['E', '!', '!!', 'dE', 'nsteps', 'nsteps_exx', 'tot_forc', 'max_forc',
                'mag_tot', 'mag_abs', 'time', 'time_exx', 'tot_mag', 'abs_mag', 'Fermi']:
        conv[key] = []
    list_pos, if_pos, ion, list_eigs = [], [], [], []

    is_first = True
    for line in lines:
        if "lattice parameter (alat)" in line:
            alat = np.float64(line.split()[-2]) * bohr_to_angstrom
        elif "number of atoms/cell" in line:
            nat = int(line.split()[-1])
        elif "number of electrons" in line:
            # TODO use nelec ?
            nelec = int(float(line.split()[-1]))    # noqa: F841
        elif "number of Kohn-Sham states" in line:
            nbnd = int(line.split()[-1])
        elif "number of atomic types" in line:
            ntyp = int(line.split()[-1])
        elif "EXX-fraction" in line:
            # TODO use is_exx ?
            is_exx = True   # noqa: F841
        elif "Starting magnetic structure" in line:
            is_mag = True
        elif "crystal axes:" in line:
            par = np.array([next(lines)[0].split()[3:6] for _ in range(3)], dtype=np.float64) * alat
        elif "total cpu time spent up to now is" in line:
            conv['time'].append(np.float64(line.split()[-2]))
        elif "total energy              =" in line:
            if "!!" in line:
                conv['!!'].append(np.float64(line.split()[-2]))
                conv['nsteps_exx'].append(sum(conv['nsteps']) - sum(conv['nsteps_exx']))
                conv['time_exx'].append(conv['time'][-1])
            elif "!" in line:
                conv['!'].append(np.float64(line.split()[-2]))
                next(lines)
                conv['dE'].append(np.float64(next(lines).split()[-2]))
                if is_mag:
                    # skip throught till total magnetization
                    while "total magnetization" not in line:
                        line = next(lines)
                    conv['tot_mag'].append(np.float64(line.split()[-3]))
                    conv['abs_mag'].append(np.float64(next(lines).split()[-3]))

            else:
                conv['E'].append(np.float64(line.split()[-2]))

        elif "has" in line:
            conv['nsteps'].append(int(line.split()[-2]))
        elif "SPIN UP" in line:
            [next(lines) for _ in range(2)]
            if "k = 0.0000 0.0000 0.0000" in line:
                # TODO need to other kpoints
                eigs, lines = _readEigs(lines, nbnd)
                list_eigs.append(eigs)

        elif "SPIN DOWN" in line:
            [next(lines) for _ in range(2)]
            if "k = 0.0000 0.0000 0.0000" in line:
                eigs, lines = _readEigs(lines, nbnd)
                list_eigs.append(eigs)

        elif "k = 0.0000 0.0000 0.0000" in line:
            eigs, lines = _readEigs(lines, nbnd)
            list_eigs.append(eigs)

        elif "Fermi" in line:
            conv['Fermi'].append(np.float64(line.split()[-2]))
        elif "Forces acting on atoms" in line:
            max_forc = -1.0
            next(lines)
            for _ in range(nat):
                this_forc = np.linalg.norm(np.fromstring(
                    next(lines).split("force =")[1], sep=' ', dtype=np.float64))
                if this_forc > max_forc:
                    max_forc = this_forc
            conv['max_forc'].append(max_forc)
        elif "Total force" in line:
            conv['tot_forc'].append(np.float64(line.split()[3]))
        elif "ATOMIC_POSITIONS" in line:
            if is_first:
                pos_units = nonalpha.sub('', line.split('ATOMIC_POSITIONS')[1]).lower()
            pos = []
            for _ in range(nat):
                nl = next(lines).split()
                if is_first:
                    ion.append(nl[0])
                    try:
                        if_pos.append(np.array([nl[4], nl[5], nl[6]], dtype=int))
                    except IndexError:
                        if_pos.append(np.ones(3, dtype=int))
                pos.append(np.array(nl[1:4], dtype=np.float64))
            if is_first:
                is_first = False
            list_pos.append(np.array(pos))

    return cls(nat, ntyp, conv, par, ion, list_pos, if_pos, pos_units, list_eigs)
