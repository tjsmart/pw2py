import numpy as np

from .._common.constants import ry_to_ev


@staticmethod
def final_energy(filename, conv_level='automatic', units='ev'):
    '''
    grab the final energy from filename and return convergence level

    filename (str):
        - quantum espresso output file

    conv_level (str):
        - 'total energy'    : return last instance of total energy (could be from inner SCF loop)
        - '!'               : return last instance of total energy with SCF loop convergence
        - '!!'              : return last instance of total energy from EXX loop convergence
        - 'Final'           : return final energy (from relax/vc-relax)
        - 'automatic'       : return highest convergence level achieved

    returns
    ---
        tuple(energy, conv_level)
            - energy(float) = energy in Ry
            - conv_level(str) = level of convergence achieved
    '''
    assert units.lower() in ['ev', 'ry',
                             'ha'], "Only eV, Ry, and Ha are supported units"
    # build dictionary of energy instances
    energy_instances = {}
    instances = ['Final', '!!', '!', 'total energy']
    conv_level = conv_level.lower()
    # check given conv_level is a valid option
    if conv_level not in instances + ['automatic']:
        return None, "invalid conv_level: '{}'".format(conv_level)
    # load file to memory, storing in lines
    with open(filename) as f:
        lines = f.readlines()
    # iterate backwards through file
    for line in reversed(lines):
        # skip lines without 'energy' in them, or line related to sum of terms
        if ('energy' not in line) or ('sum of the following terms' in line):
            continue
        # otherwise iterate through possible instances and extract energy when found
        for instance in instances:
            if instance in line and instance not in energy_instances:
                try:
                    energy_instances[instance] = float(
                        line.split('=')[1].split('R')[0])
                except IndexError:
                    # false line, i.e. no total energy was given
                    pass
        # if all instances have been found then break loop
        if instances == list(energy_instances.keys()):
            break

    # now return the appropriate energy, based on user conv_level
    if conv_level == 'automatic':
        # for automatic iterate through instances and return first available instance
        for instance in instances:
            if instance in energy_instances:
                if units.lower() == 'ev':
                    return energy_instances[instance] * ry_to_ev, instance
                elif units.lower() == 'ry':
                    return energy_instances[instance], instance
                elif units.lower() == 'ha':
                    return energy_instances[instance] * 0.5, instance

    elif conv_level not in energy_instances:
        # user specified conv_level was not found in the file
        return None, "conv_level '{}' not achieved".format(conv_level)

    else:
        # return user specified conv_level
        if units == 'eV':
            return energy_instances[conv_level] * ry_to_ev, conv_level
        else:
            return energy_instances[conv_level], conv_level

    # this should not happen
    return None, 'Unexpected error in final_energy'


def _readEigs(lines, nbnd):
    '''
    (private) parse lines for eigenvalues
    '''
    next(lines)
    eigs = []
    for _ in range(nbnd // 8 + 1):
        eigs.append(np.fromstring(next(lines), sep=' ', dtype=np.float64))

    return np.hstack(eigs), lines


def calcEigs(self):
    '''
    calculate VBM, CBM, and gap (only works with smearing)
    '''
    try:
        fermi = self.conv['Fermi'][-1]
    except KeyError:
        # TODO -- implement no smearing case
        raise ValueError(
            'calcEigs not implemented for systems without smearing')

    eigs = self.list_eigs[-1]
    vbm = max(np.where(eigs < fermi, eigs, -1E10))
    cbm = min(np.where(eigs > fermi, eigs, 1E10))

    return vbm, cbm, cbm - vbm


# def _read_qeout_data(f, prec=9, dtype=np.float64):
def _read_qeout_data(f, prec=9):
    data = []
    while True:
        fline = f.readline()
        if len(fline.strip()) == 0:
            break
        # remove two spaces at beginning, new line character at end
        fline = fline[2:-1]
        # manually split fline
        fline_split = []
        for i in range(len(fline) // prec):
            # grab next prec characters in fline
            value = fline[prec*i:prec*(i+1)]
            try:
                float(value)
                fline_split.append(value)
            except ValueError:
                fline_split.append(np.nan)
        data += fline_split
    # return np.ndarray(data, dtype=dtype)
    return data


@staticmethod
def read_bands(filename: str):
    '''
    read and return final eigenvalues from file
    returns a dict(see below) including kpt but also occ and fermi if given

    params
    ---
        filename (str)

    returns
    ---
        collection (dict):
            keys: kpt (np.ndarray), fermi (float), eig (np.ndarray), occ (np.ndarray)

    Note
    ---
        if ispin = 2, eig/occ are indexed by (spin, kpt, band)
        else, eig/occ are indexed by (kpt, band)

        likewise, kpt are indexed by spin then kpt if ispin = 2

        eig will contain np.nan if eigenvalue overflowed, i.e. exceeded
        9 characters (less than -999.9999) or (greater than 9999.9999)

        not all keys written above may be present depending on output file
    '''
    nspin = False
    occ_given = False
    fermi_given = False
    vbm_given = False
    cbm_given = False
    with open(filename) as f:
        for line in f:
            if 'End of self-consistent calculation' in line:
                kpt = []
                read_eig = []
                read_occ = []
            elif 'SPIN UP' in line:
                nspin = True
                eig_up = []
                read_eig = eig_up
                occ_up = []
                read_occ = occ_up
            elif 'SPIN DOWN' in line:
                eig_dn = []
                read_eig = eig_dn
                occ_dn = []
                read_occ = occ_dn
            elif line.startswith('          k ='):
                # read the kpt
                kptstr = line.split('=')[-1].split('(')[0][:-1]
                k = np.array(
                    [kptstr[7*i:7*(i+1)] for i in range(3)], dtype=float
                )
                kpt.append(k)
                # skip a line
                f.readline()
                # read the eigenvalues
                read_eig.append(_read_qeout_data(f))
            elif line.startswith('     occupation numbers'):
                occ_given = True
                read_occ.append(_read_qeout_data(f))
            elif 'Fermi' in line:
                fermi_given = True
                fermi = float(line.split()[4])
            elif line.startswith('     highest occupied, lowest unoccupied level (ev):'):
                vbm_given = True
                vbm = float(line.split()[-2])
                cbm_given = True
                cbm = float(line.split()[-1])
            elif line.startswith('     highest occupied level (ev):'):
                vbm_given = True
                vbm = float(line.split()[-1])

    # collect values read into dictionary 'collection'
    collection = {}
    collection['nspin'] = 2 if nspin else 1
    collection['kpt'] = np.array(kpt, dtype=float).reshape(
        collection['nspin'], len(kpt) // collection['nspin'], 3
    )
    if fermi_given:
        collection['fermi'] = fermi
    if vbm_given:
        collection['vbm'] = vbm
    if cbm_given:
        collection['cbm'] = cbm
    if nspin:
        collection['eig'] = np.array([eig_up, eig_dn], dtype=float)
        if occ_given:
            collection['occ'] = np.array([occ_up, occ_dn], dtype=float)
    else:
        collection['eig'] = np.array([read_eig], dtype=float)
        if occ_given:
            collection['occ'] = np.array([read_occ], dtype=float)

    return collection
