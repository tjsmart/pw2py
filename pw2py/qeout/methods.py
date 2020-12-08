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
    assert units.lower() in ['ev', 'ry', 'ha'], "Only eV, Ry, and Ha are supported units"
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
                    energy_instances[instance] = float(line.split('=')[1].split('R')[0])
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
        raise ValueError('calcEigs not implemented for systems without smearing')

    eigs = self.list_eigs[-1]
    vbm = max(np.where(eigs < fermi, eigs, -1E10))
    cbm = min(np.where(eigs > fermi, eigs, 1E10))

    return vbm, cbm, cbm - vbm
