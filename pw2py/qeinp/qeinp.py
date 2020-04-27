from .. import qegeo


class qeinp(qegeo):
    '''
    class for quantum espresso input file
        qedict = dictionary of namelists (nml), ATOMIC_SPECIES, CELL_PARAMETERS, ATOMIC_POSITIONS, and K_POINTS
        par = CELL_PARAMETERS data
        ion = ATOMIC_POSITIONS ions data
        pos = ATOMIC_POSITIONS position data
        if_pos = ATOMIC_POSITIONS if_pos data
        kpt = K_POINTS data
    '''

    from ._init import __init__

    from ._str import __repr__, __str__

    from ._io import from_file, write_file

    from ._properties import kpt, par_units, pos_units, ion, par, pos, A, B, C, celldm, \
        cosAB, cosAC, cosBC, ibrav, if_pos

    # TODO
    # def convert_ibrav(self, newBrav):
    #     '''
    #     convert ibrav of self to newBrav (intended for ibrav != 0 to 0 or back)
    #     '''
    #     pass
