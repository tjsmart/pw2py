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

    from .init import __init__

    from .str import __repr__, to_string

    from .io import from_file, write_file

    from .properties import kpt, par_units, pos_units, ion, par, pos, A, B, C, celldm, \
        cosAB, cosAC, cosBC, ibrav, if_pos, ntyp

    from .methods import load_geo, replace_ion, add_atom, remove_indices
