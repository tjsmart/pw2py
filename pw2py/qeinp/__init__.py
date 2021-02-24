from .. import atomgeo


class qeinp(atomgeo):
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

    from .str import __repr__, __str__, to_string

    from .io import from_file, write_file

    from .properties import kpt, par_units, pos_units, ion, par, pos, ibrav, if_pos, ntyp

    from .methods import load_geo, replace_ion, add_atom, drop_indices, build_supercell, to_atomgeo, _fix_species
