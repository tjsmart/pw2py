class atomgeo:
    '''
    class for atomic geometry

    Suggested to make object from file:
    -----
    geo = atomgeo.from_file(filename, ftype="auto")

    Attributes:
    -----

    par (np.array, float, shape = (3,3))
        - cell parameters

    par_units (str)
        - units of par, acceptable values include "angstrom", "bohr", "alat"

    nat (int)
        - number of atoms

    ion (list, str)
        - list of atom/ion names (QE format)

    pos (np.array, float, shape = (nat, 3))
        - atomic positions

    pos_units (str)
        - units of pos, acceptable values include "angstrom", "bohr", "crystal", "alat"
    '''

    # dunder methods
    from .init import __init__
    from .str import __str__, __repr__
    # properties
    from .properties import ion, nat, par, par_units, pos, pos_units
    # io methods
    from .io import from_file, write_file
    # other methods
    from .atoms import add_atom, remove_indices, replace_ion, sort_ions
    from .cell import build_supercell, shift_pos_to_unit