class qecard:
    '''
    class for handling quantum espresso cards
        > ATOMIC_SPECIES, ATOMIC_POSITIONS, K_POINTS, CELL_PARAMETERS, \
        CONSTRAINTS, OCCUPATIONS, ATOMIC_FORCES
    '''

    from .init import __init__

    from .properties import ATOMIC_SPECIES, ATOMIC_POSITIONS, K_POINTS, CELL_PARAMETERS, \
        OCCUPATIONS
    # TODO
    # CONSTRAINTS, ATOMIC_FORCES

    from .str import __repr__
