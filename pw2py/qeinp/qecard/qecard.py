class qecard():
    '''
    class for handling quantum espresso cards
        > ATOMIC_SPECIES, ATOMIC_POSITIONS, K_POINTS, CELL_PARAMETERS, \
        CONSTRAINTS, OCCUPATIONS, ATOMIC_FORCES
    '''

    from ._init import __init__

    from ._properties import ATOMIC_SPECIES, ATOMIC_POSITIONS, K_POINTS, CELL_PARAMETERS
    # TODO
    # CONSTRAINTS, OCCUPATIONS, ATOMIC_FORCES
