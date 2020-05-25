class qecard(dict):
    '''
    class for handling quantum espresso cards
        > is a subclass of dict
        > properties include:
            ATOMIC_SPECIES
            ATOMIC_POSITIONS
            K_POINTS
            CELL_PARAMETERS
            OCCUPATIONS

    to be implemented:
            CONSTRAINTS
            ATOMIC_FORCES
    '''

    from .properties import ATOMIC_SPECIES, ATOMIC_POSITIONS, K_POINTS, CELL_PARAMETERS, \
        OCCUPATIONS
    # TODO
    # CONSTRAINTS, ATOMIC_FORCES
