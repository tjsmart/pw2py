# defines all cards present in qe input data description


@property
def ATOMIC_SPECIES(self):
    return self._card['ATOMIC_SPECIES']


@ATOMIC_SPECIES.setter
def ATOMIC_SPECIES(self, ATOMIC_SPECIES):
    self._card['ATOMIC_SPECIES'] = dict(ATOMIC_SPECIES)


@property
def ATOMIC_POSITIONS(self):
    return self._card['ATOMIC_POSITIONS']


@ATOMIC_POSITIONS.setter
def ATOMIC_POSITIONS(self, ATOMIC_POSITIONS):
    if str(ATOMIC_POSITIONS).lower() not in \
            ['alat', 'bohr', 'angstrom', 'crystal', 'crystal_sg']:
        raise ValueError('Invalid value for ATOMIC_POSITIONS')
    self._card['ATOMIC_POSITIONS'] = str(ATOMIC_POSITIONS).lower()


@property
def K_POINTS(self):
    return self._card['K_POINTS']


@K_POINTS.setter
def K_POINTS(self, K_POINTS):
    if str(K_POINTS).lower() not in \
            ['tpiba', 'automatic', 'crystal', 'gamma', 'tpiba_b', 'crystal_b', 'tpiba_c', 'crystal_c']:
        raise ValueError('Invalid value for K_POINTS')
    self._card['K_POINTS'] = str(K_POINTS).lower()


@property
def CELL_PARAMETERS(self):
    return self._card['CELL_PARAMETERS']


@CELL_PARAMETERS.setter
def CELL_PARAMETERS(self, CELL_PARAMETERS):
    if str(CELL_PARAMETERS).lower() not in [' alat', 'bohr', 'angstrom']:
        raise ValueError('Invalid value for CELL_PARAMETERS')
    self._card['CELL_PARAMETERS'] = str(CELL_PARAMETERS).lower()


# CONSTRAINTS not implemented
# @property
# def CONSTRAINTS(self):
#     return self._card['CONSTRAINTS']


# @CONSTRAINTS.setter
# def CONSTRAINTS(self, CONSTRAINTS):
#     self._card['CONSTRAINTS'] = CONSTRAINTS


# OCCUPATIONS not implemented
@property
def OCCUPATIONS(self):
    return self._card['OCCUPATIONS']


@OCCUPATIONS.setter
def OCCUPATIONS(self, OCCUPATIONS):
    self._card['OCCUPATIONS'] = OCCUPATIONS


# @property
# def ATOMIC_FORCES(self):
#     return self._card['ATOMIC_FORCES']


# @ATOMIC_FORCES.setter
# def ATOMIC_FORCES(self, ATOMIC_FORCES):
#     self._card['ATOMIC_FORCES'] = ATOMIC_FORCES
