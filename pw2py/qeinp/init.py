from .qecard import qecard
from .qenml import qenml
from .. import qegeo


def __init__(self, nml, geo, card, kpt):
    ''' initialize qeinp object '''

    '''
    object to be made of 3 components:
        nml -> all of the namelists attributes
        geo -> ion, pos, par, units
        card -> dictionary holding card options (ATOMIC_POSITIONS, etc.)
    '''
    # use qegeo to instantiate many properties
    qegeo.__init__(self, geo, geo.if_pos, geo._nml)

    delattr(self, '_nml')

    self.nml = qenml(nml)
    self.card = qecard(card)
    self._kpt = kpt
