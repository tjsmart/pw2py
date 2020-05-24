import numpy as np

from .qecard import qecard
from .qenml import qenml
from .. import atomgeo


def __init__(self, nml, card, ion, par, par_units, pos, pos_units, if_pos, kpt):
    ''' initialize qeinp object '''

    '''
    object to be made of 3 components:
        nml -> all of the namelists attributes
        geo -> ion, pos, par, units
        card -> dictionary holding card options (ATOMIC_POSITIONS, etc.)
    other:
        kpt -> specified kpoint grid
        if_pos -> if_pos for relaxing/fixing atomic positions in input file
    '''
    # use atomgeo to instantiate many properties
    atomgeo.__init__(self, ion=ion, par=par, pos=pos,
                     par_units=par_units, pos_units=pos_units)
    # set nml and card attributes
    self.nml = qenml(nml)
    self.card = qecard(card)
    # additional attributes
    self._kpt = kpt
    if if_pos is None:
        self._if_pos = np.ones((self.nat, 3), dtype=np.int32)
    else:
        self._if_pos = np.array(if_pos, dtype=np.int32)
