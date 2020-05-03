from numpy import array
from warnings import warn

from .. import atomgeo


@property
def kpt(self):
    ''' return value of atomgeo._kpt '''
    return self._kpt


@kpt.setter
def kpt(self, kpt):
    ''' set value of atomgeo._kpt '''
    if array(kpt).shape != self.kpt.shape:
        raise ValueError("Passed array is not of the same shape as the original array")
    self._kpt = array(kpt, dtype=int)


# extend functionality of par_units
@property
def par_units(self):
    # use same getter method from atomgeo
    return atomgeo.par_units.fget(self)


@par_units.setter
def par_units(self, par_units):
    # use setter method from atomgeo
    atomgeo.par_units.fset(self, par_units)
    self.card.CELL_PARAMETERS = self._par_units
    if self.ibrav != 0:
        # TODO how to handle this case?
        warn('Be careful changing par_units when using ibrav != 0')


# define pos_units property
@property
def pos_units(self):
    ''' return value of atomgeo._pos_units '''
    return self._pos_units


@pos_units.setter
def pos_units(self, pos_units):
    '''
    Convert atomic positions to 'units'
    updates values of: pos, pos_units, (card.ATOMIC_POSITIONS)
    '''
    atomgeo.pos_units.fset(self, pos_units)
    self.card.ATOMIC_POSITIONS = self._pos_units


# Place holders to extend functionality of different properties acquired from atomgeo
# otherwise these lines do not actually do anything new
@property
def ion(self):
    # use same getter method from atomgeo
    return atomgeo.ion.fget(self)


@ion.setter
def ion(self, ion):
    # use setter method from atomgeo
    atomgeo.ion.fset(self, ion)


@property
def pos(self):
    # use same getter method from atomgeo
    return atomgeo.pos.fget(self)


@pos.setter
def pos(self, pos):
    # use setter method from atomgeo
    atomgeo.pos.fset(self, pos)


@property
def par(self):
    # use same getter method from atomgeo
    return atomgeo.par.fget(self)


@par.setter
def par(self, par):
    # use setter method from atomgeo
    atomgeo.par.fset(self, par)
    if self.ibrav != 0:
        # TODO how to handle this case?
        warn('Be careful changing par when using ibrav != 0')


@property
def if_pos(self):
    ''' return value of atomgeo._if_pos '''
    return self._if_pos


@if_pos.setter
def if_pos(self, if_pos):
    ''' set value of atomgeo._if_pos '''
    if array(if_pos).shape != self._if_pos.shape:
        raise ValueError("Passed array is not of the same shape as the original array")
    self._if_pos = array(if_pos, dtype=int)


# define ibrav property
@property
def ibrav(self):
    ''' return value of atomgeo._ibrav '''
    return self.nml.ibrav


@ibrav.setter
def ibrav(self, ibrav):
    ''' set value of atomgeo._ibrav '''
    if ibrav == 0:
        self.nml.ibrav = 0
        self.card.CELL_PARAMETERS = self._par_units
        for attr in ['a', 'b', 'c', 'cosab', 'cosac', 'cosbc', 'celldm']:
            try:
                self.nml._nml['system'].pop(attr)
            except KeyError:
                pass
    else:
        raise ValueError('Changing ibrav to nonzero values is not implemented')


# define lattice A, B, C property
@property
def A(self):
    ''' return value of atomgeo._A '''
    return self.nml.A


@A.setter
def A(self, A):
    ''' set value of atomgeo._A '''
    self.nml.A = float(A)


@property
def B(self):
    ''' return value of atomgeo._B '''
    return self.nml.B


@B.setter
def B(self, B):
    ''' set value of atomgeo._B '''
    self.nml.B = float(B)


@property
def C(self):
    ''' return value of atomgeo._C '''
    return self.nml.C


@C.setter
def C(self, C):
    ''' set value of atomgeo._C '''
    self.nml.C = float(C)


@property
def cosAB(self):
    ''' return value of atomgeo._cosAB '''
    return self.nml.cosAB


@cosAB.setter
def cosAB(self, cosAB):
    ''' set value of atomgeo._cosAB '''
    self.nml.cosAB = float(cosAB)


@property
def cosBC(self):
    ''' return value of atomgeo._cosBC '''
    return self.nml.cosBC


@cosBC.setter
def cosBC(self, cosBC):
    ''' set value of atomgeo._cosBC '''
    self.nml.cosBC = float(cosBC)


@property
def cosAC(self):
    ''' return value of atomgeo._cosAC '''
    return self.nml.cosAC


@cosAC.setter
def cosAC(self, cosAC):
    ''' set value of atomgeo._cosAC '''
    self.nml.cosAC = float(cosAC)


# define lattice celldm property
@property
def celldm(self):
    ''' return value of atomgeo._celldm '''
    return self.nml.celldm


@celldm.setter
def celldm(self, celldm):
    ''' set value of atomgeo._celldm '''
    if array(celldm).shape != self.nml.celldm.shape:
        raise ValueError("Passed array is not of the same shape as the original array")
    self.nml.celldm = array(celldm, dtype=float)
