from numpy import array, unique
from warnings import warn

from .. import atomgeo


# define ntyp property
@property
def ntyp(self):
    ''' return number of atoms '''
    return unique(self.ion).size


@property
def kpt(self):
    ''' return value of qeinp._kpt '''
    return self._kpt


@kpt.setter
def kpt(self, kpt):
    ''' set value of qeinp._kpt, pass None to use gamma only'''
    if kpt is None:
        self._kpt = None
        self.card.K_POINTS = 'gamma'
    else:
        assert array(kpt).shape == (2, 3), \
            "Passed array does not have the correct shape (2, 3), passed: {}".format(array(kpt).shape)
        self.card.K_POINTS = 'automatic'
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
    ''' return value of self.nml.ibrav '''
    return self.nml.ibrav


@ibrav.setter
def ibrav(self, ibrav):
    ''' set value of self.nml.ibrav but also update CELL_PARAMETERS and namelist '''
    if ibrav == 0:
        self.nml.ibrav = 0
        self.card.CELL_PARAMETERS = self.par_units
        for attr in ['a', 'b', 'c', 'cosab', 'cosac', 'cosbc', 'celldm']:
            try:
                self.nml['system'].pop(attr)
            except KeyError:
                pass
    else:
        raise NotImplementedError('Changing ibrav to nonzero values is not implemented')
