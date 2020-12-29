from copy import deepcopy
from numpy import array, cross

from .._common.resource import _convert_par, _convert_pos, _calc_rec


# Warning alat not fully supported, use with caution!
_valid_par_units = ['angstrom', 'bohr', 'alat']
_valid_pos_units = ['angstrom', 'bohr', 'crystal', 'alat']


# define ion property
@property
def ion(self):
    ''' return value of atomgeo._ion '''
    return self._ion


@ion.setter
def ion(self, ion, ion_units=None):
    ''' set value of atomgeo._ion '''
    if array(ion).shape != self._ion.shape:
        raise ValueError(
            "Passed array is not of the same shape as the original array")
    self._ion = array(ion, dtype=str)


# define par property
@property
def par(self):
    ''' return value of atomgeo._par '''
    return self._par


@par.setter
def par(self, par):
    ''' set value of atomgeo._par '''
    if array(par).shape != self.par.shape:
        raise ValueError("Cell parameters must be an array of shape (3,3)")
    self._par = array(par, dtype=float)
    self._rec = _calc_rec(self._par)


# define rec property
@property
def rec(self):
    ''' return value of atomgeo._rec '''
    return self._rec


# @rec.setter
# no rec setter, user must set par instead which will update rec


# define pos property
@property
def pos(self):
    ''' return value of atomgeo._pos '''
    return self._pos


@pos.setter
def pos(self, pos, pos_units=None):
    ''' set value of atomgeo._pos '''
    if array(pos).shape != self.pos.shape:
        raise ValueError(
            "Passed array is not of the same shape as the original array")
    self._pos = array(pos, dtype=float)


# define atoms property
@property
def atoms(self):
    ''' return atoms as dictionary with keys as ion names and values as positions '''
    return [{'ion': ion, 'pos': pos} for ion, pos in zip(self.ion, self.pos)]


# define par_units property
@property
def par_units(self):
    ''' return value of atomgeo._par_units '''
    return self._par_units


@par_units.setter
def par_units(self, units, inplace=True):
    '''
    Convert cell parameters to 'units'
    updates values of: par, par_units, (qedict['CELL_PARAMETERS'])
    '''
    units = str(units).lower()
    if units not in _valid_par_units:
        raise ValueError("Invalid value for units: {}".format(units))
    out = self if inplace else deepcopy(self)
    out.par = _convert_par(out.par, out.par_units, out_units=units)
    out._par_units = units
    return out


# define pos_units property
@property
def pos_units(self):
    ''' return value of atomgeo._pos_units '''
    return self._pos_units


@pos_units.setter
def pos_units(self, units, inplace=True):
    '''
    Convert atomic positions to 'units'
    updates values of: pos, pos_units, (qedict['ATOMIC_POSITIONS'])
    '''
    units = str(units).lower()
    if units not in _valid_pos_units:
        raise ValueError("Invalid value for units: {}".format(units))
    out = self if inplace else deepcopy(self)
    out.pos = _convert_pos(
        out.pos, out.pos_units, out_units=units, par=out.par, par_units=out.par_units)
    out._pos_units = units
    return out


# define nat property
@property
def nat(self):
    ''' return number of atoms '''
    return self._ion.size


# define vol property
@property
def vol(self):
    ''' return volume of the cell '''
    return cross(self.par[0], self.par[1]).dot(self.par[2])
