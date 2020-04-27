from numpy import array

# define if_pos property
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
    return self._ibrav


# TODO allow ibrav to be set?
# @ibrav.setter
# def ibrav(self, ibrav):
#     ''' set value of atomgeo._ibrav '''
#     self._ibrav = array(ibrav, dtype=int)


# define lattice A, B, C property
@property
def A(self):
    ''' return value of atomgeo._A '''
    return self._A


@A.setter
def A(self, A):
    ''' set value of atomgeo._A '''
    if hasattr(self, 'celldm'):
        raise ValueError("Cannot simultaneosly set A,B,C, when celldm is in use")
    self._A = float(A)


@property
def B(self):
    ''' return value of atomgeo._B '''
    return self._B


@B.setter
def B(self, B):
    ''' set value of atomgeo._B '''
    if not hasattr(self, 'A'):
        raise ValueError("Cannot set B if A is not in use")
    self._B = float(B)


@property
def C(self):
    ''' return value of atomgeo._C '''
    return self._C


@C.setter
def C(self, C):
    ''' set value of atomgeo._C '''
    if not hasattr(self, 'A'):
        raise ValueError("Cannot set C if A is not in use")
    self._C = float(C)


@property
def cosAB(self):
    ''' return value of atomgeo._cosAB '''
    return self._cosAB


@cosAB.setter
def cosAB(self, cosAB):
    ''' set value of atomgeo._cosAB '''
    if not hasattr(self, 'A'):
        raise ValueError("Cannot set cosAB if A is not in use")
    self._cosAB = float(cosAB)


@property
def cosBC(self):
    ''' return value of atomgeo._cosBC '''
    return self._cosBC


@cosBC.setter
def cosBC(self, cosBC):
    ''' set value of atomgeo._cosBC '''
    if not hasattr(self, 'A'):
        raise ValueError("Cannot set cosBC if A is not in use")
    self._cosBC = float(cosBC)


@property
def cosAC(self):
    ''' return value of atomgeo._cosAC '''
    return self._cosAC


@cosAC.setter
def cosAC(self, cosAC):
    ''' set value of atomgeo._cosAC '''
    if not hasattr(self, 'A'):
        raise ValueError("Cannot set cosAC if A is not in use")
    self._cosAC = float(cosAC)


# define lattice celldm property
@property
def celldm(self):
    ''' return value of atomgeo._celldm '''
    return self._celldm


@celldm.setter
def celldm(self, celldm):
    ''' set value of atomgeo._celldm '''
    if hasattr(self, 'A'):
        raise ValueError("Cannot simultaneosly set celldm, when A,B,C are in use")
    elif array(celldm).shape != self._celldm.shape:
        raise ValueError("Passed array is not of the same shape as the original array")
    self._celldm = array(celldm, dtype=float)
