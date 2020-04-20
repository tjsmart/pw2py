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
    if hasattr(self, 'celldm'):
        raise ValueError("Cannot simultaneosly set A,B,C, when celldm is in use")
    self._B = float(B)


@property
def C(self):
    ''' return value of atomgeo._C '''
    return self._C


@C.setter
def C(self, C):
    ''' set value of atomgeo._C '''
    if hasattr(self, 'celldm'):
        raise ValueError("Cannot simultaneosly set A,B,C, when celldm is in use")
    self._C = float(C)


# define lattice celldm property
@property
def celldm(self):
    ''' return value of atomgeo._celldm '''
    return self._celldm


@celldm.setter
def celldm(self, celldm):
    ''' set value of atomgeo._celldm '''
    if array(celldm).shape != self._celldm.shape:
        raise ValueError("Passed array is not of the same shape as the original array")
    elif hasattr(self, 'A'):
        raise ValueError("Cannot simultaneosly set celldm, when A,B,C are in use")
    self._celldm = array(celldm, dtype=float)
