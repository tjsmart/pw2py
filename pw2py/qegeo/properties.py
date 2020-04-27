from numpy import array


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
    return self._nml['ibrav']


# TODO allow ibrav to be set?
# @ibrav.setter
# def ibrav(self, ibrav):
#     ''' set value of atomgeo._ibrav '''
#     self._nml['ibrav'] = array(ibrav, dtype=int)


# define lattice A, B, C property
@property
def A(self):
    ''' return value of atomgeo._A '''
    return self._nml['a']


@A.setter
def A(self, A):
    ''' set value of atomgeo._A '''
    self._nml['a'] = float(A)


@property
def B(self):
    ''' return value of atomgeo._B '''
    return self._nml['b']


@B.setter
def B(self, B):
    ''' set value of atomgeo._B '''
    self._nml['b'] = float(B)


@property
def C(self):
    ''' return value of atomgeo._C '''
    return self._nml['c']


@C.setter
def C(self, C):
    ''' set value of atomgeo._C '''
    self._nml['c'] = float(C)


@property
def cosAB(self):
    ''' return value of atomgeo._cosAB '''
    return self._nml['cosab']


@cosAB.setter
def cosAB(self, cosAB):
    ''' set value of atomgeo._cosAB '''
    self._nml['cosab'] = float(cosAB)


@property
def cosBC(self):
    ''' return value of atomgeo._cosBC '''
    return self._nml['cosbc']


@cosBC.setter
def cosBC(self, cosBC):
    ''' set value of atomgeo._cosBC '''
    self._nml['cosbc'] = float(cosBC)


@property
def cosAC(self):
    ''' return value of atomgeo._cosAC '''
    return self._nml['cosac']


@cosAC.setter
def cosAC(self, cosAC):
    ''' set value of atomgeo._cosAC '''
    self._nml['cosac'] = float(cosAC)


# define lattice celldm property
@property
def celldm(self):
    ''' return value of atomgeo._celldm '''
    return self._nml['celldm']


@celldm.setter
def celldm(self, celldm):
    ''' set value of atomgeo._celldm '''
    if array(celldm).shape != self._nml['celldm'].shape:
        raise ValueError("Passed array is not of the same shape as the original array")
    self._nml['celldm'] = array(celldm, dtype=float)
