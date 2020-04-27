from numpy import array, ones
from .. import atomgeo


# def __init__(self, geo, if_pos=None, ibrav=None, A=None, B=None, C=None,
#              cosAB=None, cosAC=None, cosBC=None, celldm=None):
def __init__(self, geo, if_pos, nml):
    ''' initialize qegeo instance '''
    # use atomgeo init method
    atomgeo.__init__(self, ion=geo.ion, par=geo.par, pos=geo.pos, par_units=geo.par_units, pos_units=geo.pos_units)

    # set additional attributes
    if if_pos is None:
        self._if_pos = ones((self.nat, 3), dtype=int)
    else:
        self._if_pos = array(if_pos, dtype=int)

    self._nml = nml

    # if (ibrav is None) or (ibrav == 0):
    #     self._ibrav = 0
    # else:
    #     self._ibrav = int(ibrav)
    #     if (A is not None) and (celldm is not None):
    #         raise ValueError("Only A,B,C or celldm can be used for one instance.")
    #     elif (A is None) and (celldm is None):
    #         raise ValueError("ibrav !=0 so either A,B,C or celldm must be set.")
    #     elif A is not None:
    #         self._A = float(A)
    #         if B is not None:
    #             self._B = float(B)
    #         if C is not None:
    #             self._C = float(C)
    #     else:
    #         self._celldm = array(celldm, dtype=float)
