from numpy import array


def __init__(self, ion=None, par=None, pos=None, par_units=None, pos_units=None,
             if_pos=None, ibrav=None, A=None, B=None, C=None, celldm=None):
    ''' initialize qegeo instance '''

    super().__init__(ion=ion, par=par, pos=pos, par_units=par_units, pos_units=pos_units)
    self._if_pos = array(if_pos, dtype=int)
    self._ibrav = int(ibrav)
    if (A is not None) and (celldm is not None):
        raise ValueError("Only A,B,C or celldm can be used for one instance.")
    elif A is not None:
        self._A = float(A)
        if B is not None:
            self._B = float(B)
        if C is not None:
            self._C = float(C)
    else:
        self._celldm = array(celldm, dtype=float) if celldm is not None else None
