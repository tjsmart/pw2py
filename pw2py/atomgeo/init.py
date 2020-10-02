from numpy import array


def __init__(self, ion=None, par=None, pos=None, par_units=None, pos_units=None):
    ''' initialize atomgeo instance '''
    self._ion = array(ion, dtype=object)
    self._par = array(par, dtype=float)
    self._pos = array(pos, dtype=float)
    self._par_units = str(par_units).lower()
    self._pos_units = str(pos_units).lower()
