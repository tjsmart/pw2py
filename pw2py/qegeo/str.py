from f90nml import Namelist
from numpy import array_equal, ones


def __repr__(self):
    '''
    convert qegeo to str (QE format)
    '''
    # TODO need to be able to handle cases with ibrav != 0
    out = "&control\n/\n\n"
    out += "{}\n\n".format(Namelist({'system': self._nml}))
    out += "&electrons\n/\n\n"
    if self.ibrav == 0:
        out += "CELL_PARAMETERS {}\n".format(self.par_units)
        for par in self.par:
            out += "    {:16.9f}  {:16.9f}  {:16.9f}\n".format(par[0], par[1], par[2])
        out += "\n"
    out += "ATOMIC_POSITIONS {}\n".format(self.pos_units)
    for ion, pos, if_pos in zip(self.ion, self.pos, self.if_pos):
        out += "    {:5s}  {:16.9f}  {:16.9f}  {:16.9f}".format(ion, pos[0], pos[1], pos[2])
        if array_equal(if_pos, ones(3, dtype=int)):
            out += "\n"
        else:
            try:
                out += "    {}  {}  {}\n".format(if_pos[0], if_pos[1], if_pos[2])
            except IndexError:
                None
    out += "\n"

    return out


def __str__(self):
    return repr(self)
