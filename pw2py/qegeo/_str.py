from numpy import array_equal, ones


def __str__(self):
    '''
    convert qegeo to str (QE format)
    '''
    # TODO need to be able to handle cases with ibrav != 0
    ntyp = len(set(self.ion))
    out = "&control\n/\n&system\n    ibrav = 0\n    ntyp = {}\n    nat = {}\n/\n&electrons\n/\n"\
        .format(ntyp, str(self.nat))
    out += "CELL_PARAMETERS {}\n".format(self.par_units)
    for par in self.par:
        out += "    {:16.9f}  {:16.9f}  {:16.9f}\n".format(par[0], par[1], par[2])
    out += "ATOMIC_POSITIONS {}\n".format(self.pos_units)
    for ion, pos, if_pos in zip(self.ion, self.pos, self.if_pos):
        out += "    {:5s}  {:16.9f}  {:16.9f}  {:16.9f}\n".format(ion, pos[0], pos[1], pos[2])
        if array_equal(if_pos, ones(3, dtype=int)):
            out += "\n"
        else:
            try:
                out += "    {}  {}  {}\n".format(if_pos[0], if_pos[1], if_pos[2])
            except IndexError:
                None

    return out
