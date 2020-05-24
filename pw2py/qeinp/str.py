from numpy import array_equal, ones


def to_string(self):
    '''
    convert qeinp to string
    '''
    out = str(self.nml)
    out += "\n\n"
    out += "ATOMIC_SPECIES\n"
    for spec in self.card.ATOMIC_SPECIES:
        out += "    {:5s}  {:8.4f}  {}\n".format(spec, *self.card.ATOMIC_SPECIES[spec])
    out += "\n"
    if self.nml.ibrav == 0:
        out += "CELL_PARAMETERS {}\n".format(self.card.CELL_PARAMETERS)
        for par in self.par:
            out += "    {:16.9f}  {:16.9f}  {:16.9f}\n".format(*par)
        out += "\n"
    out += "ATOMIC_POSITIONS {}\n".format(self.card.ATOMIC_POSITIONS)
    for ion, pos, if_pos in zip(self.ion, self.pos, self.if_pos):
        out += "    {:5s}  {:16.9f}  {:16.9f}  {:16.9f}".format(ion, pos[0], pos[1], pos[2])
        if array_equal(if_pos, ones(3, dtype=int)):
            out += "\n"
        else:
            out += "    {:1d}  {:1d}  {:1d}\n".format(*if_pos)
    out += "\n"
    if 'K_POINTS' in self.card._card:
        out += "K_POINTS {}\n".format(self.card.K_POINTS)
        if self.card.K_POINTS == "automatic":
            out += "    {}  {}  {}    {}  {}  {}\n".format(*self.kpt[0], *self.kpt[1])
        out += "\n"
    if 'OCCUPATIONS' in self.card._card:
        out += "OCCUPATIONS\n"
        row_length = 8
        for ispin in range(self.card.OCCUPATIONS.shape[0]):
            iband = 0
            while iband < self.nml.nbnd:
                row_count = 0
                while row_count < row_length:
                    if row_count == 0:
                        out += "  "
                    out += "  {:6.4f}".format(self.card.OCCUPATIONS[ispin, iband])
                    iband += 1
                    row_count += 1
                    if iband == self.nml.nbnd:
                        break
                out += "\n"
        out += "\n"

    return out


def __repr__(self):
    return to_string(self)
