from collections import Counter
from copy import deepcopy


def to_string(self, ftype='qeinp'):
    '''
    convert atomgeo to str (default is qeinp format)
    '''
    ftype = ftype.lower()

    if ftype == 'qeinp':
        # TODO need to be able to handle cases with ibrav != 0 (maybe this should be qegeo specific)
        _ntyp = len(set(self.ion))
        out = "&control\n/\n\n&system\n    ibrav = 0\n    ntyp = {}\n    nat = {}\n/\n\n&electrons\n/\n\n"\
            .format(_ntyp, str(self.nat))
        out += "CELL_PARAMETERS {}\n".format(self.par_units)
        for par in self.par:
            out += "    {:16.9f}  {:16.9f}  {:16.9f}\n".format(par[0], par[1], par[2])
        out += "\n"
        out += "ATOMIC_POSITIONS {}\n".format(self.pos_units)
        for ion, pos in zip(self.ion, self.pos):
            # print(ion, pos)
            out += "    {:5s}  {:16.9f}  {:16.9f}  {:16.9f}\n".format(ion, pos[0], pos[1], pos[2])
        out += "\n"

    elif ftype == 'vasp':
        geo = self.sort_ions(inplace=False)
        # header
        out = "Generated from the PW2PY python module\n"
        # alat
        # out += "{}\n".format(float(alat))
        out += "{}\n".format('1.0')
        # par
        for par in geo.par:
            out += "    {:16.9f}  {:16.9f}  {:16.9f}\n".format(par[0], par[1], par[2])
        # collect unique ions and count of each
        ionAndCount = Counter(geo.ion)
        # ions
        out += "     {}\n".format("     ".join(ionAndCount.keys()))
        # count
        out += "     {}\n".format("     ".join(map(str, ionAndCount.values())))
        # pos type
        if geo.pos_units == "angstrom":
            out += "{}\n".format("cartesian")
        elif geo.pos_units == "crystal":
            out += "{}\n".format("direct")
        # pos
        for pos in geo.pos:
            out += "    {:16.9f}  {:16.9f}  {:16.9f}\n".format(pos[0], pos[1], pos[2])

    elif ftype == 'jdftx':
        _if_pos = 1
        geo = deepcopy(self)
        # change par to bohr
        geo.par_units = 'bohr'
        # change pos angstrom to bohr
        if geo.pos_units == 'angstrom':
            geo.pos_units = 'bohr'
        # write header
        out = "# Generated from the PW2PY python module\n"
        # write lattice
        out += "lattice \\\n"
        for i, par in enumerate(geo.par.T):
            if i == 0 or i == 1:
                out += "    {:16.9f}  {:16.9f}  {:16.9f} \\\n".format(par[0], par[1], par[2])
            elif i == 2:
                out += "    {:16.9f}  {:16.9f}  {:16.9f}\n".format(par[0], par[1], par[2])
        # write ion and pos
        for ion, pos in zip(geo.ion, geo.pos):
            out += "ion  {:5s}  {:16.9f}  {:16.9f}  {:16.9f} {}\n".format(ion, pos[0], pos[1], pos[2], _if_pos)

    elif ftype == 'xyz':
        geo = deepcopy(self)
        geo.pos_units = 'angstrom'
        # nat
        out = "{}\n".format(geo.nat)
        # description
        out += "Generated from the PW2PY python module\n"
        for ion, pos in zip(geo.ion, geo.pos):
            out += "    {:5s}  {:16.9f}  {:16.9f}  {:16.9f}\n".format(ion, pos[0], pos[1], pos[2])

    elif ftype == 'xsf':
        geo = deepcopy(self)
        out = "# Generated from the PW2PY python module\n"
        out += "DIM-GROUP\n"
        out += "    {:5d}    {:3d}\n".format(3, 1)
        # par
        geo.par_units = 'angstrom'
        for _section in ["PRIMVEC\n", "CONVVEC\n"]:
            out += _section
            for par in geo.par:
                out += "    {:16.9f}  {:16.9f}  {:16.9f}\n".format(par[0], par[1], par[2])
        # ion and pos
        geo.pos_units = 'angstrom'
        for _section in ["PRIMCOORD\n", "CONVCOORD\n"]:
            out += _section
            out += "    {:5d}    {:3d}\n".format(geo.nat, 1)
            for ion, pos in zip(geo.ion, geo.pos):
                out += "    {:5s}  {:16.9f}  {:16.9f}  {:16.9f}\n".format(ion, pos[0], pos[1], pos[2])

    else:
        raise ValueError('Value of ftype not recognized: {}'.format(ftype))

    return out


def __str__(self, ftype='qeinp'):
    '''
    default str method of atomgeo
    '''
    return self.to_string(ftype=ftype)


def __repr__(self):
    '''
    define atomgeo representation (str)
    '''
    out = "atomgeo("
    out += "ion = " + repr(self.ion) + ","
    out += "par = " + repr(self.par) + ","
    out += "pos = " + repr(self.pos) + ","
    out += "par_units = " + repr(self.par_units) + ","
    out += "pos_units = " + repr(self.pos_units) + ")"
    return out
