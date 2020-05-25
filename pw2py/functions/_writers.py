#  functions for writing different file types

from copy import deepcopy
import numpy as np


def write_xsf(geo, filename, force=None, grid=None):
    '''
    Write geo to xsf file
    '''
    # check inputs
    if force is not None:
        # cast force to the appropriate shape
        force = np.array(force).reshape((geo.nat, 3))
    if grid is not None:
        # TODO
        raise NotImplementedError('Grid is not yet supported, coming soon!')

    # check units, if necessary make copy
    if not geo.pos_units == 'angstrom' or not geo.par_units == 'angstrom':
        geo_c = deepcopy(geo)
        if not geo.pos_units == 'angstrom':
            geo_c.pos_units = 'angstrom'
        if not geo.par_units == 'angstrom':
            geo_c.par_units = 'angstrom'
    else:
        geo_c = geo

    with open(filename, 'w') as f:
        # start building output string
        f.write("# Generated from the PW2PY python module\n")
        f.write("DIM-GROUP\n")
        f.write("    {:5d}    {:3d}\n".format(3, 1))
        for _section in ["PRIMVEC\n", "CONVVEC\n"]:
            f.write(_section)
            for par in geo_c.par:
                f.write("    {:16.9f}  {:16.9f}  {:16.9f}\n".format(par[0], par[1], par[2]))
        # ion and pos (possibly force)
        for _section in ["PRIMCOORD\n", "CONVCOORD\n"]:
            f.write(_section)
            f.write("    {:5d}    {:3d}\n".format(geo_c.nat, 1))
            if force is None:
                for ion, pos in zip(geo_c.ion, geo_c.pos):
                    f.write("    {:5s}  {:16.9f}  {:16.9f}  {:16.9f}\n".format(ion, pos[0], pos[1], pos[2]))
            else:
                for ion, pos, frc in zip(geo_c.ion, geo_c.pos, force):
                    f.write("    {:5s}  {:16.9f}  {:16.9f}  {:16.9f}  {:16.9f}  {:16.9f}  {:16.9f}\n".format(
                            ion, pos[0], pos[1], pos[2], frc[0], frc[1], frc[2]
                            ))
