#  functions for writing different file types

from copy import deepcopy
import numpy as np


def write_reshaped_grid(f, grid, grid_order='C', fmt='%13.6e', columns=6):
    '''
    reshape grid and write for xsf file
    '''
    # TODO
    if grid_order != 'C':
        raise NotImplementedError('F order not implemented')
    # reorder grid and write to f
    icol = 0
    for k in range(grid.shape[2] + 1):
        for j in range(grid.shape[1] + 1):
            for i in range(grid.shape[0] + 1):
                try:
                    f.write(fmt % (grid[i, j, k]))
                except IndexError:
                    i %= grid.shape[0]
                    j %= grid.shape[1]
                    k %= grid.shape[2]
                    f.write(fmt % (grid[i, j, k]))
                icol += 1
                if icol == 6:
                    f.write('\n')
                    icol = 0
                else:
                    f.write(' ')


def write_xsf(geo, filename, force=None, grid=None, grid_order='C'):
    '''
    Write geo to xsf file
    '''
    # check inputs
    if force is not None:
        # cast force to the appropriate shape
        force = np.array(force).reshape((geo.nat, 3))
    if grid is not None:
        # make sure grid is an array of 3 dimensions
        grid = np.array(grid)
        assert grid.ndim == 3, 'grid must be an array with 3 dimensions: {}'.format(grid.ndim)

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
        f.write("CRYSTAL\n")
        f.write("PRIMVEC\n")
        for par in geo_c.par:
            f.write("    {:16.9f}  {:16.9f}  {:16.9f}\n".format(*par))
        # ion and pos (possibly force)
        f.write('PRIMCOORD\n')
        f.write("    {:5d}    {:3d}\n".format(geo_c.nat, 1))
        if force is None:
            for ion, pos in zip(geo_c.ion, geo_c.pos):
                f.write("    {:5s}  {:16.9f}  {:16.9f}  {:16.9f}\n".format(ion, *pos))
        else:
            for ion, pos, frc in zip(geo_c.ion, geo_c.pos, force):
                f.write("    {:5s}  {:16.9f}  {:16.9f}  {:16.9f}  {:16.9f}  {:16.9f}  {:16.9f}\n".
                        format(ion, *pos, *frc))
        # write grid (if it exists)
        if grid is not None:
            # grid header
            f.write('BEGIN_BLOCK_DATAGRID_3D\n3D_PWSCF\nDATAGRID_3D_UNKNOWN\n')
            f.write('    {:5d}  {:5d}  {:5d}\n'.format(*(np.array(grid.shape) + 1)))
            f.write('    {:6.4f}  {:6.4f}  {:6.4f}\n'.format(0, 0, 0))  # what is this line for?
            for par in geo_c.par:
                f.write('    {:16.9f}  {:16.9f}  {:16.9f}\n'.format(*par))
            # grid data
            write_reshaped_grid(f, grid, grid_order=grid_order)
            f.write('END_DATAGRID_3D\nEND_BLOCK_DATAGRID_3D\n')
