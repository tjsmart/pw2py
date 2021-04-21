
import numpy as np

from ..element import element
# from ..atomgeo import atomgeo


def _skip_lines(f, i):
    [f.readline() for _ in range(i)]


def read_xsf_grid(filename: str) -> np.ndarray:
    '''
    read grid from xsf file
    '''
    with open(filename, 'r') as f:
        for line in f:
            if line.strip().startswith('BEGIN_BLOCK_DATAGRID_3D'):
                _skip_lines(f, 2)
                shape = np.fromstring(f.readline(), sep=' ', dtype=int)
                _skip_lines(f, 4)
                data = []
                while True:
                    line = f.readline().strip()
                    if line.startswith('END'):
                        break
                    data += line.split()
    data = np.array(data, np.float64).reshape(shape, order='F')
    return np.array(data, order='C')


# def read_cub(filename: str) -> (atomgeo, np.ndarray):
# TODO Need to fix this
def read_cub(filename: str):
    '''
    read grid from cub file
    '''
    with open(filename, 'r') as f:
        _skip_lines(f, 2)
        line_dat = np.fromstring(f.readline(), sep=' ', dtype=float)
        nat = int(line_dat[0])
        assert np.allclose(line_dat[1:4], np.zeros(3)), \
            "Sorry, not implemented"
        # ------------------------------------------------------------
        # read fft and cell parameters
        fft = []
        par = []
        for _ in range(3):
            line_dat = np.fromstring(f.readline(), sep=' ', dtype=float)
            fft.append(line_dat[0])
            par.append(line_dat[1:4])
        fft = np.array(fft, dtype=int)
        par = np.array(par)
        par *= fft.reshape((3, 1))
        # ------------------------------------------------------------
        # read atoms (ion name and positions)
        ion = []
        pos = []
        for _ in range(nat):
            line_dat = np.fromstring(f.readline(), sep=' ', dtype=float)
            num = int(line_dat[0])
            ion.append(element(num).symbol)
            pos.append(line_dat[2:5])
    # ------------------------------------------------------------
    # read grid
    skiprows = 6 + nat
    grid = np.loadtxt(filename, skiprows=skiprows)
    grid = grid.reshape(fft)
    # return an atomgeo object and the grid
    # return atomgeo(ion=ion, par=par, par_units='bohr', pos=pos, pos_units='bohr'), grid
    geo = dict(ion=ion, par=par, par_units='bohr', pos=pos, pos_units='bohr')
    return geo, grid
