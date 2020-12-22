
import numpy as np


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
