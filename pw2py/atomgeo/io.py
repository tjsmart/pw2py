from ..functions import read_geo


@classmethod
def from_file(cls, filename, ftype='auto', xyz_par=20):
    '''
    generate atomgeo object from file

    input
    ---
    filename (str)
        - name/path to file to be read

    ftype (str)
        - specify filetype otherwise ftype='auto' will use filename to detect file type
        - acceptable values: ['auto', 'vasp', 'qeinp', 'qeout', 'jdftx', 'xsf', 'xyz']

    xyz_par (float)
        - optional, only used when reading from xyz, par is set to a box of with dimension xyz_par in angstrom

    returns
    ---
    atomgeo (atomgeo object)

    notes
    ---
    reads geometry file (filename) using the ..functions.read_geo method
    '''
    ion, par, par_units, pos, pos_units = read_geo(filename)

    return cls(ion=ion, par=par, pos=pos, pos_units=pos_units, par_units=par_units)  # pylint: disable=E1102


def write_file(self, filename, ftype='auto'):
    '''
    write atomgeo to file (QE format)
    '''
    if ftype == 'auto':
        if any([filename.lower().endswith(ext.lower()) for ext in ["OUTCAR", "CONTCAR", "vasp"]]):
            ftype = 'vasp'
        elif filename.lower().endswith("xyz"):
            ftype = 'xyz'
        elif filename.lower().endswith("xsf"):
            ftype = 'xsf'
        elif any([filename.lower().endswith(ext.lower()) for ext in ["pos", "jdftx"]]):
            ftype = 'jdftx'
        else:
            ftype = 'qeinp'
    elif ftype not in ['vasp', 'xyz', 'xsf', 'jdftx', 'qeinp']:
        raise ValueError('Value of ftype not recognized')

    with open(filename, 'w') as f:
        f.write(self.__str__(ftype=ftype))

    return None
