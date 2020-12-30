import numpy as np

from ..functions import read_geo
from ..qesave import qesave


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


def write_file(self, filename, ftype='auto', sort=True):
    '''
    write atomgeo to file (QE format)
    '''
    if ftype == 'auto':
        if any([filename.lower().endswith(ext.lower()) for ext in ["POSCAR", "CONTCAR", "vasp"]]):
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
        f.write(self.__str__(ftype=ftype, sort=sort))

    return None


def _read_qesave_geo_61(save):
    '''
    Read atomic geometry from 6.1 save folder
    '''
    # ------------------------------------------------------------
    # read cell parameters
    latt_child = save.root.find('CELL').find('DIRECT_LATTICE_VECTORS')
    par_units = latt_child.find(
        'UNITS_FOR_DIRECT_LATTICE_VECTORS').attrib['UNITS'].lower()
    par = [latt_child.find(f'a{i+1}').text.split() for i in range(3)]
    # ------------------------------------------------------------
    # read ions and positions
    ion_child = save.root.find('IONS')
    nat = int(ion_child.find('NUMBER_OF_ATOMS').text)
    pos_units = ion_child.find(
        'UNITS_FOR_ATOMIC_POSITIONS').attrib['UNITS'].lower()
    ion, pos = [], []
    for i in range(nat):
        atom_attrib = ion_child.find(f'ATOM.{i+1}').attrib
        ion.append(atom_attrib['SPECIES'])
        pos.append(atom_attrib['tau'].split())
    # ------------------------------------------------------------
    # return dictionary
    return dict(ion=ion, par=par, pos=pos, pos_units=pos_units, par_units=par_units)


def _read_qesave_geo_66(save):
    '''
    Read atomic geometry from 6.6 save folder
    '''
    # ------------------------------------------------------------
    # all units are hartree
    par_units = 'bohr'
    pos_units = 'bohr'
    # ------------------------------------------------------------
    # read ion and positions
    as_child = save.root.find('input').find('atomic_structure')
    ion, pos = [], []
    for a_child in as_child.find('atomic_positions').findall('atom'):
        ion.append(a_child.attrib['name'])
        pos.append(a_child.text.split())
    # ------------------------------------------------------------
    # read cell parameters
    cell_child = as_child.find('cell')
    par = [cell_child.find(f'a{i+1}').text.split() for i in range(3)]
    # ------------------------------------------------------------
    # return dictionary
    return dict(ion=ion, par=par, pos=pos, pos_units=pos_units, par_units=par_units)


@classmethod
def from_save(cls, prefix: qesave, path: str = '.', version: str = None):
    '''
    generate atomgeo object qe savefolder
    '''
    save = qesave(prefix=prefix, path=path)
    version = save.version if version is None else version
    if version == '6.1':
        return cls(**_read_qesave_geo_61(save))
    else:
        return cls(**_read_qesave_geo_66(save))
