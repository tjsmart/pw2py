from copy import deepcopy
import numpy as np
from pandas import DataFrame
from warnings import warn
import itertools as it

from .. import element
from ..functions._writers import write_xsf


def _cell_translations(self) -> np.ndarray:
    '''
    cell translations by +1/-1 of cell parameters
    '''
    return np.dot([t for t in it.product([-1, 0, 1], repeat=3)], self.par)


def elements(self):
    '''
    return elemental names of each ion (stripping number i.e. 'Fe2' will return 'Fe')
    '''
    return np.array([''.join(filter(str.isalpha, ion)) for ion in self.ion])


def mass(self, units='au'):
    '''
    return mass of each ion

    TODO implement different unit conversions, currently will only return au
    '''
    units = units.lower()
    if units != 'au':
        raise NotImplementedError('Only au is currently supported')
    mass_dict = element.load_db('symbol_to_mass')
    return np.array([mass_dict[ion] for ion in self.elements()])


def replace_ion(self, old_ion, new_ion):
    '''
    replace ion species 'old_ion' with 'new_ion'
    '''
    np.place(self.ion, self.ion == old_ion, new_ion)


def add_atom(self, atoms):
    '''
    append atoms to self.ion and self.pos
    atoms = ion, pos
    '''
    ion, pos = atoms
    ion = np.array(ion, dtype=str)
    ion = ion.flatten()
    pos = np.array(pos, dtype=float)
    pos = pos.reshape((ion.size, 3))
    self._ion = np.append(self._ion, ion)
    self._pos = np.append(self._pos, pos).reshape((self.nat, 3))


def drop_indices(self, indices, inplace=True):
    '''
    remove list of indices from self.ion and self.pos
    '''
    out = self if inplace else deepcopy(self)
    filt = [index not in indices for index in range(self.nat)]
    out._ion = self._ion[filt]
    out._pos = self._pos[filt]
    if not inplace:
        return out


def sort_ions(self, inplace=False):
    '''
    sort ions into categories and reorder corresponding positions
    '''
    out = self if inplace else deepcopy(self)
    # calculate indices of sorted ion array
    sorted_index = np.argsort(self.ion)
    # update ion and pos
    out.ion = self.ion[sorted_index]
    out.pos = self.pos[sorted_index]
    if not inplace:
        return out


def shift_pos_to_unit(self, inplace=True):
    '''
    Shift all atoms to be within the unit cell dimensions (i.e. between 0<1 in crystal)
    '''
    out = self if inplace else deepcopy(self)
    # store original units, then change to crystal
    _save_units = deepcopy(self.pos_units)
    out.pos_units = 'crystal'
    # loop through positions (px,py,pz), shift to all be in range of 0 <= p < 1
    out.pos = [[p - (p // 1) for p in pos] for pos in out.pos]

    # restore units
    out.pos_units = _save_units
    if not inplace:
        return out


def build_supercell(self, P, inplace=True):
    '''
    Build supercell

    Only simple supercells are implemented wherein the shape of the cell cannot be changed

    P (np.array, int, size = (3,1))
    '''
    # format input transformation vector
    P = np.array(P, dtype=int).reshape((3, 1))
    # deepcopy if not inplace
    out = self if inplace else deepcopy(self)
    # store pos units
    _save_units = deepcopy(self.pos_units)
    # convert to cystal and shift all to fit in unit
    # (in principle this can be skipped but it is anticipated to necessary for complicated transfroms)
    out.pos_units = 'crystal'
    out.shift_pos_to_unit()

    # rescale lattice parameters
    out._par = np.multiply(out._par, P.reshape(-1, 1))

    # rescale atomic positions
    out._pos = np.divide(out._pos, P.reshape(1, -1))

    # generate images of ion and positions
    inv_P = 1 / P.reshape(1, -1)
    ion_images, pos_images = [], []
    for i in range(int(P[0])):
        for j in range(int(P[1])):
            for k in range(int(P[2])):
                if i == j == k == 0:
                    continue
                ion_images.append(out._ion)
                pos_images.append(out._pos + inv_P * [i, j, k])

    # append images
    out._ion = np.append(out._ion, ion_images)
    out._pos = np.append(out._pos, pos_images).reshape((out.nat, 3))

    # restore original pos_units
    out.pos_units = _save_units

    if not inplace:
        return out


def nearest_neighbor(self, site_id, N=1, include_site=False, return_type='index', with_boundaries=True):
    '''
    Calculate nearest neighbors of site_id

    input
    ----
        site_id (int)
            - index (base 0) of the site to calculate nearest neighbors of
        N (optional) (int)
            - number of neighbors to return data for (note if include_site=True,
                total number of elements returned is N+1, otherwise it is N)
        include_site (optional) (bool)
            - if True than the site_id atom is returned as well, default is False
        return_type (optional) (str)
            - possible options include:
                - 'df': return dataframe with columns = ['ion', 'pos', 'dist']
                - 'index': list of indices (ints)
                - 'dist': list of distances in angstrom (floats)
                - 'atoms': list of dictionaries with keys 'ion' and 'pos' (dicts)
        with_boundaries (optional) (bool)
            - if True than boundary conditions are considered when computing distances

    returns
    ----
        see return_type above
    '''
    if not isinstance(site_id, (int, np.integer)):
        raise TypeError(
            'site_id should be an int, not type: {}'.format(type(site_id)))
    if not isinstance(N, (int, np.integer)):
        raise TypeError('N should be an int, not type: {}'.format(type(N)))
    if self.pos_units != 'angstrom':
        # geo is a copy self
        geo = deepcopy(self)
        # convert pos to angstrom
        geo.pos_units = 'angstrom'
    else:
        # geo is a pointer to self
        geo = self
    # build dataframe
    atoms_df = DataFrame(geo.atoms)
    # create column for distances to site_id
    site_pos = atoms_df.loc[site_id, 'pos']
    if not with_boundaries:
        atoms_df['dist'] = atoms_df['pos'].apply(
            lambda x: np.linalg.norm(x - site_pos))
    else:
        # calculate translations by +/-1 of cell parameters
        translations = self._cell_translations()
        # calculate images of the site
        site_images = site_pos + translations
        atoms_df['dist'] = atoms_df['pos'].apply(
            lambda x: np.linalg.norm(x - site_images, axis=1).min()
        )
    # sort by distance
    atoms_df.sort_values(['dist', 'ion'], inplace=True)
    # if including site_id, start_index = 0, else start_index = 1
    start_index = int(not include_site)
    # parse neighbors
    if N == -1:
        atoms_df = atoms_df[start_index:]
    else:
        atoms_df = atoms_df[start_index:N+1]
    # determine answer
    if return_type == 'df':
        # return dataframe
        return atoms_df
    elif return_type == 'index':
        # return list of indices
        answer = list(atoms_df.index)
    elif return_type == 'dist':
        # return distance values as list
        answer = list(atoms_df['dist'])
    elif return_type == 'atoms':
        # transpose so ion and pos are keys, convert to dict, then return just values as list of dicts
        answer = list(atoms_df.loc[['ion', 'pos']].T.to_dict().values())
    else:
        raise ValueError('invalid value for return_type')
    # if len(answer) == 1 return single value, else return list
    if len(answer) == 1:
        return answer[0]
    else:
        return answer


def calc_dR(self, geo, units='angstrom', suppress_warnings=False):
    '''
    Calculate deltaR_i = self.pos - geo.pos

    if necessary creates copies of self and geo to convert them to angstrom
    '''
    units = units.lower()
    # wrong usage checks
    # assert str(type(geo)) == "<class 'pw2py.atomgeo.atomgeo'>", 'geo is not of an instance of pw2py.atomgeo!'
    assert self.nat == geo.nat, 'number of atoms from self and geo do not match: {} != {}'.format(
        self.nat, geo.nat)
    assert units in ['angstrom',
                     'bohr'], 'Only angstrom and bohr are supported'

    # warning checks
    if not suppress_warnings:
        if not np.array_equal(self.par, geo.par):
            warn('Cell parameters do not match!')
        if not np.array_equal(self.elements(), geo.elements()):
            warn('Elements of calculation do not match!')

    # if needed copy and convert to pos to supplied units
    if self.pos_units != units:
        self_c = deepcopy(self)
        self_c.pos_units = units
    else:
        self_c = self

    if geo.pos_units != units:
        geo_c = deepcopy(geo)
        geo_c.pos_units = units
    else:
        geo_c = geo

    return self_c.pos - geo_c.pos


def calc_dR2(self, geo, units='angstrom', suppress_warnings=False):
    '''
    Calculate deltaR2_i = (self.pos - geo.pos)**2

    returns np.sum(np.power(self.calc_dR(geo), 2), axis=1)
    '''
    return np.sum(np.power(self.calc_dR(
        geo, units=units, suppress_warnings=suppress_warnings
    ), 2), axis=1)


def calc_dQ2(self, geo, pos_units='angstrom', mass_units='au', suppress_warnings=False):
    '''
    Calculate deltaQ2 = mass_i * (deltaR_i)**2

    (deltaR_i)**2 = calc_dR2(self, geo)
    '''
    # calculate deltaR_i
    dR2 = calc_dR2(
        self, geo, suppress_warnings=suppress_warnings, units=pos_units)

    # check that self and geo are comprised of the same elements
    assert np.array_equal(self.elements(), geo.elements()), \
        'Elements of calculation do not match! Mass needs to be unambiguous'

    return np.multiply(self.mass(units=mass_units), dR2)


def calc_dQ(self, geo, pos_units='angstrom', mass_units='au', suppress_warnings=False):
    '''
    Calculate deltaQ = sqrt(sum_i (deltaQ_i)**2)

    returns np.sqrt(np.sum(self.calc_dQ2(geo)))
    '''
    return np.sqrt(np.sum(
        self.calc_dQ2(geo, pos_units=pos_units, mass_units=mass_units,
                      suppress_warnings=suppress_warnings)
    ))


def dQ_field_2_xsf(self, geo, filename, suppress_warnings=False):
    '''
    Calculate deltaQ_i = sqrt(mass_i) * (self.pos_i - geo.pos_i)

    Write to xsf file (filename)
    '''
    # calculate deltaQ_field
    deltaQ_field = np.multiply(
        np.power(self.mass().reshape(self.nat, 1), 0.5),
        self.calc_dR(geo, suppress_warnings=suppress_warnings)
    )
    # write to xsf file
    write_xsf(filename, self, force=deltaQ_field)


def calc_distance(self, x: np.ndarray) -> float:
    '''
    Calculate the distance from self.pos and vector x, considering boundary conditions
    '''
    if self.pos_units != 'angstrom':
        # geo is a copy self
        geo = deepcopy(self)
        # convert pos to angstrom
        geo.pos_units = 'angstrom'
    else:
        # geo is a pointer to self
        geo = self
    # check users input X
    X = np.array(x)
    assert X.shape == (3,), \
        "x must be vector in 3D with shape (3,), passed: {}".format(x.shape)
    translations = geo._cell_translations()
    X_images = X + translations
    distances = np.array([
        np.linalg.norm(pos - X_images, axis=1).min() for pos in geo.pos
    ])
    return distances


def composition(self) -> dict:
    '''
    Return the composition of ions
    '''
    return {k: v for k, v in zip(*np.unique(self.ion, return_counts=True))}


def copy(self):
    '''
    Returns a deepcopy
    '''
    return deepcopy(self)
