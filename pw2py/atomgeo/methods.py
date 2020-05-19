from copy import deepcopy
import numpy as np
from pandas import DataFrame


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


def remove_indices(self, indices):
    '''
    remove list of indices from self.ion and self.pos
    '''
    filt = [index not in indices for index in range(self.nat)]
    self._ion = self._ion[filt]
    self._pos = self._pos[filt]


def sort_ions(self, inplace=False):
    '''
    sort ions into categories and reorder corresponding positions
    '''
    out = self if inplace else deepcopy(self)

    # TODO do this without making a dataframe may be better?
    # save ion and pos to dataframe
    columns = ['ion'] + ['pos{}'.format(i) for i in range(3)]
    df = DataFrame(columns=columns)
    df.ion = self.ion
    for i in range(3):
        df['pos{}'.format(i)] = self.pos[:, i]

    # sort dataframe by ion column
    df.sort_values(by=['ion'], inplace=True)

    # update ion and pos
    out.ion = list(df['ion'])
    out.pos = np.array(df.loc[:, 'pos0':'pos2'])

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


def nearest_neighbor(self, site_id, N=1, include_site=False, return_type='index'):
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

    returns
    ----
        see return_type above
    '''
    if not isinstance(N, int):
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
    atoms_df['dist'] = atoms_df['pos'].apply(lambda x: np.linalg.norm(x - atoms_df.loc[site_id, 'pos']))
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
