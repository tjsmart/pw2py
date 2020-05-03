from copy import deepcopy
from numpy import array, append, divide, multiply, place
from pandas import DataFrame


def replace_ion(self, old_ion, new_ion):
    '''
    replace ion species 'old_ion' with 'new_ion'
    '''
    place(self.ion, self.ion == old_ion, new_ion)

    # TODO for qeinp need to edit ATOMIC_SPECIES


def add_atom(self, atoms):
    '''
    append atoms to self.ion and self.pos
    atoms = ion, pos
    '''
    ion, pos = atoms
    ion = array(ion, dtype=str)
    ion = ion.flatten()
    pos = array(pos, dtype=float)
    pos = pos.reshape((ion.size, 3))
    self._ion = append(self._ion, ion)
    self._pos = append(self._pos, pos).reshape((self.nat, 3))

    # TODO for qeinp need to increment nat, ntyp (potentially), ATOMIC_SPECIES


def remove_indices(self, indices):
    '''
    remove list of indices from self.ion and self.pos
    '''
    filt = [index not in indices for index in range(self.nat)]
    self._ion = self._ion[filt]
    self._pos = self._pos[filt]

    # TODO for qeinp need to update nat, ntyp, ATOMIC_SPECIES


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
    out.pos = array(df.loc[:, 'pos0':'pos2'])

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

    P (array, int, size = (3,1))
    '''
    # format input transformation vector
    P = array(P, dtype=int).reshape((3, 1))
    # deepcopy if not inplace
    out = self if inplace else deepcopy(self)
    # store pos units
    _save_units = deepcopy(self.pos_units)
    # convert to cystal and shift all to fit in unit
    # (in principle this can be skipped but it is anticipated to necessary for complicated transfroms)
    out.pos_units = 'crystal'
    out.shift_pos_to_unit()

    # rescale lattice parameters
    out._par = multiply(out._par, P.reshape(-1, 1))

    # rescale atomic positions
    out._pos = divide(out._pos, P.reshape(1, -1))

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
    out._ion = append(out._ion, ion_images)
    out._pos = append(out._pos, pos_images).reshape((out.nat, 3))

    # restore original pos_units
    out.pos_units = _save_units

    if not inplace:
        return out
