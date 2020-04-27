from copy import deepcopy
from numpy import array, append, place
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
    return out
