from copy import deepcopy
from numpy import append, array, divide, multiply


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
