import numpy as np
from warnings import warn

from .. import atomgeo
from .. import element
from .._common.constants import bohr_to_angstrom
from .._common.pseudos import load_pseudo_dict, determine_pseudo_type


def _pseudos_from_ATOMIC_SPECIES(self):
    '''
    determine type of pseudopotentials being used (pulling from first ions file name)
    '''
    pseudo_name = list(self.card.ATOMIC_SPECIES.values())[0][1]
    library, xc_type = determine_pseudo_type(pseudo_name)
    # load dictionary of pseudopotential file names, filter down to needed library/xc_type
    if library == 'dojo':
        pseudos = load_pseudo_dict()[library]
    else:
        pseudos = load_pseudo_dict()[library][xc_type]

    return pseudos


def _fix_species(self):
    '''
    for methods which edit ion or edit atomic species this method will fix
        them so they are consistent
    '''
    # store attribute ion and store ion from atomic species
    attr_ion = np.unique(self.ion)
    card_ion = list(self.card.ATOMIC_SPECIES.keys())

    # first add ion missing from atomic species
    pseudos = _pseudos_from_ATOMIC_SPECIES(self)
    for ion in attr_ion:
        if ion not in card_ion:
            symbol = ''.join(filter(str.isalpha, ion))
            self.card['ATOMIC_SPECIES'][ion] = [
                element.request(symbol, 'mass'), pseudos[symbol]]

    # second remove ion in atomic species not in ion
    for ion in card_ion:
        if ion not in attr_ion:
            self.card['ATOMIC_SPECIES'].pop(ion)


# TODO
# def convert_ibrav(self, newBrav):
#     '''
#     convert ibrav of self to newBrav (intended for ibrav != 0 to 0 or back)
#     '''
#     pass


def load_geo(self, geo, load='default', keep_if_pos=False):
    '''
    return new qeinp object with geometry from geo

    load (str)
        -- can be default, cell, or atoms
    keep_if_pos (bool)
        -- if true, and geo does not have attribute 'if_pos' then keep if_pos of self
    '''
    load = load.lower()
    if load not in ['default', 'cell', 'atoms']:
        raise ValueError('Unrecognized value for load: {}'.format(load))

    if not isinstance(geo, atomgeo):
        raise ValueError(
            'geo must be of type pw2py.atomgeo, passed type: {}'.format(type(geo)))

    if load == 'default' or load == 'atoms':
        self._pos = geo.pos
        self._pos_units = geo.pos_units
        self.card['ATOMIC_POSITIONS'] = geo.pos_units
        try:
            self._if_pos = geo.if_pos
        except AttributeError:
            if not keep_if_pos:
                # use default value for if_pos
                self._if_pos = np.ones((geo.nat, 3))

        if not np.array_equal(self._ion, geo.ion):
            self._ion = geo.ion
            # fix atomic species
            _fix_species(self)

        # resync nat and ntyp in nml
        self.nml['system']['nat'] = self.nat
        self.nml['system']['ntyp'] = self.ntyp

    if load == 'default' or load == 'cell':
        # TODO move this into an implementation of setting ibrav from 0 (i.e. put ibrav setter)
        # Try to preserve value of ibrav
        if self.ibrav in [1, 2, 3]:
            # check if new par == old par * constant
            if np.allclose(geo.par/geo.par[0, 0], self.par/self.par[0, 0]):
                if 'a' in self.nml['system']:
                    self.nml['system']['a'] = geo.par[0, 0]
                    if geo.par_units == 'bohr':
                        self.nml['system']['a'] *= bohr_to_angstrom
                elif 'celldm' in self.nml['system']:
                    self.nml['system']['celldm'][0] = geo.par[0, 0]
                    if geo.par_units == 'angstrom':
                        self.nml['system']['celldm'][0] /= bohr_to_angstrom
                else:
                    warn('Unable to preserve ibrav due to missing nml items')
                    self.ibrav = 0
            else:
                self.ibrav = 0

        elif self.ibrav in [4, 6]:
            # check if par[:2] = old par[:2] * constant and par[2] = old par[2] * constant
            if np.allclose(geo.par[:2]/geo.par[0, 0], self.par[:2]/self.par[0, 0]) and \
                    np.allclose(geo.par[2]/geo.par[2, 2], self.par[2]/self.par[2, 2]):
                if 'a' in self.nml['system']:
                    self.nml['system']['a'] = geo.par[0, 0]
                    if geo.par_units == 'bohr':
                        self.nml['system']['a'] *= bohr_to_angstrom
                    self.nml['system']['c'] = geo.par[2, 2]
                    if geo.par_units == 'bohr':
                        self.nml['system']['c'] *= bohr_to_angstrom
                elif 'celldm' in self.nml['system']:
                    self.nml['system']['celldm'][0] = geo.par[0, 0]
                    if geo.par_units == 'angstrom':
                        self.nml['system']['celldm'][0] /= bohr_to_angstrom
                    self.nml['system']['celldm'][2] = geo.par[2, 2] / \
                        geo.par[0, 0]
                else:
                    warn('Unable to preserve ibrav due to missing nml items')
                    self.ibrav = 0
            else:
                self.ibrav = 0

        elif self.ibrav == 8:
            # check if par[i] = old par[i] for i in range(3)
            if all([np.allclose(geo.par[i]/geo.par[i, i], self.par[i]/self.par[i, i]) for i in range(3)]):
                if 'a' in self.nml['system']:
                    self.nml['system']['a'] = geo.par[0, 0]
                    if geo.par_units == 'bohr':
                        self.nml['system']['a'] *= bohr_to_angstrom
                    self.nml['system']['b'] = geo.par[1, 1]
                    if geo.par_units == 'bohr':
                        self.nml['system']['b'] *= bohr_to_angstrom
                    self.nml['system']['c'] = geo.par[2, 2]
                    if geo.par_units == 'bohr':
                        self.nml['system']['c'] *= bohr_to_angstrom
                elif 'celldm' in self.nml['system']:
                    self.nml['system']['celldm'][0] = geo.par[0, 0]
                    if geo.par_units == 'angstrom':
                        self.nml['system']['celldm'][0] /= bohr_to_angstrom
                    self.nml['system']['celldm'][1] = geo.par[1, 1] / \
                        geo.par[0, 0]
                    self.nml['system']['celldm'][2] = geo.par[2, 2] / \
                        geo.par[0, 0]
                else:
                    warn('Unable to preserve ibrav due to missing nml items')
                    self.ibrav = 0
            else:
                self.ibrav = 0

        elif self.ibrav != 0:
            warn('Value of ibrav not preserved because it is not implemented in load_pos: {}'.format(
                self.ibrav))
            self.ibrav = 0

        # now that the nasty part here is what you'd expect
        # update par, par_units,
        self._par = geo.par
        self._par_units = geo.par_units
        self.card['CELL_PARAMETERS'] = geo.par_units


def replace_ion(self, old_ion, new_ion):
    '''
    replace ion species 'old_ion' with 'new_ion'

    calls atomgeo.replace_ion(self, old_ion, new_ion)
    '''
    atomgeo.replace_ion(self, old_ion, new_ion)

    # resync ntyp in nml
    self.nml['system']['ntyp'] = self.ntyp

    _fix_species(self)


def add_atom(self, atoms, if_pos=[1, 1, 1]):
    '''
    append atoms to self.ion and self.pos
    atoms = ion, pos

    calls atomgeo.add_atom(self, atoms)
    '''
    atomgeo.add_atom(self, atoms)

    # append if_pos
    init_nat = len(self.if_pos)
    append_if_pos = np.ones((self.nat - init_nat, 3), dtype=int)
    append_if_pos[:] = if_pos
    self._if_pos = np.append(self.if_pos, append_if_pos).reshape(self.nat, 3)

    # resync nat and ntyp in nml
    self.nml['system']['nat'] = self.nat
    self.nml['system']['ntyp'] = self.ntyp

    # fix ATOMIC_SPECIES
    _fix_species(self)


def remove_indices(self, indices):
    '''
    remove list of indices from self.ion and self.pos

    calls atomgeo.remove_indices(self, indices)
    '''
    atomgeo.remove_indices(self, indices)

    # resync nat and ntyp in nml
    self.nml['system']['nat'] = self.nat
    self.nml['system']['ntyp'] = self.ntyp

    # fix ATOMIC_SPECIES
    _fix_species(self)


def build_supercell(self, P, inplace=True, if_pos=[1, 1, 1]):
    '''
    Build supercell

    Only simple supercells are implemented wherein the shape of the cell cannot be changed

    P (np.array, int, size = (3,1))
    '''
    atomgeo.build_supercell(self, P, inplace=inplace)

    # append if_pos
    init_nat = len(self.if_pos)
    append_if_pos = np.ones((self.nat - init_nat, 3), dtype=int)
    append_if_pos[:] = if_pos
    self._if_pos = np.append(self.if_pos, append_if_pos).reshape(self.nat, 3)

    # resync nat and ntyp in nml
    self.nml['system']['nat'] = self.nat
    self.nml['system']['ntyp'] = self.ntyp

    # fix ATOMIC_SPECIES
    _fix_species(self)


def to_atomgeo(self):
    '''
    Return atomgeo object
    '''
    return atomgeo(
        ion=self.ion,
        par=self.par,
        pos=self.pos,
        par_units=self.par_units,
        pos_units=self.pos_units
    )
