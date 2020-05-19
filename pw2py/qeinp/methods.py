from mendeleev import element
import numpy as np
from warnings import warn

from .. import atomgeo, qegeo
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
            self.card._card['ATOMIC_SPECIES'][ion] = [element(ion).mass, pseudos[ion]]

    # second remove ion in atomic species not in ion
    for ion in card_ion:
        if ion not in attr_ion:
            self.card._card['ATOMIC_SPECIES'].pop(ion)


# TODO
# def convert_ibrav(self, newBrav):
#     '''
#     convert ibrav of self to newBrav (intended for ibrav != 0 to 0 or back)
#     '''
#     pass


def load_geo(self, geo, load='default'):
    '''
    return new qeinp object with geometry from geo

    load (str)
        -- can be default, cell, or atoms
    '''
    load = load.lower()
    if load not in ['default', 'cell', 'atoms']:
        raise ValueError('Unrecognized value for load: {}'.format(load))

    if isinstance(geo, atomgeo):
        # convert atomgeo to qegeo
        geo = qegeo(geo, None, None)
    elif not isinstance(geo, qegeo):
        raise ValueError('geo must be of type atomgeo or qegeo, passed type: {}'.format(type(geo)))

    if load == 'default' or load == 'atoms':
        self._pos = geo.pos
        self._pos_units = geo.pos_units
        self.card._card['ATOMIC_POSITIONS'] = geo.pos_units
        self._if_pos = geo.if_pos

        if not np.array_equal(self._ion, geo.ion):
            self._ion = geo.ion
            # fix atomic species
            _fix_species(self)

        # resync nat and ntyp in nml
        self.nml._nml['system']['nat'] = self.nat
        self.nml._nml['system']['ntyp'] = self.ntyp

    if load == 'default' or load == 'cell':
        # TODO move this into an implementation of setting ibrav from 0 (i.e. put ibrav setter)
        # Try to preserve value of ibrav
        if self.ibrav in [1, 2, 3]:
            # check if new par == old par * constant
            if np.allclose(geo.par/geo.par[0, 0], self.par/self.par[0, 0]):
                if 'A' in self.nml._nml['system']:
                    self.nml._nml['system']['A'] = geo.par[0, 0]
                    if geo.par_units == 'bohr':
                        self.nml._nml['system']['A'] *= bohr_to_angstrom
                elif 'celldm' in self.nml._nml['system']:
                    self.nml._nml['system']['celldm'][0] = geo.par[0, 0]
                    if geo.par_units == 'angstrom':
                        self.nml._nml['system']['celldm'][0] /= bohr_to_angstrom
                else:
                    warn('Unable to preserve ibrav due to missing nml items')
                    self.ibrav = 0
            else:
                self.ibrav = 0

        elif self.ibrav in [4, 6]:
            # check if par[:2] = old par[:2] * constant and par[2] = old par[2] * constant
            if np.allclose(geo.par[:2]/geo.par[0, 0], self.par[:2]/self.par[0, 0]) and \
                    np.allclose(geo.par[2]/geo.par[2, 2], self.par[2]/self.par[2, 2]):
                if 'A' in self.nml._nml['system']:
                    self.nml._nml['system']['A'] = geo.par[0, 0]
                    if geo.par_units == 'bohr':
                        self.nml._nml['system']['A'] *= bohr_to_angstrom
                    self.nml._nml['system']['C'] = geo.par[2, 2]
                    if geo.par_units == 'bohr':
                        self.nml._nml['system']['C'] *= bohr_to_angstrom
                elif 'celldm' in self.nml._nml['system']:
                    self.nml._nml['system']['celldm'][0] = geo.par[0, 0]
                    if geo.par_units == 'angstrom':
                        self.nml._nml['system']['celldm'][0] /= bohr_to_angstrom
                    self.nml._nml['system']['celldm'][2] = geo.par[2, 2] / geo.par[0, 0]
                else:
                    warn('Unable to preserve ibrav due to missing nml items')
                    self.ibrav = 0
            else:
                self.ibrav = 0

        elif self.ibrav == 8:
            # check if par[i] = old par[i] for i in range(3)
            if all([np.allclose(geo.par[i]/geo.par[i, i], self.par[i]/self.par[i, i]) for i in range(3)]):
                if 'A' in self.nml._nml['system']:
                    self.nml._nml['system']['A'] = geo.par[0, 0]
                    if geo.par_units == 'bohr':
                        self.nml._nml['system']['A'] *= bohr_to_angstrom
                    self.nml._nml['system']['B'] = geo.par[1, 1]
                    if geo.par_units == 'bohr':
                        self.nml._nml['system']['B'] *= bohr_to_angstrom
                    self.nml._nml['system']['C'] = geo.par[2, 2]
                    if geo.par_units == 'bohr':
                        self.nml._nml['system']['C'] *= bohr_to_angstrom
                elif 'celldm' in self.nml._nml['system']:
                    self.nml._nml['system']['celldm'][0] = geo.par[0, 0]
                    if geo.par_units == 'angstrom':
                        self.nml._nml['system']['celldm'][0] /= bohr_to_angstrom
                    self.nml._nml['system']['celldm'][1] = geo.par[1, 1] / geo.par[0, 0]
                    self.nml._nml['system']['celldm'][2] = geo.par[2, 2] / geo.par[0, 0]
                else:
                    warn('Unable to preserve ibrav due to missing nml items')
                    self.ibrav = 0
            else:
                self.ibrav = 0

        elif self.ibrav != 0:
            warn('Value of ibrav not preserved because it is not implemented in load_pos: {}'.format(self.ibrav))
            self.ibrav = 0

        # now that the nasty part here is what you'd expect
        # update par, par_units,
        self._par = geo.par
        self._par_units = geo.par_units
        self.card._card['CELL_PARAMETERS'] = geo.par_units


def replace_ion(self, old_ion, new_ion):
    '''
    replace ion species 'old_ion' with 'new_ion'

    calls atomgeo.replace_ion(self, old_ion, new_ion)
    '''
    atomgeo.replace_ion(self, old_ion, new_ion)
    _fix_species(self)


def add_atom(self, atoms):
    '''
    append atoms to self.ion and self.pos
    atoms = ion, pos

    calls atomgeo.add_atom(self, atoms)
    '''
    atomgeo.add_atom(self, atoms)

    # resync nat and ntyp in nml
    self.nml._nml['system']['nat'] = self.nat
    self.nml._nml['system']['ntyp'] = self.ntyp

    # fix ATOMIC_SPECIES
    _fix_species(self)


def remove_indices(self, indices):
    '''
    remove list of indices from self.ion and self.pos

    calls atomgeo.remove_indices(self, indices)
    '''
    atomgeo.remove_indices(self, indices)

    # resync nat and ntyp in nml
    self.nml._nml['system']['nat'] = self.nat
    self.nml._nml['system']['ntyp'] = self.ntyp

    # fix ATOMIC_SPECIES
    _fix_species(self)
