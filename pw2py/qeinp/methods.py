from .. import atomgeo, qegeo


def load_geo(self, geo, load='all'):
    '''
    return new qeinp object with geometry from geo

    load (str)
        -- can be all, par, or atoms
    '''
    load = load.lower()
    if load not in ['all', 'par', 'atoms']:
        raise ValueError('Unrecognized value for load: {}'.format(load))

    if isinstance(geo, atomgeo):
        # convert atomgeo to qegeo
        geo = qegeo(geo, None, None)
    elif not isinstance(geo, qegeo):
        raise ValueError('geo must be of type atomgeo or qegeo, passed type: {}'.format(type(geo)))

    if load == 'all' or load == 'atoms':
        self._ion = geo.ion
        self._pos = geo.pos
        self._pos_units = geo.pos_units
        self._if_pos = geo.if_pos
    if load == 'all' or load == 'par':
        self._par = geo.par
        self._par_units = geo.par_units
        self.ibrav = 0
