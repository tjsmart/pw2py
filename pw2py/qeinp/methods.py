from . import qeinp
from .. import atomgeo, qegeo


def load_geo(self, geo):
    '''
    return new qeinp object with geometry from geo
    '''
    if isinstance(geo, atomgeo):
        # convert to qegeo
        geo = qegeo(geo, None, None)
    elif not isinstance(geo, qegeo):
        raise ValueError('geo must be of type atomgeo or qegeo, passed type: {}'.format(type(geo)))

    return qeinp(self.nml, geo, self.card, self.kpt)
