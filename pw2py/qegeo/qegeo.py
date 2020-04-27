from .. import atomgeo


class qegeo(atomgeo):
    '''
    same as atomgeo with if_pos and ibrav capabilities
    '''

    from ._init import __init__
    from ._str import __repr__, __str__
    from ._properties import A, B, C, celldm, ibrav, if_pos
    from ._io import from_file, write_file
