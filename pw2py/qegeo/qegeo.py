from .. import atomgeo


class qegeo(atomgeo):
    '''
    same as atomgeo with if_pos and ibrav capabilities
    '''

    from .init import __init__
    from .str import __repr__, __str__
    from .properties import A, B, C, celldm, ibrav, if_pos
    from .io import from_file, write_file
