class bands:
    '''
    class for bands

    typical usage is to create an instance using the 'from_file' methods:

        bands = pw.bands.from_file(filename)
    '''

    from .io import from_file, from_save
    from .properties import eig, kpt, fermi, occ, is_metallic, vbm, cbm, gap, \
        nspin, nk, nbnd, weight
    from .init import __init__, _check_is_metallic, _calc_bandedge
    # from .str import __repr__
    from .methods import shift_bands
