class qeout:
    '''
    class for quantum espresso output file

    Suggested to make object from file:
    -----

    out = qeout.from_file(filename)

    Attributes:
    -----

    geo (pw2py.atomgeo)
        - atomic geometry (will read *.in if found)
    bands (pw2py.bands)
        - kohn-sham bands
    final_energy (tuple)
        - tuple, of energy and convergence level
    '''

    from .init import __init__

    from .io import from_file

    from .methods import final_energy
