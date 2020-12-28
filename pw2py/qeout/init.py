def __init__(self, **kw):
    '''
    initialize qeout from keyword arguments

    input
    ---
        kw: dict, with keys:
            bands           : pw2py.bands, kohn-sham bands
            geo             : pw2py.atomgo, atomic geometry
            final_energy    : tuple(float, str), final total energy (in eV) and then convergence level

    returns
    ---
        qeout object
    '''
    if 'bands' in kw:
        self.bands = kw['bands']
    if 'geo' in kw:
        self.geo = kw['geo']
    if 'final_energy' in kw:
        self.final_energy = kw['final_energy']
