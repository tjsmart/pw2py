'''
module housing several functions for reading qe *.save folder. More details to come!

Example usage of loaders for qe 6.6:

    >>> import pw2py as pw
    >>> gk, evc = pw.qesave.loaders_66.read_wfc_file('wfcup1.dat')

'''

# from .loaders import (  # noqa : F401
#     read_charge_density, read_data_xml, read_eigenvalues, read_gkvectors,
#     read_gvectors, read_spin_polarization, read_wavefunction
# )

from . import loaders_61  # noqa: F401
from . import loaders_66  # noqa: F401


'''
class qesave:
    from .init import __init__
'''
