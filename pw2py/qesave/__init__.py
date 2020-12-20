'''
module housing several functions for reading qe *.save folder. More details to come!

Example usage of loaders for qe 6.6:

    >>> from pw2py.qesave.load66 import read_wavefunction
    >>> gk, evc = read_wavefunction('./path/to/savefolder')

'''

from . import load61  # noqa: F401
from . import load66  # noqa: F401


class qesave:
    from .init import __init__
