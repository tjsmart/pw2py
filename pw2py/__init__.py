'''
pw2py
=====
Provides
  1. An interface to between several atomic geometry formats via a python class
  2. Several routines to manipulate atomic geometry
  3. Interface with quantum espresso input and output
Documentation
----------------------------
Documentation is available at `the pw2py git repository <https://gitlab.com/tjsmart/pw2py>`.
Usage
----------------------------
Import `pw2py` as `pw`::
  >>> import pw2py as pw
Create class object from file (`atomgeo`)::
  >>> geo = atomgeo.from_file('filename')
Available subpackages
---------------------
atomgeo
    class for reading/manipulating/writing atomic geometry
qegeo
    subclass of atomgeo with extended functionality for quantum espresso
qeinp
    subclass of qegeo with functionality for handling full quantum espresso input file
qeout
    class for processing quantum espresso output file
qesave
    module with functions for reading quantum espresso *.save folder
'''

from .atomgeo import atomgeo  # noqa: F401
from .qegeo import qegeo  # noqa: F401
from .qeinp import qeinp  # noqa: F401
from .qeout import qeout  # noqa: F401
from . import qesave  # noqa: F401
