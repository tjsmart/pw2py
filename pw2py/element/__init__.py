class element:
    '''
    class for elements

    can make instatiate an element by its mass (float), symbol (string), or atomic number (int)

    Attributes:
    -----

    mass (float)
        - atomic mass in a.u.

    symbol (string)
        - atomic symbol (e.g. 'H', 'O', ...)

    number (int)
        - atomic number
    '''

    # properties
    from .properties import mass, symbol, number
    # dunder methods
    from .init import __init__
    from .str import __repr__
    # other methods
    from .methods import request, symbols
    from .element_db import load_db

