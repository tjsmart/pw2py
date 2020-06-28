def __init__(self, arg):
    ''' initialize element instance '''
    if isinstance(arg, int):
        # then arg is atomic number
        self.number = arg
    elif isinstance(arg, float):
        # then arg is mass
        self.mass = arg
    elif isinstance(arg, str):
        # then arg is symbol
        self.symbol = arg
    else:
        raise TypeError("Can only instatiate element from int, float, or str, passed type: {}".format(type(arg)))

