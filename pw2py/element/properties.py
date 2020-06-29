# define mass property
@property
def mass(self):
    ''' return value of element._mass '''
    return self._mass


@mass.setter
def mass(self, arg):
    ''' set value of element._mass '''
    self._mass = arg
    self._number = self.request(arg, 'number')
    self._symbol = self.request(arg, 'symbol')


# define number property
@property
def number(self):
    ''' return value of element._number '''
    return self._number


@number.setter
def number(self, arg):
    ''' set value of element._number '''
    self._number = arg
    self._mass = self.request(arg, 'mass')
    self._symbol = self.request(arg, 'symbol')


# define symbol property
@property
def symbol(self):
    ''' return value of element._symbol '''
    return self._symbol


@symbol.setter
def symbol(self, arg):
    ''' set value of element._symbol '''
    self._symbol = arg
    self._mass = self.request(arg, 'mass')
    self._number = self.request(arg, 'number')
