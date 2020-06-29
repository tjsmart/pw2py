from .element_db import load_db


@staticmethod
def request(arg, attr, inp='auto'):
    '''
    request elemental attribute (mass, symbol, or number) from another provided attribute (arg)

    input
    ----

        arg
            - can be float (mass), str (symbol), or int (number)
        attr (str)
            - name of requested attribute, can be 'mass', 'symbol', or 'number'
        inp (str, optional)
            - specify input type ('mass', 'symbol', or 'number'), default='auto'

    return
    ---
        mass(float), symbol(str), or number(int)
    '''
    assert attr in ['number', 'symbol', 'mass'], "attr must be 'number', 'symbol' or 'mass', passed: {}".format(attr)
    if inp == 'auto':
        # here arg must be int, str, or float
        if isinstance(arg, int):
            inp = 'number'
        elif isinstance(arg, str):
            inp = 'symbol'
        elif isinstance(arg, float):
            inp = 'mass'
        else:
            raise TypeError("Input argument must be int, float, or str, passed type: {}".format(type(arg)))
    else:
        # here inp must be number, symbol, or mass
        if inp == 'number':
            arg = int(arg)
        elif inp == 'symbol':
            arg = str(arg)
        elif inp == 'mass':
            arg = float(arg)
        else:
            raise TypeError("inp must be 'number', 'symbol' or 'mass', passed: {}".format(inp))
    # return request
    if attr == inp:
        return arg
    if inp == 'mass':
        arg_str = '{:.4f}'.format(arg)
    elif inp == 'number':
        arg_str = str(arg)
    elif inp == 'symbol':
        arg_str = arg
    db_name = inp + '_to_' + attr
    db = load_db(db_name)
    return db[arg_str]


@staticmethod
def symbols():
    '''
    returns a list of all atomic symbols
    '''
    return load_db('symbols')
