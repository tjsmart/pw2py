#!/usr/bin/env python3

import os
import mendeleev as mlv
import pickle as pkl


# file where mass will be stored
dump_file = os.path.join(os.path.dirname(__file__), 'mass.pkl')


def generate_mass_dict():
    '''
    creates dictionary relating element symbols (e.g., 'H' or 'O') and element mass (e.g., 1.008, 15.999)
    returns {element.symbol: element.mass for element in mendeleev.get_all_elements()}
    '''
    return {element.symbol: element.mass for element in mlv.get_all_elements()}


def load_mass_dict():
    '''
    Load dictionary of masses for each element using mendeleev (keys are elemental symbols)
    '''
    if not os.path.exists(dump_file):
        main()
    with open(dump_file, 'rb') as f:
        mass = pkl.load(f)
    return mass


def main():
    mass = generate_mass_dict()
    with open(dump_file, 'wb') as f:
        pkl.dump(mass, f)


if __name__ == "__main__":
    main()
