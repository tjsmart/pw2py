#!/usr/bin/env python3

import os
import pickle as pkl
import pandas as pd

# Ignore ImportError for users without mendeleev
try:
    import mendeleev
except ImportError:
    pass


db_folder = os.path.join(os.path.dirname(__file__), 'db')

def get_file_path(db_name):
    return os.path.join(db_folder, '{}.pkl'.format(db_name))


def generate_db():
    '''
    creates dictionaries relating element symbols (e.g., 'H' or 'O'), element mass (e.g., 1.008, 15.999), and element numbers

    for mass key, use format .4f

    also creates list for symbols
    '''
    if not os.path.exists(db_folder):
        os.mkdir(db_folder)
    assert mendeleev, "Mendeleev module is reuired to regenerate database"
    # cache mendeleev elements
    all_elements = mendeleev.get_all_elements()
    # make each database
    with open(get_file_path('symbol_to_mass'), 'wb') as f:
        pkl.dump({element.symbol: element.atomic_weight for element in all_elements}, f)
    with open(get_file_path('mass_to_symbol'), 'wb') as f:
        pkl.dump({'{:.4f}'.format(element.atomic_weight): element.symbol for element in all_elements}, f)
    with open(get_file_path('symbol_to_number'), 'wb') as f:
        pkl.dump({element.symbol: element.atomic_number for element in all_elements}, f)
    with open(get_file_path('number_to_symbol'), 'wb') as f:
        pkl.dump({str(element.atomic_number): element.symbol for element in all_elements}, f)
    with open(get_file_path('mass_to_number'), 'wb') as f:
        pkl.dump({'{:.4f}'.format(element.atomic_weight): element.atomic_number for element in all_elements}, f)
    with open(get_file_path('number_to_mass'), 'wb') as f:
        pkl.dump({str(element.atomic_number): element.atomic_weight for element in all_elements}, f)
    with open(get_file_path('symbols'), 'wb') as f:
        pkl.dump([element.symbol for element in all_elements], f)


def load_db(db_name):
    '''
    load appropriate dictionary database
    '''
    db_file_path = get_file_path(db_name)
    if not os.path.exists(db_file_path):
        generate_db()
    with open(db_file_path, 'rb') as f:
        return pkl.load(f)


if __name__ == "__main__":
    generate_db()
