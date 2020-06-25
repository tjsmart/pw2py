#!/usr/bin/env python3

import os
# import mendeleev as mlv
import pickle as pkl


# file where pseudos will be stored
dump_file = os.path.join(os.path.dirname(__file__), 'pseudos.pkl')


def generate_pseudo_dict(pseudo_path='check environment'):
    '''
    Generate dictionary of pseudopotential filenames organized hierarchically by
        pseudos[library]:
            library in ['dojo', 'oncv', 'gbrv']
        pseudos[library][xc_type]:
            xc_type in ['pbe', 'lda', 'pbesol']
                NOTE: for oncv, no pbesol
                NOTE: for dojo, see below
        pseudos[library][xc_type][element]:
            element in ['H', 'C', 'O', ...]

        NOTE: for dojo all pseudopotentials are named identically so simply use pseudos[library][element]
    '''

    # path to pseudopotential files
    if pseudo_path == 'check environment':
        pseudo_path = os.environ['ESPRESSO_PSEUDO']
    if not os.path.isdir(pseudo_path):
        raise FileNotFoundError('Provide the value of pseudo_path with option -p')

    # get list of element symbols from mendeleev
    raise NotImplementedError('Not implemented without mendeleev')
    # elements = [e.symbol for e in mlv.get_all_elements()]

    def checkList(elem, check_list):
        # check obtained list
        if len(check_list) == 1:
            target = check_list[0]
        elif len(check_list) == 0:
            # print("There are no options for {}, skipping it.".format(elem))
            target = None
        else:
            print("There are multiple options for {}".format(elem))
            target = None
            while target not in check_list:
                target = input("Please select a valid option (type it in):\n   {}\n".format(check_list))
        return target

    pseudos = {}
    libraries = ['oncv', 'gbrv', 'dojo']
    xc_types = ['pbe', 'lda', 'pbesol']

    for library in libraries:
        if library == 'dojo':
            pseudos[library] = {}
            for element in elements:
                pseudos[library][element] = '{}.upf'.format(element)
        elif library == 'oncv':
            pseudos[library] = {}
            for xc_type in xc_types:
                if xc_type == 'pbesol':
                    continue
                pseudos[library][xc_type] = {}
                for element in elements:
                    prefix = element + "_{}_{}".format(library.upper(), xc_type.upper())
                    pp_list = [pp for pp in os.listdir(pseudo_path) if pp.startswith(prefix)]
                    if len(pp_list) == 2:
                        # default to picking 1.1 for oncv
                        pp_list = [pp for pp in pp_list if "1.1" in pp]
                    pseudos[library][xc_type][element] = checkList(element, pp_list)
        elif library == 'gbrv':
            pseudos[library] = {}
            for xc_type in xc_types:
                pseudos[library][xc_type] = {}
                for element in elements:
                    if element == 'Hf':
                        # special case use Hf plus4
                        prefix = element.lower() + "_{}_plus4".format(xc_type)
                    else:
                        prefix = element.lower() + "_{}_".format(xc_type)
                    pp_list = [pp for pp in os.listdir(pseudo_path) if pp.startswith(prefix)]
                    pseudos[library][xc_type][element] = checkList(element, pp_list)

    return pseudos


def load_pseudo_dict():
    '''
    Load dictionary of pseudopotential filenames organized hierarchically by
        pseudos[library]:
            library in ['dojo', 'oncv', 'gbrv']
        pseudos[library][xc_type]:
            xc_type in ['pbe', 'lda', 'pbesol']
                NOTE: for oncv, no pbesol
                NOTE: for dojo, see below
        pseudos[library][xc_type][element]:
            element in ['H', 'C', 'O', ...]

        NOTE: for dojo all pseudopotentials are name identically so simply use pseudos[library][element]
    '''
    if not os.path.exists(dump_file):
        main()
    with open(dump_file, 'rb') as f:
        pseudos = pkl.load(f)
    return pseudos


def determine_pseudo_type(pseudo_file: str):
    ''' return library name and xc_type based on pseudo_file name '''
    if 'ONCV' in pseudo_file:
        library = 'oncv'
        if 'PBE' in pseudo_file:
            xc_type = 'pbe'
        elif 'LDA' in pseudo_file:
            xc_type = 'lda'
        else:
            raise ValueError('Unable to determine xc_type for oncv file: {}'.format(pseudo_file))
    elif pseudo_file.endswith('uspp.F.UPF'):
        library = 'gbrv'
        found_xc = False
        for xc_type in ['pbe', 'lda', 'pbesol']:
            if '_{}_'.format(xc_type) in pseudo_file:
                found_xc = True
                break
        if not found_xc:
            raise ValueError('Unable to determine xc_type for gbrv file: {}'.format(pseudo_file))
    elif pseudo_file.endswith('.upf'):
        # pseudo is assumed to be dojo
        library = 'dojo'
        ion = os.path.splitext(pseudo_file)[0]
#         try:
#             # try passing ion to mendeleev to ensure is is indeed formatted this way
#             mlv.element(ion)
#         except:  # noqa: E722
#             # raises sqlalchemy.orm.exc.NoResultFound (rather than import this/require install just catch anything)
#             raise ValueError('Tried to parse pseudo_file as oncv, gbrv, and dojo but failed: {}'.format(pseudo_file))
        xc_type = None

    return library, xc_type


def main():
    pseudos = generate_pseudo_dict()
    with open(dump_file, 'wb') as f:
        pkl.dump(pseudos, f)


if __name__ == "__main__":
    main()
