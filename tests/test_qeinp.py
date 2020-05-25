#!/usr/bin/env python3

import unittest
import numpy as np
import pw2py as pw
from f90nml import Namelist
from copy import deepcopy

import test_atomgeo


class qi_io(unittest.TestCase):
    ''' test from_file method of qeinp class '''

    def run_inp_tests(self, inp, lrelaxed):
        ''' run tests to check input was read correctly '''
        # test ion, par, pos, par_units, pos_units (defined in test_atomgeo)
        test_atomgeo.ag_io.run_inp_tests(self, inp, lrelaxed)
        # extended testing for qeinp
        # check kpt
        self.assertIsNone(np.testing.assert_array_equal(
            inp.kpt, np.array([[8, 8, 8], [0, 0, 0]], dtype=int)
        ))
        # check if_pos
        self.assertIsNone(np.testing.assert_array_equal(
            inp.if_pos, np.ones((inp.nat, 3), dtype=int)
        ))

    def test_from_file_qeinp(self):
        ''' test creating atomgeo from qe input file '''
        inp = pw.qeinp.from_file('input_files/og.in')
        self.run_inp_tests(inp, False)
        # check nml
        self.assertEqual(
            inp.nml, Namelist([
                ('control',
                 Namelist([('calculation', 'vc-relax'), ('wf_collect', True)])
                 ),
                ('system',
                 Namelist([('ibrav', 1), ('a', 6.017), ('nat', 5), ('ntyp', 3), ('ecutwfc', 80)])
                 ),
                ('electrons',
                 Namelist([('conv_thr', 1e-08)])
                 ),
                ('ions',
                 Namelist()
                 ),
                ('cell',
                 Namelist()
                 )
            ])
        )
        # check card
        self.assertEqual(
            inp.card, {
                'ATOMIC_SPECIES': {
                    'Cs': [1.0, 'Cs_ONCV_PBE-1.0.upf'],
                    'Pb': [1.0, 'Pb_ONCV_PBE-1.0.upf'],
                    'Br': [1.0, 'Br_ONCV_PBE-1.0.upf']
                },
                'ATOMIC_POSITIONS': 'crystal',
                'K_POINTS': 'automatic'
            }
        )

    def test_from_file_qeinp_0(self):
        ''' test creating qeinp from qe input file w/ ibrav = 0'''
        inp = pw.qeinp.from_file('input_files/og_0.in')
        self.run_inp_tests(inp, False)
        # check nml
        self.assertEqual(
            inp.nml, Namelist([
                ('control',
                 Namelist([('calculation', 'vc-relax'), ('wf_collect', True)])
                 ),
                ('system',
                 Namelist([('ibrav', 0), ('nat', 5), ('ntyp', 3), ('ecutwfc', 80)])
                 ),
                ('electrons',
                 Namelist([('conv_thr', 1e-08)])
                 ),
                ('ions',
                 Namelist()
                 ),
                ('cell',
                 Namelist()
                 )
            ])
        )
        # check card
        self.assertEqual(
            inp.card, {
                'ATOMIC_SPECIES': {
                    'Cs': [1.0, 'Cs_ONCV_PBE-1.0.upf'],
                    'Pb': [1.0, 'Pb_ONCV_PBE-1.0.upf'],
                    'Br': [1.0, 'Br_ONCV_PBE-1.0.upf']
                },
                'CELL_PARAMETERS': 'angstrom',
                'ATOMIC_POSITIONS': 'crystal',
                'K_POINTS': 'automatic'
            }
        )

    def test_from_file_qeinp_ang(self):
        ''' test creating qeinp from qe input file w/ ibrav = 0'''
        inp = pw.qeinp.from_file('input_files/og_ang.in')
        inp.pos_units = 'crystal'
        self.run_inp_tests(inp, False)
        # check nml
        self.assertEqual(
            inp.nml, Namelist([
                ('control',
                 Namelist([('calculation', 'vc-relax'), ('wf_collect', True)])
                 ),
                ('system',
                 Namelist([('ibrav', 1), ('a', 6.017), ('nat', 5), ('ntyp', 3), ('ecutwfc', 80)])
                 ),
                ('electrons',
                 Namelist([('conv_thr', 1e-08)])
                 ),
                ('ions',
                 Namelist()
                 ),
                ('cell',
                 Namelist()
                 )
            ])
        )
        # check card
        self.assertEqual(
            inp.card, {
                'ATOMIC_SPECIES': {
                    'Cs': [1.0, 'Cs_ONCV_PBE-1.0.upf'],
                    'Pb': [1.0, 'Pb_ONCV_PBE-1.0.upf'],
                    'Br': [1.0, 'Br_ONCV_PBE-1.0.upf']
                },
                'ATOMIC_POSITIONS': 'crystal',
                'K_POINTS': 'automatic'
            }
        )


class qi_methods(unittest.TestCase):
    ''' test various methods of qeinp class '''

    starting_inp = pw.qeinp.from_file('input_files/og.in')

    def test_replace_ion(self):
        ''' test qeinp.replace_ion method '''
        inp = deepcopy(qi_methods.starting_inp)
        inp.replace_ion('Br', 'Cl')
        self.assertIsNone(np.testing.assert_array_equal(
            inp.ion, ['Cs', 'Pb', 'Cl', 'Cl', 'Cl']
        ))
        self.assertEqual(
            inp.card.ATOMIC_SPECIES, {
                'Cs': [1.0, 'Cs_ONCV_PBE-1.0.upf'],
                'Pb': [1.0, 'Pb_ONCV_PBE-1.0.upf'],
                'Cl': [35.45, 'Cl_ONCV_PBE-1.1.upf']
            }
        )
        self.assertEqual(inp.ntyp, 3)
        self.assertEqual(inp.nat, 5)

    def test_add_atom(self):
        ''' test qeinp.add_atom method '''
        inp = deepcopy(qi_methods.starting_inp)
        add_ion = ['H', 'H', 'O']
        add_pos = [[0, 0, 0], [1, 1, 1], [2, 2, 2]]
        inp.add_atom([add_ion, add_pos])
        self.assertIsNone(np.testing.assert_array_equal(
            inp.ion, list(qi_methods.starting_inp.ion) + add_ion
        ))
        self.assertIsNone(np.testing.assert_array_equal(
            inp.pos, list(qi_methods.starting_inp.pos) + add_pos
        ))
        # check ATOMIC SPECIES
        self.assertEqual(
            inp.card.ATOMIC_SPECIES, {
                'Cs': [1.0, 'Cs_ONCV_PBE-1.0.upf'],
                'Pb': [1.0, 'Pb_ONCV_PBE-1.0.upf'],
                'Br': [1.0, 'Br_ONCV_PBE-1.0.upf'],
                'H': [1.008, 'H_ONCV_PBE-1.0.upf'],
                'O': [15.999, 'O_ONCV_PBE-1.0.upf']
            }
        )
        # check ntyp
        self.assertEqual(inp.ntyp, 5)
        self.assertEqual(inp.ntyp, inp.nml.ntyp)
        # check nat
        self.assertEqual(inp.nat, 8)
        self.assertEqual(inp.nat, inp.nml.nat)
        # check if_pos
        self.assertIsNone(np.testing.assert_array_equal(
            inp.if_pos, np.ones((inp.nat, 3), dtype=int)
        ))


if __name__ == '__main__':
    unittest.main()
