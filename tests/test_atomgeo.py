#!/usr/bin/env python3

import unittest
import numpy as np
import pw2py as pw
from copy import deepcopy


expected = {
    'ion': ['Cs', 'Pb', 'Br', 'Br', 'Br'],
    'par_initial': 6.017 * np.eye(3),
    'par_final': 6.00537872 * np.eye(3),
    'pos': [
        [0., 0., 0.],
        [0.5, 0.5, 0.5],
        [0., 0.5, 0.5],
        [0.5, 0., 0.5],
        [0.5, 0.5, 0.]
    ],
    'par_units': 'angstrom',
    'pos_units': 'crystal'
}


class ag_io(unittest.TestCase):
    ''' test from_file method of atomgeo class '''

    # precision for asserting array equivalence
    decimal = 5

    def run_inp_tests(self, geo, lrelaxed):
        ''' run tests to check input was read correctly '''
        # test ion
        self.assertIsNone(np.testing.assert_array_equal(geo.ion, expected['ion']))
        # test par
        if lrelaxed:
            self.assertIsNone(np.testing.assert_array_almost_equal(
                geo.par, expected['par_final'], decimal=ag_io.decimal
            ))
        else:
            self.assertIsNone(np.testing.assert_array_equal(geo.par, expected['par_initial']))
        # test pos
        self.assertIsNone(np.testing.assert_array_almost_equal(
            geo.pos, expected['pos'], decimal=ag_io.decimal
        ))
        # test par_units
        self.assertEqual(geo.par_units, expected['par_units'])
        # test pos_units
        self.assertEqual(geo.pos_units, expected['pos_units'])

    def test_from_file_qeinp(self):
        ''' test creating atomgeo from qe input file '''
        self.run_inp_tests(pw.atomgeo.from_file('input_files/og.in'), False)

    def test_from_file_qeinp_0(self):
        ''' test creating atomgeo from qe input file w/ ibrav = 0'''
        self.run_inp_tests(pw.atomgeo.from_file('input_files/og_0.in'), False)

    def test_from_file_qeinp_ang(self):
        ''' test creating atomgeo from qe input file w/ ibrav = 0'''
        geo = pw.atomgeo.from_file('input_files/og_ang.in')
        geo.pos_units = 'crystal'
        self.run_inp_tests(geo, False)

    def test_from_file_vasp(self):
        ''' test creating atomgeo from vasp file '''
        self.run_inp_tests(pw.atomgeo.from_file('input_files/og.vasp'), True)

    def test_from_file_qeout(self):
        ''' test creating atomgeo from qeout file '''
        self.run_inp_tests(pw.atomgeo.from_file('input_files/og.out'), True)

    def test_from_file_vasp_cart(self):
        ''' test creating atomgeo from vasp file with cartesian then converting to crystal '''
        geo = pw.atomgeo.from_file('input_files/og_cart.vasp')
        geo.pos_units = 'crystal'
        self.run_inp_tests(geo, True)

    def test_from_file_xsf(self):
        ''' test creating atomgeo from xsf file '''
        geo = pw.atomgeo.from_file('input_files/og.xsf')
        geo.pos_units = 'crystal'
        self.run_inp_tests(geo, True)

    def test_from_file_xyz(self):
        ''' test creating atomgeo from xyz file '''
        geo = pw.atomgeo.from_file('input_files/og.xyz', xyz_par=expected['par_final'])
        geo.pos_units = 'crystal'
        self.run_inp_tests(geo, True)


class ag_methods(unittest.TestCase):
    ''' test various methods of atomgeo class '''

    starting_geo = pw.atomgeo.from_file('input_files/og.in')

    def test_replace_ion(self):
        ''' test atomgeo.replace_ion method '''
        geo = deepcopy(ag_methods.starting_geo)
        geo.replace_ion('Br', 'Cl')
        self.assertIsNone(np.testing.assert_array_equal(
            geo.ion, ['Cs', 'Pb', 'Cl', 'Cl', 'Cl']
        ))
        # check nat
        self.assertEqual(geo.nat, 5)

    def test_add_atom(self):
        ''' test atomgeo.add_atom method '''
        geo = deepcopy(ag_methods.starting_geo)
        add_ion = ['H', 'H', 'O']
        add_pos = [[0, 0, 0], [1, 1, 1], [2, 2, 2]]
        geo.add_atom([add_ion, add_pos])
        self.assertIsNone(np.testing.assert_array_equal(
            geo.ion, list(ag_methods.starting_geo.ion) + add_ion
        ))
        self.assertIsNone(np.testing.assert_array_equal(
            geo.pos, list(ag_methods.starting_geo.pos) + add_pos
        ))
        # check nat
        self.assertEqual(geo.nat, 8)

    def test_remove_indices(self):
        ''' test atomgeo.remove_indices method '''
        geo = deepcopy(ag_methods.starting_geo)
        indices = [0, 4]
        not_indices = [i for i in range(geo.nat) if i not in indices]
        geo.remove_indices(indices)
        self.assertEqual(geo.nat, ag_methods.starting_geo.nat - len(indices))
        self.assertIsNone(np.testing.assert_array_equal(
            geo.ion, ag_methods.starting_geo.ion[not_indices]
        ))
        self.assertIsNone(np.testing.assert_array_equal(
            geo.pos, ag_methods.starting_geo.pos[not_indices]
        ))

    def test_sort_ions(self):
        ''' test atomgeo.sort_ions method '''
        geo = ag_methods.starting_geo.sort_ions(inplace=False)
        self.assertIsNone(np.testing.assert_array_equal(
            geo.ion, ['Br', 'Br', 'Br', 'Cs', 'Pb']
        ))
        self.assertIsNone(np.testing.assert_array_equal(
            geo.pos, [
                [0., 0.5, 0.5],
                [0.5, 0., 0.5],
                [0.5, 0.5, 0.],
                [0., 0., 0.],
                [0.5, 0.5, 0.5]
            ]
        ))

    # def test_shift_pos_to_unit(self):
    #     ''' test atomgeo.shift_pos_to_unit method '''
    #     TODO
    #     pass

    def test_build_supercell(self):
        ''' test atomgeo.shift_pos_to_unit method '''
        geo = ag_methods.starting_geo.build_supercell([2, 4, 3], inplace=False)
        self.assertIsNone(np.testing.assert_array_equal(
            geo.ion, ['Cs', 'Pb', 'Br', 'Br', 'Br', 'Cs', 'Pb', 'Br', 'Br', 'Br', 'Cs',
                      'Pb', 'Br', 'Br', 'Br', 'Cs', 'Pb', 'Br', 'Br', 'Br', 'Cs', 'Pb',
                      'Br', 'Br', 'Br', 'Cs', 'Pb', 'Br', 'Br', 'Br', 'Cs', 'Pb', 'Br',
                      'Br', 'Br', 'Cs', 'Pb', 'Br', 'Br', 'Br', 'Cs', 'Pb', 'Br', 'Br',
                      'Br', 'Cs', 'Pb', 'Br', 'Br', 'Br', 'Cs', 'Pb', 'Br', 'Br', 'Br',
                      'Cs', 'Pb', 'Br', 'Br', 'Br', 'Cs', 'Pb', 'Br', 'Br', 'Br', 'Cs',
                      'Pb', 'Br', 'Br', 'Br', 'Cs', 'Pb', 'Br', 'Br', 'Br', 'Cs', 'Pb',
                      'Br', 'Br', 'Br', 'Cs', 'Pb', 'Br', 'Br', 'Br', 'Cs', 'Pb', 'Br',
                      'Br', 'Br', 'Cs', 'Pb', 'Br', 'Br', 'Br', 'Cs', 'Pb', 'Br', 'Br',
                      'Br', 'Cs', 'Pb', 'Br', 'Br', 'Br', 'Cs', 'Pb', 'Br', 'Br', 'Br',
                      'Cs', 'Pb', 'Br', 'Br', 'Br', 'Cs', 'Pb', 'Br', 'Br', 'Br']
        ))
        self.assertIsNone(np.testing.assert_array_almost_equal(
            geo.pos, [
                [0., 0., 0.],
                [0.25, 0.125, 0.16666667],
                [0., 0.125, 0.16666667],
                [0.25, 0., 0.16666667],
                [0.25, 0.125, 0.],
                [0., 0., 0.33333333],
                [0.25, 0.125, 0.5],
                [0., 0.125, 0.5],
                [0.25, 0., 0.5],
                [0.25, 0.125, 0.33333333],
                [0., 0., 0.66666667],
                [0.25, 0.125, 0.83333333],
                [0., 0.125, 0.83333333],
                [0.25, 0., 0.83333333],
                [0.25, 0.125, 0.66666667],
                [0., 0.25, 0.],
                [0.25, 0.375, 0.16666667],
                [0., 0.375, 0.16666667],
                [0.25, 0.25, 0.16666667],
                [0.25, 0.375, 0.],
                [0., 0.25, 0.33333333],
                [0.25, 0.375, 0.5],
                [0., 0.375, 0.5],
                [0.25, 0.25, 0.5],
                [0.25, 0.375, 0.33333333],
                [0., 0.25, 0.66666667],
                [0.25, 0.375, 0.83333333],
                [0., 0.375, 0.83333333],
                [0.25, 0.25, 0.83333333],
                [0.25, 0.375, 0.66666667],
                [0., 0.5, 0.],
                [0.25, 0.625, 0.16666667],
                [0., 0.625, 0.16666667],
                [0.25, 0.5, 0.16666667],
                [0.25, 0.625, 0.],
                [0., 0.5, 0.33333333],
                [0.25, 0.625, 0.5],
                [0., 0.625, 0.5],
                [0.25, 0.5, 0.5],
                [0.25, 0.625, 0.33333333],
                [0., 0.5, 0.66666667],
                [0.25, 0.625, 0.83333333],
                [0., 0.625, 0.83333333],
                [0.25, 0.5, 0.83333333],
                [0.25, 0.625, 0.66666667],
                [0., 0.75, 0.],
                [0.25, 0.875, 0.16666667],
                [0., 0.875, 0.16666667],
                [0.25, 0.75, 0.16666667],
                [0.25, 0.875, 0.],
                [0., 0.75, 0.33333333],
                [0.25, 0.875, 0.5],
                [0., 0.875, 0.5],
                [0.25, 0.75, 0.5],
                [0.25, 0.875, 0.33333333],
                [0., 0.75, 0.66666667],
                [0.25, 0.875, 0.83333333],
                [0., 0.875, 0.83333333],
                [0.25, 0.75, 0.83333333],
                [0.25, 0.875, 0.66666667],
                [0.5, 0., 0.],
                [0.75, 0.125, 0.16666667],
                [0.5, 0.125, 0.16666667],
                [0.75, 0., 0.16666667],
                [0.75, 0.125, 0.],
                [0.5, 0., 0.33333333],
                [0.75, 0.125, 0.5],
                [0.5, 0.125, 0.5],
                [0.75, 0., 0.5],
                [0.75, 0.125, 0.33333333],
                [0.5, 0., 0.66666667],
                [0.75, 0.125, 0.83333333],
                [0.5, 0.125, 0.83333333],
                [0.75, 0., 0.83333333],
                [0.75, 0.125, 0.66666667],
                [0.5, 0.25, 0.],
                [0.75, 0.375, 0.16666667],
                [0.5, 0.375, 0.16666667],
                [0.75, 0.25, 0.16666667],
                [0.75, 0.375, 0.],
                [0.5, 0.25, 0.33333333],
                [0.75, 0.375, 0.5],
                [0.5, 0.375, 0.5],
                [0.75, 0.25, 0.5],
                [0.75, 0.375, 0.33333333],
                [0.5, 0.25, 0.66666667],
                [0.75, 0.375, 0.83333333],
                [0.5, 0.375, 0.83333333],
                [0.75, 0.25, 0.83333333],
                [0.75, 0.375, 0.66666667],
                [0.5, 0.5, 0.],
                [0.75, 0.625, 0.16666667],
                [0.5, 0.625, 0.16666667],
                [0.75, 0.5, 0.16666667],
                [0.75, 0.625, 0.],
                [0.5, 0.5, 0.33333333],
                [0.75, 0.625, 0.5],
                [0.5, 0.625, 0.5],
                [0.75, 0.5, 0.5],
                [0.75, 0.625, 0.33333333],
                [0.5, 0.5, 0.66666667],
                [0.75, 0.625, 0.83333333],
                [0.5, 0.625, 0.83333333],
                [0.75, 0.5, 0.83333333],
                [0.75, 0.625, 0.66666667],
                [0.5, 0.75, 0.],
                [0.75, 0.875, 0.16666667],
                [0.5, 0.875, 0.16666667],
                [0.75, 0.75, 0.16666667],
                [0.75, 0.875, 0.],
                [0.5, 0.75, 0.33333333],
                [0.75, 0.875, 0.5],
                [0.5, 0.875, 0.5],
                [0.75, 0.75, 0.5],
                [0.75, 0.875, 0.33333333],
                [0.5, 0.75, 0.66666667],
                [0.75, 0.875, 0.83333333],
                [0.5, 0.875, 0.83333333],
                [0.75, 0.75, 0.83333333],
                [0.75, 0.875, 0.66666667]
            ]
        ))
        self.assertIsNone(np.testing.assert_array_almost_equal(
            geo.par, [
                [12.034,  0.,  0.],
                [0., 24.068,  0.],
                [0.,  0., 18.051]
            ]
        ))


if __name__ == '__main__':
    unittest.main()
