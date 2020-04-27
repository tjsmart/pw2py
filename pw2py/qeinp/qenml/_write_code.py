#!/usr/bin/env python3


def write_code(parameter, namelist):
    print("\
@property\n\
def {0}(self):\n\
    return self.qedict['nml']['{1}']['{0}']\n\
\n\
\n\
@{0}.setter\n\
def {0}(self, {0}):\n\
    self.qedict['nml']['{1}']['{0}'] = {0}\n\n".format(parameter, namelist))


namelists = {
    'control': ['calculation', 'title', 'verbosity', 'restart_mode', 'wf_collect', 'nstep', 'iprint', 'tstress',
                'tprnfor', 'dt', 'outdir', 'wfcdir', 'prefix', 'lkpoint_dir', 'max_seconds', 'etot_conv_thr',
                'forc_conv_thr', 'disk_io', 'pseudo_dir', 'tefield', 'dipfield', 'lelfield', 'nberrycyc',
                'lorbm', 'lberry', 'gdir', 'nppstr', 'lfcpopt', 'monopole'],

    'system': ['ibrav', 'celldm', 'A', 'B', 'C', 'cosAB', 'cosAC', 'cosBC', 'nat', 'ntyp', 'nbnd', 'tot_charge',
               'tot_magnetization', 'starting_magnetization', 'ecutwfc', 'ecutrho', 'ecutfock', 'nr1', 'nr2',
               'nr3', 'nr1s', 'nr2s', 'nr3s', 'nosym', 'nosym_evc', 'noinv', 'no_t_rev', 'force_symmorphic',
               'use_all_frac', 'occupations', 'one_atom_occupations', 'starting_spin_angle', 'degauss',
               'smearing', 'nspin', 'noncolin', 'ecfixed', 'qcutz', 'q2sigma', 'input_dft', 'exx_fraction',
               'screening_parameter', 'exxdiv_treatment', 'x_gamma_extrapolation', 'ecutvcut', 'nqx1', 'nqx2',
               'nqx3', 'lda_plus_u', 'lda_plus_u_kind', 'Hubbard_U', 'Hubbard_J0', 'Hubbard_alpha', 'Hubbard_beta',
               'Hubbard_J', 'starting_ns_eigenvalue', 'U_projection_type', 'edir', 'emaxpos', 'eopreg', 'eamp',
               'angle1', 'angle2', 'constrained_magnetization', 'fixed_magnetization', 'lambda_qe', 'report',
               'lspinorb', 'assume_isolated', 'esm_bc', 'esm_w', 'esm_efield', 'esm_nfit', 'fcp_mu',
               'vdw_corr', 'london', 'london_s6', 'london_c6', 'london_rvdw', 'london_rcut', 'ts_vdw_econv_thr',
               'ts_vdw_isolated', 'xdm', 'xdm_a1', 'xdm_a2', 'space_group', 'uniqueb', 'origin_choice',
               'rhombohedral', 'zmon', 'realxz', 'block', 'block_1', 'block_2', 'block_height'],

    'electrons': ['electron_maxstep', 'scf_must_converge', 'conv_thr', 'adaptive_thr', 'conv_thr_init',
                  'conv_thr_multi', 'mixing_mode', 'mixing_beta', 'mixing_ndim', 'mixing_fixed_ns',
                  'diagonalization', 'ortho_para', 'diago_thr_init', 'diago_cg_maxiter',
                  'diago_david_ndim', 'diago_full_acc', 'efield', 'efield_cart', 'efield_phase',
                  'startingpot', 'startingwfc', 'tqr'],

    'ions': ['ion_dynamics', 'ion_positions', 'pot_extrapolation', 'wfc_extrapolation', 'remove_rigid_rot',
             'ion_temperature', 'tempw', 'tolp', 'delta_t', 'nraise', 'refold_pos', 'upscale', 'bfgs_ndim',
             'trust_radius_max', 'trust_radius_min', 'trust_radius_ini', 'w_1', 'w_2'],

    'cell': ['cell_dynamics', 'press', 'wmass', 'cell_factor', 'press_conv_thr', 'cell_dofree']
}


if __name__ == "__main__":
    ''' script for printing out attributes for qeinp '''

    for namelist in namelists:
        for parameter in namelists[namelist]:
            write_code(parameter, namelist)
