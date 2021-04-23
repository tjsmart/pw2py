# defines all properties from qe input data description


'''
######################################################################################
#################################### &control ########################################
######################################################################################
'''


@property
def calculation(self):
    return self['control']['calculation']


@calculation.setter
def calculation(self, calculation):
    if str(calculation).lower() not in ['scf', 'nscf', 'bands', 'relax', 'md', 'vc-relax', 'vc-md']:
        raise ValueError('Invalid value for calculation')
    self['control']['calculation'] = str(calculation).lower()


@property
def title(self):
    return self['control']['title']


@title.setter
def title(self, title):
    self['control']['title'] = str(title)


@property
def verbosity(self):
    return self['control']['verbosity']


@verbosity.setter
def verbosity(self, verbosity):
    if str(verbosity).lower() not in ['high', 'low']:
        raise ValueError('Invalid value for verbosity')
    self['control']['verbosity'] = str(verbosity).lower()


@property
def restart_mode(self):
    return self['control']['restart_mode']


@restart_mode.setter
def restart_mode(self, restart_mode):
    if str(restart_mode).lower() not in ['from_scratch', 'restart']:
        raise ValueError('Invalid value for restart_mode')
    self['control']['restart_mode'] = str(restart_mode).lower()


@property
def wf_collect(self):
    return self['control']['wf_collect']


@wf_collect.setter
def wf_collect(self, wf_collect):
    self['control']['wf_collect'] = bool(wf_collect)


@property
def nstep(self):
    return self['control']['nstep']


@nstep.setter
def nstep(self, nstep):
    self['control']['nstep'] = int(nstep)


@property
def iprint(self):
    return self['control']['iprint']


@iprint.setter
def iprint(self, iprint):
    self['control']['iprint'] = int(iprint)


@property
def tstress(self):
    return self['control']['tstress']


@tstress.setter
def tstress(self, tstress):
    self['control']['tstress'] = bool(tstress)


@property
def tprnfor(self):
    return self['control']['tprnfor']


@tprnfor.setter
def tprnfor(self, tprnfor):
    self['control']['tprnfor'] = bool(tprnfor)


@property
def dt(self):
    return self['control']['dt']


@dt.setter
def dt(self, dt):
    self['control']['dt'] = float(dt)


@property
def outdir(self):
    return self['control']['outdir']


@outdir.setter
def outdir(self, outdir):
    self['control']['outdir'] = str(outdir)


@property
def wfcdir(self):
    return self['control']['wfcdir']


@wfcdir.setter
def wfcdir(self, wfcdir):
    self['control']['wfcdir'] = str(wfcdir)


@property
def prefix(self):
    return self['control']['prefix']


@prefix.setter
def prefix(self, prefix):
    self['control']['prefix'] = str(prefix)


@property
def lkpoint_dir(self):
    return self['control']['lkpoint_dir']


@lkpoint_dir.setter
def lkpoint_dir(self, lkpoint_dir):
    self['control']['lkpoint_dir'] = bool(lkpoint_dir)


@property
def max_seconds(self):
    return self['control']['max_seconds']


@max_seconds.setter
def max_seconds(self, max_seconds):
    self['control']['max_seconds'] = float(max_seconds)


@property
def etot_conv_thr(self):
    return self['control']['etot_conv_thr']


@etot_conv_thr.setter
def etot_conv_thr(self, etot_conv_thr):
    self['control']['etot_conv_thr'] = float(etot_conv_thr)


@property
def forc_conv_thr(self):
    return self['control']['forc_conv_thr']


@forc_conv_thr.setter
def forc_conv_thr(self, forc_conv_thr):
    self['control']['forc_conv_thr'] = float(forc_conv_thr)


@property
def disk_io(self):
    return self['control']['disk_io']


@disk_io.setter
def disk_io(self, disk_io):
    if str(disk_io).lower() not in ['high', 'medium', 'low', 'none']:
        raise ValueError('Invalid value for disk_io')
    self['control']['disk_io'] = str(disk_io).lower()


@property
def pseudo_dir(self):
    return self['control']['pseudo_dir']


@pseudo_dir.setter
def pseudo_dir(self, pseudo_dir):
    self['control']['pseudo_dir'] = str(pseudo_dir)


@property
def tefield(self):
    return self['control']['tefield']


@tefield.setter
def tefield(self, tefield):
    self['control']['tefield'] = bool(tefield)


@property
def dipfield(self):
    return self['control']['dipfield']


@dipfield.setter
def dipfield(self, dipfield):
    self['control']['dipfield'] = bool(dipfield)


@property
def lelfield(self):
    return self['control']['lelfield']


@lelfield.setter
def lelfield(self, lelfield):
    self['control']['lelfield'] = bool(lelfield)


@property
def nberrycyc(self):
    return self['control']['nberrycyc']


@nberrycyc.setter
def nberrycyc(self, nberrycyc):
    self['control']['nberrycyc'] = int(nberrycyc)


@property
def lorbm(self):
    return self['control']['lorbm']


@lorbm.setter
def lorbm(self, lorbm):
    self['control']['lorbm'] = bool(lorbm)


@property
def lberry(self):
    return self['control']['lberry']


@lberry.setter
def lberry(self, lberry):
    self['control']['lberry'] = bool(lberry)


@property
def gdir(self):
    return self['control']['gdir']


@gdir.setter
def gdir(self, gdir):
    self['control']['gdir'] = int(gdir)


@property
def nppstr(self):
    return self['control']['nppstr']


@nppstr.setter
def nppstr(self, nppstr):
    self['control']['nppstr'] = int(nppstr)


@property
def lfcpopt(self):
    return self['control']['lfcpopt']


@lfcpopt.setter
def lfcpopt(self, lfcpopt):
    self['control']['lfcpopt'] = bool(lfcpopt)


@property
def monopole(self):
    return self['control']['monopole']


@monopole.setter
def monopole(self, monopole):
    self['control']['monopole'] = bool(monopole)


'''
######################################################################################
#################################### &system #########################################
######################################################################################
'''

# TODO
# these require more thinking, perhaps should go in main properties file
# the rest of the properties in this file are dummer than these
@property
def ibrav(self):
    return self['system']['ibrav']


@ibrav.setter
def ibrav(self, ibrav):
    self['system']['ibrav'] = int(ibrav)


@property
def celldm(self):
    return self['system']['celldm']


@celldm.setter
def celldm(self, celldm):
    self['system']['celldm'] = celldm


@property
def a(self):
    return self['system']['a']


@a.setter
def a(self, a):
    self['system']['a'] = a


@property
def b(self):
    return self['system']['b']


@b.setter
def b(self, b):
    self['system']['b'] = b


@property
def c(self):
    return self['system']['c']


@c.setter
def c(self, c):
    self['system']['c'] = c


@property
def cosab(self):
    return self['system']['cosab']


@cosab.setter
def cosab(self, cosab):
    self['system']['cosab'] = cosab


@property
def cosac(self):
    return self['system']['cosac']


@cosac.setter
def cosac(self, cosac):
    self['system']['cosac'] = cosac


@property
def cosbc(self):
    return self['system']['cosbc']


@cosbc.setter
def cosbc(self, cosbc):
    self['system']['cosbc'] = cosbc


@property
def nat(self):
    return self['system']['nat']


@nat.setter
def nat(self, nat):
    self['system']['nat'] = nat


@property
def ntyp(self):
    return self['system']['ntyp']


@ntyp.setter
def ntyp(self, ntyp):
    self['system']['ntyp'] = ntyp


@property
def nbnd(self):
    return self['system']['nbnd']


@nbnd.setter
def nbnd(self, nbnd):
    self['system']['nbnd'] = int(nbnd)


@property
def tot_charge(self):
    return self['system']['tot_charge']


@tot_charge.setter
def tot_charge(self, tot_charge):
    self['system']['tot_charge'] = float(tot_charge)


@property
def tot_magnetization(self):
    return self['system']['tot_magnetization']


@tot_magnetization.setter
def tot_magnetization(self, tot_magnetization):
    self['system']['tot_magnetization'] = float(tot_magnetization)


@property
def starting_magnetization(self):
    return self['system']['starting_magnetization']


@starting_magnetization.setter
def starting_magnetization(self, starting_magnetization):
    self['system']['starting_magnetization'] = [float(sm) for sm in starting_magnetization]


@property
def ecutwfc(self):
    return self['system']['ecutwfc']


@ecutwfc.setter
def ecutwfc(self, ecutwfc):
    self['system']['ecutwfc'] = float(ecutwfc)


@property
def ecutrho(self):
    return self['system']['ecutrho']


@ecutrho.setter
def ecutrho(self, ecutrho):
    self['system']['ecutrho'] = float(ecutrho)


@property
def ecutfock(self):
    return self['system']['ecutfock']


@ecutfock.setter
def ecutfock(self, ecutfock):
    self['system']['ecutfock'] = float(ecutfock)


@property
def nr1(self):
    return self['system']['nr1']


@nr1.setter
def nr1(self, nr1):
    self['system']['nr1'] = int(nr1)


@property
def nr2(self):
    return self['system']['nr2']


@nr2.setter
def nr2(self, nr2):
    self['system']['nr2'] = int(nr2)


@property
def nr3(self):
    return self['system']['nr3']


@nr3.setter
def nr3(self, nr3):
    self['system']['nr3'] = int(nr3)


@property
def nr1s(self):
    return self['system']['nr1s']


@nr1s.setter
def nr1s(self, nr1s):
    self['system']['nr1s'] = int(nr1s)


@property
def nr2s(self):
    return self['system']['nr2s']


@nr2s.setter
def nr2s(self, nr2s):
    self['system']['nr2s'] = int(nr2s)


@property
def nr3s(self):
    return self['system']['nr3s']


@nr3s.setter
def nr3s(self, nr3s):
    self['system']['nr3s'] = int(nr3s)


@property
def nosym(self):
    return self['system']['nosym']


@nosym.setter
def nosym(self, nosym):
    self['system']['nosym'] = bool(nosym)


@property
def nosym_evc(self):
    return self['system']['nosym_evc']


@nosym_evc.setter
def nosym_evc(self, nosym_evc):
    self['system']['nosym_evc'] = bool(nosym_evc)


@property
def noinv(self):
    return self['system']['noinv']


@noinv.setter
def noinv(self, noinv):
    self['system']['noinv'] = bool(noinv)


@property
def no_t_rev(self):
    return self['system']['no_t_rev']


@no_t_rev.setter
def no_t_rev(self, no_t_rev):
    self['system']['no_t_rev'] = bool(no_t_rev)


@property
def force_symmorphic(self):
    return self['system']['force_symmorphic']


@force_symmorphic.setter
def force_symmorphic(self, force_symmorphic):
    self['system']['force_symmorphic'] = bool(force_symmorphic)


@property
def use_all_frac(self):
    return self['system']['use_all_frac']


@use_all_frac.setter
def use_all_frac(self, use_all_frac):
    self['system']['use_all_frac'] = bool(use_all_frac)


@property
def occupations(self):
    return self['system']['occupations']


@occupations.setter
def occupations(self, occupations):
    if str(occupations).lower() not in \
            ['smearing', 'tetrahedra', 'tetrahedra_lin', 'tetrahedra_opt', 'fixed', 'from_input']:
        raise ValueError('Invalid value for occupations')
    self['system']['occupations'] = str(occupations).lower()


@property
def one_atom_occupations(self):
    return self['system']['one_atom_occupations']


@one_atom_occupations.setter
def one_atom_occupations(self, one_atom_occupations):
    self['system']['one_atom_occupations'] = bool(one_atom_occupations)


@property
def starting_spin_angle(self):
    return self['system']['starting_spin_angle']


@starting_spin_angle.setter
def starting_spin_angle(self, starting_spin_angle):
    self['system']['starting_spin_angle'] = bool(starting_spin_angle)


@property
def degauss(self):
    return self['system']['degauss']


@degauss.setter
def degauss(self, degauss):
    self['system']['degauss'] = float(degauss)


@property
def smearing(self):
    return self['system']['smearing']


@smearing.setter
def smearing(self, smearing):
    if str(smearing).lower() not in \
            ['gaussian', 'gauss', 'methfessel-paxton', 'm-p', 'mp',
                'marzari-vanderbilt', 'cold', 'm-v', 'mv', 'fermi-dirac', 'f-d', 'fd']:
        raise ValueError('Invalid value for smearing')
    self['system']['smearing'] = str(smearing).lower()


@property
def nspin(self):
    return self['system']['nspin']


@nspin.setter
def nspin(self, nspin):
    if int(nspin) not in [1, 2, 4]:
        raise ValueError('Invalid value for nspin')
    self['system']['nspin'] = int(nspin)


@property
def noncolin(self):
    return self['system']['noncolin']


@noncolin.setter
def noncolin(self, noncolin):
    self['system']['noncolin'] = bool(noncolin)


@property
def ecfixed(self):
    return self['system']['ecfixed']


@ecfixed.setter
def ecfixed(self, ecfixed):
    self['system']['ecfixed'] = float(ecfixed)


@property
def qcutz(self):
    return self['system']['qcutz']


@qcutz.setter
def qcutz(self, qcutz):
    self['system']['qcutz'] = float(qcutz)


@property
def q2sigma(self):
    return self['system']['q2sigma']


@q2sigma.setter
def q2sigma(self, q2sigma):
    self['system']['q2sigma'] = float(q2sigma)


@property
def input_dft(self):
    return self['system']['input_dft']


@input_dft.setter
def input_dft(self, input_dft):
    self['system']['input_dft'] = str(input_dft)


@property
def exx_fraction(self):
    return self['system']['exx_fraction']


@exx_fraction.setter
def exx_fraction(self, exx_fraction):
    self['system']['exx_fraction'] = float(exx_fraction)


@property
def screening_parameter(self):
    return self['system']['screening_parameter']


@screening_parameter.setter
def screening_parameter(self, screening_parameter):
    self['system']['screening_parameter'] = float(screening_parameter)


@property
def exxdiv_treatment(self):
    return self['system']['exxdiv_treatment']


@exxdiv_treatment.setter
def exxdiv_treatment(self, exxdiv_treatment):
    if str(exxdiv_treatment).lower() not in ['gygi-baldereschi', 'vcut_spherical', 'vcut_ws', 'none']:
        raise ValueError('Invalid value for exxdiv_treatment')
    self['system']['exxdiv_treatment'] = str(exxdiv_treatment).lower()


@property
def x_gamma_extrapolation(self):
    return self['system']['x_gamma_extrapolation']


@x_gamma_extrapolation.setter
def x_gamma_extrapolation(self, x_gamma_extrapolation):
    self['system']['x_gamma_extrapolation'] = bool(x_gamma_extrapolation)


@property
def ecutvcut(self):
    return self['system']['ecutvcut']


@ecutvcut.setter
def ecutvcut(self, ecutvcut):
    self['system']['ecutvcut'] = float(ecutvcut)


@property
def nqx1(self):
    return self['system']['nqx1']


@nqx1.setter
def nqx1(self, nqx1):
    self['system']['nqx1'] = int(nqx1)


@property
def nqx2(self):
    return self['system']['nqx2']


@nqx2.setter
def nqx2(self, nqx2):
    self['system']['nqx2'] = int(nqx2)


@property
def nqx3(self):
    return self['system']['nqx3']


@nqx3.setter
def nqx3(self, nqx3):
    self['system']['nqx3'] = int(nqx3)


@property
def lda_plus_u(self):
    return self['system']['lda_plus_u']


@lda_plus_u.setter
def lda_plus_u(self, lda_plus_u):
    self['system']['lda_plus_u'] = bool(lda_plus_u)


@property
def lda_plus_u_kind(self):
    return self['system']['lda_plus_u_kind']


@lda_plus_u_kind.setter
def lda_plus_u_kind(self, lda_plus_u_kind):
    self['system']['lda_plus_u_kind'] = int(lda_plus_u_kind)


@property
def hubbard_u(self):
    return self['system']['hubbard_u']


@hubbard_u.setter
def hubbard_u(self, hubbard_u):
    self['system']['hubbard_u'] = [float(h) for h in list(hubbard_u)]


@property
def hubbard_j0(self):
    return self['system']['hubbard_j0']


@hubbard_j0.setter
def hubbard_j0(self, hubbard_j0):
    self['system']['hubbard_j0'] = [float(h) for h in list(hubbard_j0)]


@property
def hubbard_alpha(self):
    return self['system']['hubbard_alpha']


@hubbard_alpha.setter
def hubbard_alpha(self, hubbard_alpha):
    self['system']['hubbard_alpha'] = [float(h) for h in list(hubbard_alpha)]


@property
def hubbard_beta(self):
    return self['system']['hubbard_beta']


@hubbard_beta.setter
def hubbard_beta(self, hubbard_beta):
    self['system']['hubbard_beta'] = [float(h) for h in list(hubbard_beta)]


@property
def hubbard_j(self):
    return self['system']['hubbard_j']


@hubbard_j.setter
def hubbard_j(self, hubbard_j):
    # TODO type not understood leaving this open
    self['system']['hubbard_j'] = hubbard_j


@property
def starting_ns_eigenvalue(self):
    return self['system']['starting_ns_eigenvalue']


@starting_ns_eigenvalue.setter
def starting_ns_eigenvalue(self, starting_ns_eigenvalue):
    # TODO type not understand leaving this open
    self['system']['starting_ns_eigenvalue'] = starting_ns_eigenvalue


@property
def u_projection_type(self):
    return self['system']['u_projection_type']


@u_projection_type.setter
def u_projection_type(self, u_projection_type):
    if str(u_projection_type).lower() not in \
            ['atomic', 'ortho-atomic', 'norm-atomic', 'file', 'pseudo']:
        raise ValueError('Invalid value for u_projection_type')
    self['system']['u_projection_type'] = str(u_projection_type).lower()


@property
def edir(self):
    return self['system']['edir']


@edir.setter
def edir(self, edir):
    self['system']['edir'] = int(edir)


@property
def emaxpos(self):
    return self['system']['emaxpos']


@emaxpos.setter
def emaxpos(self, emaxpos):
    self['system']['emaxpos'] = float(emaxpos)


@property
def eopreg(self):
    return self['system']['eopreg']


@eopreg.setter
def eopreg(self, eopreg):
    self['system']['eopreg'] = float(eopreg)


@property
def eamp(self):
    return self['system']['eamp']


@eamp.setter
def eamp(self, eamp):
    self['system']['eamp'] = float(eamp)


@property
def angle1(self):
    return self['system']['angle1']


@angle1.setter
def angle1(self, angle1):
    self['system']['angle1'] = [float(a) for a in list(angle1)]


@property
def angle2(self):
    return self['system']['angle2']


@angle2.setter
def angle2(self, angle2):
    self['system']['angle2'] = [float(a) for a in list(angle2)]


@property
def constrained_magnetization(self):
    return self['system']['constrained_magnetization']


@constrained_magnetization.setter
def constrained_magnetization(self, constrained_magnetization):
    if str(constrained_magnetization).lower() not in \
            ['none', 'total', 'atomic', 'total direction', 'atomic direction']:
        raise ValueError('Invalid value for constrained_magnetization')
    self['system']['constrained_magnetization'] = str(constrained_magnetization).lower()


@property
def fixed_magnetization(self):
    return self['system']['fixed_magnetization']


@fixed_magnetization.setter
def fixed_magnetization(self, fixed_magnetization):
    self['system']['fixed_magnetization'] = \
        [float(fm) for fm in list(fixed_magnetization)[0:3]]


# changed lambda to lambda_qe
@property
def lambda_qe(self):
    return self['system']['lambda']


@lambda_qe.setter
def lambda_qe(self, lambda_qe):
    self['system']['lambda'] = float(lambda_qe)


@property
def report(self):
    return self['system']['report']


@report.setter
def report(self, report):
    self['system']['report'] = int(report)


@property
def lspinorb(self):
    return self['system']['lspinorb']


@lspinorb.setter
def lspinorb(self, lspinorb):
    self['system']['lspinorb'] = bool(lspinorb)


@property
def assume_isolated(self):
    return self['system']['assume_isolated']


@assume_isolated.setter
def assume_isolated(self, assume_isolated):
    if str(assume_isolated).lower() not in \
            ['none', 'makov-payne', 'm-p', 'mp', 'martyna-tuckerman', 'm-t', 'mt', 'esm']:
        raise ValueError('Invalid value for assume_isolated')
    self['system']['assume_isolated'] = str(assume_isolated).lower()


@property
def esm_bc(self):
    return self['system']['esm_bc']


@esm_bc.setter
def esm_bc(self, esm_bc):
    if str(esm_bc).lower() not in \
            ['pbc', 'bc1', 'bc2', 'bc3']:
        raise ValueError('Invalid value for esm_bc')
    self['system']['esm_bc'] = str(esm_bc).lower()


@property
def esm_w(self):
    return self['system']['esm_w']


@esm_w.setter
def esm_w(self, esm_w):
    self['system']['esm_w'] = float(esm_w)


@property
def esm_efield(self):
    return self['system']['esm_efield']


@esm_efield.setter
def esm_efield(self, esm_efield):
    self['system']['esm_efield'] = float(esm_efield)


@property
def esm_nfit(self):
    return self['system']['esm_nfit']


@esm_nfit.setter
def esm_nfit(self, esm_nfit):
    self['system']['esm_nfit'] = int(esm_nfit)


@property
def fcp_mu(self):
    return self['system']['fcp_mu']


@fcp_mu.setter
def fcp_mu(self, fcp_mu):
    self['system']['fcp_mu'] = float(fcp_mu)


@property
def vdw_corr(self):
    return self['system']['vdw_corr']


@vdw_corr.setter
def vdw_corr(self, vdw_corr):
    if str(vdw_corr).lower() not in \
            ['grimme-d2', 'dft-d', 'ts', 'ts-vdw', 'tkatchenko-scheffler', 'xdm']:
        raise ValueError('Invalid value for vdw_corr')
    self['system']['vdw_corr'] = str(vdw_corr).lower()


# Obsolete
# @property
# def london(self):
#     return self['system']['london']


# @london.setter
# def london(self, london):
#     self['system']['london'] = london


@property
def london_s6(self):
    return self['system']['london_s6']


@london_s6.setter
def london_s6(self, london_s6):
    self['system']['london_s6'] = float(london_s6)


@property
def london_c6(self):
    return self['system']['london_c6']


@london_c6.setter
def london_c6(self, london_c6):
    self['system']['london_c6'] = float(london_c6)


@property
def london_rvdw(self):
    return self['system']['london_rvdw']


@london_rvdw.setter
def london_rvdw(self, london_rvdw):
    self['system']['london_rvdw'] = float(london_rvdw)


@property
def london_rcut(self):
    return self['system']['london_rcut']


@london_rcut.setter
def london_rcut(self, london_rcut):
    self['system']['london_rcut'] = float(london_rcut)


@property
def ts_vdw_econv_thr(self):
    return self['system']['ts_vdw_econv_thr']


@ts_vdw_econv_thr.setter
def ts_vdw_econv_thr(self, ts_vdw_econv_thr):
    self['system']['ts_vdw_econv_thr'] = float(ts_vdw_econv_thr)


@property
def ts_vdw_isolated(self):
    return self['system']['ts_vdw_isolated']


@ts_vdw_isolated.setter
def ts_vdw_isolated(self, ts_vdw_isolated):
    self['system']['ts_vdw_isolated'] = bool(ts_vdw_isolated)


# Obsolete
# @property
# def xdm(self):
#     return self['system']['xdm']


# @xdm.setter
# def xdm(self, xdm):
#     self['system']['xdm'] = xdm


@property
def xdm_a1(self):
    return self['system']['xdm_a1']


@xdm_a1.setter
def xdm_a1(self, xdm_a1):
    self['system']['xdm_a1'] = float(xdm_a1)


@property
def xdm_a2(self):
    return self['system']['xdm_a2']


@xdm_a2.setter
def xdm_a2(self, xdm_a2):
    self['system']['xdm_a2'] = float(xdm_a2)


@property
def space_group(self):
    return self['system']['space_group']


@space_group.setter
def space_group(self, space_group):
    self['system']['space_group'] = int(space_group)


@property
def uniqueb(self):
    return self['system']['uniqueb']


@uniqueb.setter
def uniqueb(self, uniqueb):
    self['system']['uniqueb'] = bool(uniqueb)


@property
def origin_choice(self):
    return self['system']['origin_choice']


@origin_choice.setter
def origin_choice(self, origin_choice):
    self['system']['origin_choice'] = int(origin_choice)


@property
def rhombohedral(self):
    return self['system']['rhombohedral']


@rhombohedral.setter
def rhombohedral(self, rhombohedral):
    self['system']['rhombohedral'] = bool(rhombohedral)


@property
def zmon(self):
    return self['system']['zmon']


@zmon.setter
def zmon(self, zmon):
    self['system']['zmon'] = float(zmon)


@property
def realxz(self):
    return self['system']['realxz']


@realxz.setter
def realxz(self, realxz):
    self['system']['realxz'] = bool(realxz)


@property
def block(self):
    return self['system']['block']


@block.setter
def block(self, block):
    self['system']['block'] = bool(block)


@property
def block_1(self):
    return self['system']['block_1']


@block_1.setter
def block_1(self, block_1):
    self['system']['block_1'] = float(block_1)


@property
def block_2(self):
    return self['system']['block_2']


@block_2.setter
def block_2(self, block_2):
    self['system']['block_2'] = float(block_2)


@property
def block_height(self):
    return self['system']['block_height']


@block_height.setter
def block_height(self, block_height):
    self['system']['block_height'] = float(block_height)


'''
######################################################################################
################################## &electrons ########################################
######################################################################################
'''


@property
def electron_maxstep(self):
    return self['electrons']['electron_maxstep']


@electron_maxstep.setter
def electron_maxstep(self, electron_maxstep):
    self['electrons']['electron_maxstep'] = int(electron_maxstep)


@property
def scf_must_converge(self):
    return self['electrons']['scf_must_converge']


@scf_must_converge.setter
def scf_must_converge(self, scf_must_converge):
    self['electrons']['scf_must_converge'] = bool(scf_must_converge)


@property
def conv_thr(self):
    return self['electrons']['conv_thr']


@conv_thr.setter
def conv_thr(self, conv_thr):
    self['electrons']['conv_thr'] = float(conv_thr)


@property
def adaptive_thr(self):
    return self['electrons']['adaptive_thr']


@adaptive_thr.setter
def adaptive_thr(self, adaptive_thr):
    self['electrons']['adaptive_thr'] = bool(adaptive_thr)


@property
def conv_thr_init(self):
    return self['electrons']['conv_thr_init']


@conv_thr_init.setter
def conv_thr_init(self, conv_thr_init):
    self['electrons']['conv_thr_init'] = float(conv_thr_init)


@property
def conv_thr_multi(self):
    return self['electrons']['conv_thr_multi']


@conv_thr_multi.setter
def conv_thr_multi(self, conv_thr_multi):
    self['electrons']['conv_thr_multi'] = float(conv_thr_multi)


@property
def mixing_mode(self):
    return self['electrons']['mixing_mode']


@mixing_mode.setter
def mixing_mode(self, mixing_mode):
    if str(mixing_mode).lower() == 'plain':
        self['system']['mixing_mode'] = str(mixing_mode).lower()
    elif str(mixing_mode) in ['local-TF', 'TF']:
        # don't lower, these are not lower-cased
        self['system']['mixing_mode'] = str(mixing_mode)
    else:
        raise ValueError('Invalid value for mixing_mode')


@property
def mixing_beta(self):
    return self['electrons']['mixing_beta']


@mixing_beta.setter
def mixing_beta(self, mixing_beta):
    self['electrons']['mixing_beta'] = float(mixing_beta)


@property
def mixing_ndim(self):
    return self['electrons']['mixing_ndim']


@mixing_ndim.setter
def mixing_ndim(self, mixing_ndim):
    self['electrons']['mixing_ndim'] = int(mixing_ndim)


@property
def mixing_fixed_ns(self):
    return self['electrons']['mixing_fixed_ns']


@mixing_fixed_ns.setter
def mixing_fixed_ns(self, mixing_fixed_ns):
    self['electrons']['mixing_fixed_ns'] = int(mixing_fixed_ns)


@property
def diagonalization(self):
    return self['electrons']['diagonalization']


@diagonalization.setter
def diagonalization(self, diagonalization):
    if str(diagonalization).lower() not in ['david', 'cg', 'david-serial', 'cg-serial']:
        raise ValueError('Invalid value for diagonalization')
    self['system']['diagonalization'] = str(diagonalization).lower()


# Obsolete
# @property
# def ortho_para(self):
#     return self['electrons']['ortho_para']


# @ortho_para.setter
# def ortho_para(self, ortho_para):
#     self['electrons']['ortho_para'] = int(ortho_para)


@property
def diago_thr_init(self):
    return self['electrons']['diago_thr_init']


@diago_thr_init.setter
def diago_thr_init(self, diago_thr_init):
    self['electrons']['diago_thr_init'] = float(diago_thr_init)


@property
def diago_cg_maxiter(self):
    return self['electrons']['diago_cg_maxiter']


@diago_cg_maxiter.setter
def diago_cg_maxiter(self, diago_cg_maxiter):
    self['electrons']['diago_cg_maxiter'] = int(diago_cg_maxiter)


@property
def diago_david_ndim(self):
    return self['electrons']['diago_david_ndim']


@diago_david_ndim.setter
def diago_david_ndim(self, diago_david_ndim):
    self['electrons']['diago_david_ndim'] = int(diago_david_ndim)


@property
def diago_full_acc(self):
    return self['electrons']['diago_full_acc']


@diago_full_acc.setter
def diago_full_acc(self, diago_full_acc):
    self['electrons']['diago_full_acc'] = bool(diago_full_acc)


@property
def efield(self):
    return self['electrons']['efield']


@efield.setter
def efield(self, efield):
    self['electrons']['efield'] = float(efield)


@property
def efield_cart(self):
    return self['electrons']['efield_cart']


@efield_cart.setter
def efield_cart(self, efield_cart):
    self['electrons']['efield_cart'] = [float(ec) for ec in list(efield_cart)[0:3]]


@property
def efield_phase(self):
    return self['electrons']['efield_phase']


@efield_phase.setter
def efield_phase(self, efield_phase):
    if str(efield_phase).lower() not in ['read', 'write', 'none']:
        raise ValueError('Invalid value for efield_phase')
    self['system']['efield_phase'] = str(efield_phase).lower()


@property
def startingpot(self):
    return self['electrons']['startingpot']


@startingpot.setter
def startingpot(self, startingpot):
    if str(startingpot).lower() not in ['atomic', 'file']:
        raise ValueError('Invalid value for startingpot')
    self['system']['startingpot'] = str(startingpot).lower()


@property
def startingwfc(self):
    return self['electrons']['startingwfc']


@startingwfc.setter
def startingwfc(self, startingwfc):
    if str(startingwfc).lower() not in ['atomic', 'atomic+random', 'random', 'file']:
        raise ValueError('Invalid value for startingwfc')
    self['system']['startingwfc'] = str(startingwfc).lower()


@property
def tqr(self):
    return self['electrons']['tqr']


@tqr.setter
def tqr(self, tqr):
    self['electrons']['tqr'] = bool(tqr)


'''
######################################################################################
##################################### &ions ##########################################
######################################################################################
'''


@property
def ion_dynamics(self):
    return self['ions']['ion_dynamics']


@ion_dynamics.setter
def ion_dynamics(self, ion_dynamics):
    if str(ion_dynamics).lower() not in ['bfgs', 'damp', 'verlet', 'langevin',
                                         'langevin-smc', 'beeman']:
        raise ValueError('Invalid value for ion_dynamics')
    self['system']['ion_dynamics'] = str(ion_dynamics).lower()


@property
def ion_positions(self):
    return self['ions']['ion_positions']


@ion_positions.setter
def ion_positions(self, ion_positions):
    if str(ion_positions).lower() not in ['default', 'from_input']:
        raise ValueError('Invalid value for ion_positions')
    self['system']['ion_positions'] = str(ion_positions).lower()


@property
def pot_extrapolation(self):
    return self['ions']['pot_extrapolation']


@pot_extrapolation.setter
def pot_extrapolation(self, pot_extrapolation):
    if str(pot_extrapolation).lower() not in ['none', 'atomic', 'first_order', 'second_order']:
        raise ValueError('Invalid value for pot_extrapolation')
    self['system']['pot_extrapolation'] = str(pot_extrapolation).lower()


@property
def wfc_extrapolation(self):
    return self['ions']['wfc_extrapolation']


@wfc_extrapolation.setter
def wfc_extrapolation(self, wfc_extrapolation):
    if str(wfc_extrapolation).lower() not in ['none', 'first_order', 'second_order']:
        raise ValueError('Invalid value for wfc_extrapolation')
    self['system']['wfc_extrapolation'] = str(wfc_extrapolation).lower()


@property
def remove_rigid_rot(self):
    return self['ions']['remove_rigid_rot']


@remove_rigid_rot.setter
def remove_rigid_rot(self, remove_rigid_rot):
    self['ions']['remove_rigid_rot'] = bool(remove_rigid_rot)


@property
def ion_temperature(self):
    return self['ions']['ion_temperature']


@ion_temperature.setter
def ion_temperature(self, ion_temperature):
    if str(ion_temperature).lower() in ['rescaling', 'rescaling-v', 'berendsen',
                                        'andersen', 'initial', 'not_controlled']:
        self['system']['ion_temperature'] = str(ion_temperature).lower()
    elif str(ion_temperature) in ['rescaling-T', 'reduce-T']:
        # don't lower, these are not lower-cased
        self['system']['ion_temperature'] = str(ion_temperature)
    else:
        raise ValueError('Invalid value for ion_temperature')


@property
def tempw(self):
    return self['ions']['tempw']


@tempw.setter
def tempw(self, tempw):
    self['ions']['tempw'] = float(tempw)


@property
def tolp(self):
    return self['ions']['tolp']


@tolp.setter
def tolp(self, tolp):
    self['ions']['tolp'] = float(tolp)


@property
def delta_t(self):
    return self['ions']['delta_t']


@delta_t.setter
def delta_t(self, delta_t):
    self['ions']['delta_t'] = float(delta_t)


@property
def nraise(self):
    return self['ions']['nraise']


@nraise.setter
def nraise(self, nraise):
    self['ions']['nraise'] = int(nraise)


@property
def refold_pos(self):
    return self['ions']['refold_pos']


@refold_pos.setter
def refold_pos(self, refold_pos):
    self['ions']['refold_pos'] = bool(refold_pos)


@property
def upscale(self):
    return self['ions']['upscale']


@upscale.setter
def upscale(self, upscale):
    self['ions']['upscale'] = float(upscale)


@property
def bfgs_ndim(self):
    return self['ions']['bfgs_ndim']


@bfgs_ndim.setter
def bfgs_ndim(self, bfgs_ndim):
    self['ions']['bfgs_ndim'] = int(bfgs_ndim)


@property
def trust_radius_max(self):
    return self['ions']['trust_radius_max']


@trust_radius_max.setter
def trust_radius_max(self, trust_radius_max):
    self['ions']['trust_radius_max'] = float(trust_radius_max)


@property
def trust_radius_min(self):
    return self['ions']['trust_radius_min']


@trust_radius_min.setter
def trust_radius_min(self, trust_radius_min):
    self['ions']['trust_radius_min'] = float(trust_radius_min)


@property
def trust_radius_ini(self):
    return self['ions']['trust_radius_ini']


@trust_radius_ini.setter
def trust_radius_ini(self, trust_radius_ini):
    self['ions']['trust_radius_ini'] = float(trust_radius_ini)


@property
def w_1(self):
    return self['ions']['w_1']


@w_1.setter
def w_1(self, w_1):
    self['ions']['w_1'] = float(w_1)


@property
def w_2(self):
    return self['ions']['w_2']


@w_2.setter
def w_2(self, w_2):
    self['ions']['w_2'] = float(w_2)


'''
######################################################################################
##################################### &cell ##########################################
######################################################################################
'''


@property
def cell_dynamics(self):
    return self['cell']['cell_dynamics']


@cell_dynamics.setter
def cell_dynamics(self, cell_dynamics):
    if str(cell_dynamics).lower() not in ['none', 'sd', 'damp-pr', 'damp-w', 'bfgs', 'pr', 'w']:
        raise ValueError('Invalid value for cell_dynamics')
    self['system']['cell_dynamics'] = str(cell_dynamics).lower()


@property
def press(self):
    return self['cell']['press']


@press.setter
def press(self, press):
    self['cell']['press'] = float(press)


@property
def wmass(self):
    return self['cell']['wmass']


@wmass.setter
def wmass(self, wmass):
    self['cell']['wmass'] = float(wmass)


@property
def cell_factor(self):
    return self['cell']['cell_factor']


@cell_factor.setter
def cell_factor(self, cell_factor):
    self['cell']['cell_factor'] = float(cell_factor)


@property
def press_conv_thr(self):
    return self['cell']['press_conv_thr']


@press_conv_thr.setter
def press_conv_thr(self, press_conv_thr):
    self['cell']['press_conv_thr'] = float(press_conv_thr)


@property
def cell_dofree(self):
    return self['cell']['cell_dofree']


@cell_dofree.setter
def cell_dofree(self, cell_dofree):
    if str(cell_dofree).lower() not in ['all', 'x', 'y', 'z',
                                        'xy', 'xz', 'yz', 'xyz', 'shape', 'volume', '2Dxy', '2Dshape']:
        raise ValueError('Invalid value for cell_dofree')
    self['system']['cell_dofree'] = str(cell_dofree).lower()
