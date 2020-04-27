# defines all properties from qe input data description


'''
######################################################################################
#################################### &control ########################################
######################################################################################
'''


@property
def calculation(self):
    return self._nml['control']['calculation']


@calculation.setter
def calculation(self, calculation):
    if str(calculation).lower() not in ['scf', 'nscf', 'bands', 'relax', 'md', 'vc-relax', 'vc-md']:
        raise ValueError('Invalid value for calculation')
    self._nml['control']['calculation'] = str(calculation).lower()


@property
def title(self):
    return self._nml['control']['title']


@title.setter
def title(self, title):
    self._nml['control']['title'] = str(title)


@property
def verbosity(self):
    return self._nml['control']['verbosity']


@verbosity.setter
def verbosity(self, verbosity):
    if str(verbosity).lower() not in ['high', 'low']:
        raise ValueError('Invalid value for verbosity')
    self._nml['control']['verbosity'] = str(verbosity).lower()


@property
def restart_mode(self):
    return self._nml['control']['restart_mode']


@restart_mode.setter
def restart_mode(self, restart_mode):
    if str(restart_mode).lower() not in ['from_scratch', 'restart']:
        raise ValueError('Invalid value for restart_mode')
    self._nml['control']['restart_mode'] = str(restart_mode).lower()


@property
def wf_collect(self):
    return self._nml['control']['wf_collect']


@wf_collect.setter
def wf_collect(self, wf_collect):
    self._nml['control']['wf_collect'] = bool(wf_collect)


@property
def nstep(self):
    return self._nml['control']['nstep']


@nstep.setter
def nstep(self, nstep):
    self._nml['control']['nstep'] = int(nstep)


@property
def iprint(self):
    return self._nml['control']['iprint']


@iprint.setter
def iprint(self, iprint):
    self._nml['control']['iprint'] = int(iprint)


@property
def tstress(self):
    return self._nml['control']['tstress']


@tstress.setter
def tstress(self, tstress):
    self._nml['control']['tstress'] = bool(tstress)


@property
def tprnfor(self):
    return self._nml['control']['tprnfor']


@tprnfor.setter
def tprnfor(self, tprnfor):
    self._nml['control']['tprnfor'] = bool(tprnfor)


@property
def dt(self):
    return self._nml['control']['dt']


@dt.setter
def dt(self, dt):
    self._nml['control']['dt'] = float(dt)


@property
def outdir(self):
    return self._nml['control']['outdir']


@outdir.setter
def outdir(self, outdir):
    self._nml['control']['outdir'] = str(outdir)


@property
def wfcdir(self):
    return self._nml['control']['wfcdir']


@wfcdir.setter
def wfcdir(self, wfcdir):
    self._nml['control']['wfcdir'] = str(wfcdir)


@property
def prefix(self):
    return self._nml['control']['prefix']


@prefix.setter
def prefix(self, prefix):
    self._nml['control']['prefix'] = str(prefix)


@property
def lkpoint_dir(self):
    return self._nml['control']['lkpoint_dir']


@lkpoint_dir.setter
def lkpoint_dir(self, lkpoint_dir):
    self._nml['control']['lkpoint_dir'] = bool(lkpoint_dir)


@property
def max_seconds(self):
    return self._nml['control']['max_seconds']


@max_seconds.setter
def max_seconds(self, max_seconds):
    self._nml['control']['max_seconds'] = float(max_seconds)


@property
def etot_conv_thr(self):
    return self._nml['control']['etot_conv_thr']


@etot_conv_thr.setter
def etot_conv_thr(self, etot_conv_thr):
    self._nml['control']['etot_conv_thr'] = float(etot_conv_thr)


@property
def forc_conv_thr(self):
    return self._nml['control']['forc_conv_thr']


@forc_conv_thr.setter
def forc_conv_thr(self, forc_conv_thr):
    self._nml['control']['forc_conv_thr'] = float(forc_conv_thr)


@property
def disk_io(self):
    return self._nml['control']['disk_io']


@disk_io.setter
def disk_io(self, disk_io):
    if str(disk_io).lower() not in ['high', 'medium', 'low', 'none']:
        raise ValueError('Invalid value for disk_io')
    self._nml['control']['disk_io'] = str(disk_io).lower()


@property
def pseudo_dir(self):
    return self._nml['control']['pseudo_dir']


@pseudo_dir.setter
def pseudo_dir(self, pseudo_dir):
    self._nml['control']['pseudo_dir'] = str(pseudo_dir)


@property
def tefield(self):
    return self._nml['control']['tefield']


@tefield.setter
def tefield(self, tefield):
    self._nml['control']['tefield'] = bool(tefield)


@property
def dipfield(self):
    return self._nml['control']['dipfield']


@dipfield.setter
def dipfield(self, dipfield):
    self._nml['control']['dipfield'] = bool(dipfield)


@property
def lelfield(self):
    return self._nml['control']['lelfield']


@lelfield.setter
def lelfield(self, lelfield):
    self._nml['control']['lelfield'] = bool(lelfield)


@property
def nberrycyc(self):
    return self._nml['control']['nberrycyc']


@nberrycyc.setter
def nberrycyc(self, nberrycyc):
    self._nml['control']['nberrycyc'] = int(nberrycyc)


@property
def lorbm(self):
    return self._nml['control']['lorbm']


@lorbm.setter
def lorbm(self, lorbm):
    self._nml['control']['lorbm'] = bool(lorbm)


@property
def lberry(self):
    return self._nml['control']['lberry']


@lberry.setter
def lberry(self, lberry):
    self._nml['control']['lberry'] = bool(lberry)


@property
def gdir(self):
    return self._nml['control']['gdir']


@gdir.setter
def gdir(self, gdir):
    self._nml['control']['gdir'] = int(gdir)


@property
def nppstr(self):
    return self._nml['control']['nppstr']


@nppstr.setter
def nppstr(self, nppstr):
    self._nml['control']['nppstr'] = int(nppstr)


@property
def lfcpopt(self):
    return self._nml['control']['lfcpopt']


@lfcpopt.setter
def lfcpopt(self, lfcpopt):
    self._nml['control']['lfcpopt'] = bool(lfcpopt)


@property
def monopole(self):
    return self._nml['control']['monopole']


@monopole.setter
def monopole(self, monopole):
    self._nml['control']['monopole'] = bool(monopole)


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
    return self._nml['system']['ibrav']


@ibrav.setter
def ibrav(self, ibrav):
    self._nml['system']['ibrav'] = int(ibrav)


@property
def celldm(self):
    return self._nml['system']['celldm']


@celldm.setter
def celldm(self, celldm):
    self._nml['system']['celldm'] = celldm


@property
def A(self):
    return self._nml['system']['A']


@A.setter
def A(self, A):
    self._nml['system']['A'] = A


@property
def B(self):
    return self._nml['system']['B']


@B.setter
def B(self, B):
    self._nml['system']['B'] = B


@property
def C(self):
    return self._nml['system']['C']


@C.setter
def C(self, C):
    self._nml['system']['C'] = C


@property
def cosAB(self):
    return self._nml['system']['cosAB']


@cosAB.setter
def cosAB(self, cosAB):
    self._nml['system']['cosAB'] = cosAB


@property
def cosAC(self):
    return self._nml['system']['cosAC']


@cosAC.setter
def cosAC(self, cosAC):
    self._nml['system']['cosAC'] = cosAC


@property
def cosBC(self):
    return self._nml['system']['cosBC']


@cosBC.setter
def cosBC(self, cosBC):
    self._nml['system']['cosBC'] = cosBC


@property
def nat(self):
    return self._nml['system']['nat']


@nat.setter
def nat(self, nat):
    self._nml['system']['nat'] = nat


@property
def ntyp(self):
    return self._nml['system']['ntyp']


@ntyp.setter
def ntyp(self, ntyp):
    self._nml['system']['ntyp'] = ntyp


@property
def nbnd(self):
    return self._nml['system']['nbnd']


@nbnd.setter
def nbnd(self, nbnd):
    self._nml['system']['nbnd'] = int(nbnd)


@property
def tot_charge(self):
    return self._nml['system']['tot_charge']


@tot_charge.setter
def tot_charge(self, tot_charge):
    self._nml['system']['tot_charge'] = float(tot_charge)


@property
def tot_magnetization(self):
    return self._nml['system']['tot_magnetization']


@tot_magnetization.setter
def tot_magnetization(self, tot_magnetization):
    self._nml['system']['tot_magnetization'] = float(tot_magnetization)


@property
def starting_magnetization(self):
    return self._nml['system']['starting_magnetization']


@starting_magnetization.setter
def starting_magnetization(self, starting_magnetization):
    self._nml['system']['starting_magnetization'] = [float(sm) for sm in starting_magnetization]


@property
def ecutwfc(self):
    return self._nml['system']['ecutwfc']


@ecutwfc.setter
def ecutwfc(self, ecutwfc):
    self._nml['system']['ecutwfc'] = float(ecutwfc)


@property
def ecutrho(self):
    return self._nml['system']['ecutrho']


@ecutrho.setter
def ecutrho(self, ecutrho):
    self._nml['system']['ecutrho'] = float(ecutrho)


@property
def ecutfock(self):
    return self._nml['system']['ecutfock']


@ecutfock.setter
def ecutfock(self, ecutfock):
    self._nml['system']['ecutfock'] = float(ecutfock)


@property
def nr1(self):
    return self._nml['system']['nr1']


@nr1.setter
def nr1(self, nr1):
    self._nml['system']['nr1'] = int(nr1)


@property
def nr2(self):
    return self._nml['system']['nr2']


@nr2.setter
def nr2(self, nr2):
    self._nml['system']['nr2'] = int(nr2)


@property
def nr3(self):
    return self._nml['system']['nr3']


@nr3.setter
def nr3(self, nr3):
    self._nml['system']['nr3'] = int(nr3)


@property
def nr1s(self):
    return self._nml['system']['nr1s']


@nr1s.setter
def nr1s(self, nr1s):
    self._nml['system']['nr1s'] = int(nr1s)


@property
def nr2s(self):
    return self._nml['system']['nr2s']


@nr2s.setter
def nr2s(self, nr2s):
    self._nml['system']['nr2s'] = int(nr2s)


@property
def nr3s(self):
    return self._nml['system']['nr3s']


@nr3s.setter
def nr3s(self, nr3s):
    self._nml['system']['nr3s'] = int(nr3s)


@property
def nosym(self):
    return self._nml['system']['nosym']


@nosym.setter
def nosym(self, nosym):
    self._nml['system']['nosym'] = bool(nosym)


@property
def nosym_evc(self):
    return self._nml['system']['nosym_evc']


@nosym_evc.setter
def nosym_evc(self, nosym_evc):
    self._nml['system']['nosym_evc'] = bool(nosym_evc)


@property
def noinv(self):
    return self._nml['system']['noinv']


@noinv.setter
def noinv(self, noinv):
    self._nml['system']['noinv'] = bool(noinv)


@property
def no_t_rev(self):
    return self._nml['system']['no_t_rev']


@no_t_rev.setter
def no_t_rev(self, no_t_rev):
    self._nml['system']['no_t_rev'] = bool(no_t_rev)


@property
def force_symmorphic(self):
    return self._nml['system']['force_symmorphic']


@force_symmorphic.setter
def force_symmorphic(self, force_symmorphic):
    self._nml['system']['force_symmorphic'] = bool(force_symmorphic)


@property
def use_all_frac(self):
    return self._nml['system']['use_all_frac']


@use_all_frac.setter
def use_all_frac(self, use_all_frac):
    self._nml['system']['use_all_frac'] = bool(use_all_frac)


@property
def occupations(self):
    return self._nml['system']['occupations']


@occupations.setter
def occupations(self, occupations):
    if str(occupations).lower() not in \
            ['smearing', 'tetrahedra', 'tetrahedra_lin', 'tetrahedra_opt', 'fixed', 'vc-from_input']:
        raise ValueError('Invalid value for occupations')
    self._nml['system']['occupations'] = str(occupations).lower()


@property
def one_atom_occupations(self):
    return self._nml['system']['one_atom_occupations']


@one_atom_occupations.setter
def one_atom_occupations(self, one_atom_occupations):
    self._nml['system']['one_atom_occupations'] = bool(one_atom_occupations)


@property
def starting_spin_angle(self):
    return self._nml['system']['starting_spin_angle']


@starting_spin_angle.setter
def starting_spin_angle(self, starting_spin_angle):
    self._nml['system']['starting_spin_angle'] = bool(starting_spin_angle)


@property
def degauss(self):
    return self._nml['system']['degauss']


@degauss.setter
def degauss(self, degauss):
    self._nml['system']['degauss'] = float(degauss)


@property
def smearing(self):
    return self._nml['system']['smearing']


@smearing.setter
def smearing(self, smearing):
    if str(smearing).lower() not in \
            ['gaussian', 'gauss', 'methfessel-paxton', 'm-p', 'mp',
                'marzari-vanderbilt', 'cold', 'm-v', 'mv', 'fermi-dirac', 'f-d', 'fd']:
        raise ValueError('Invalid value for smearing')
    self._nml['system']['smearing'] = str(smearing).lower()


@property
def nspin(self):
    return self._nml['system']['nspin']


@nspin.setter
def nspin(self, nspin):
    if int(nspin) not in [1, 2, 4]:
        raise ValueError('Invalid value for nspin')
    self._nml['system']['nspin'] = int(nspin)


@property
def noncolin(self):
    return self._nml['system']['noncolin']


@noncolin.setter
def noncolin(self, noncolin):
    self._nml['system']['noncolin'] = bool(noncolin)


@property
def ecfixed(self):
    return self._nml['system']['ecfixed']


@ecfixed.setter
def ecfixed(self, ecfixed):
    self._nml['system']['ecfixed'] = float(ecfixed)


@property
def qcutz(self):
    return self._nml['system']['qcutz']


@qcutz.setter
def qcutz(self, qcutz):
    self._nml['system']['qcutz'] = float(qcutz)


@property
def q2sigma(self):
    return self._nml['system']['q2sigma']


@q2sigma.setter
def q2sigma(self, q2sigma):
    self._nml['system']['q2sigma'] = float(q2sigma)


@property
def input_dft(self):
    return self._nml['system']['input_dft']


@input_dft.setter
def input_dft(self, input_dft):
    self._nml['system']['input_dft'] = str(input_dft)


@property
def exx_fraction(self):
    return self._nml['system']['exx_fraction']


@exx_fraction.setter
def exx_fraction(self, exx_fraction):
    self._nml['system']['exx_fraction'] = float(exx_fraction)


@property
def screening_parameter(self):
    return self._nml['system']['screening_parameter']


@screening_parameter.setter
def screening_parameter(self, screening_parameter):
    self._nml['system']['screening_parameter'] = float(screening_parameter)


@property
def exxdiv_treatment(self):
    return self._nml['system']['exxdiv_treatment']


@exxdiv_treatment.setter
def exxdiv_treatment(self, exxdiv_treatment):
    if str(exxdiv_treatment).lower() not in ['gygi-baldereschi', 'vcut_spherical', 'vcut_ws', 'none']:
        raise ValueError('Invalid value for exxdiv_treatment')
    self._nml['system']['exxdiv_treatment'] = str(exxdiv_treatment).lower()


@property
def x_gamma_extrapolation(self):
    return self._nml['system']['x_gamma_extrapolation']


@x_gamma_extrapolation.setter
def x_gamma_extrapolation(self, x_gamma_extrapolation):
    self._nml['system']['x_gamma_extrapolation'] = bool(x_gamma_extrapolation)


@property
def ecutvcut(self):
    return self._nml['system']['ecutvcut']


@ecutvcut.setter
def ecutvcut(self, ecutvcut):
    self._nml['system']['ecutvcut'] = float(ecutvcut)


@property
def nqx1(self):
    return self._nml['system']['nqx1']


@nqx1.setter
def nqx1(self, nqx1):
    self._nml['system']['nqx1'] = int(nqx1)


@property
def nqx2(self):
    return self._nml['system']['nqx2']


@nqx2.setter
def nqx2(self, nqx2):
    self._nml['system']['nqx2'] = int(nqx2)


@property
def nqx3(self):
    return self._nml['system']['nqx3']


@nqx3.setter
def nqx3(self, nqx3):
    self._nml['system']['nqx3'] = int(nqx3)


@property
def lda_plus_u(self):
    return self._nml['system']['lda_plus_u']


@lda_plus_u.setter
def lda_plus_u(self, lda_plus_u):
    self._nml['system']['lda_plus_u'] = bool(lda_plus_u)


@property
def lda_plus_u_kind(self):
    return self._nml['system']['lda_plus_u_kind']


@lda_plus_u_kind.setter
def lda_plus_u_kind(self, lda_plus_u_kind):
    self._nml['system']['lda_plus_u_kind'] = int(lda_plus_u_kind)


@property
def Hubbard_U(self):
    return self._nml['system']['Hubbard_U']


@Hubbard_U.setter
def Hubbard_U(self, Hubbard_U):
    self._nml['system']['Hubbard_U'] = [float(h) for h in list(Hubbard_U)]


@property
def Hubbard_J0(self):
    return self._nml['system']['Hubbard_J0']


@Hubbard_J0.setter
def Hubbard_J0(self, Hubbard_J0):
    self._nml['system']['Hubbard_J0'] = [float(h) for h in list(Hubbard_J0)]


@property
def Hubbard_alpha(self):
    return self._nml['system']['Hubbard_alpha']


@Hubbard_alpha.setter
def Hubbard_alpha(self, Hubbard_alpha):
    self._nml['system']['Hubbard_alpha'] = [float(h) for h in list(Hubbard_alpha)]


@property
def Hubbard_beta(self):
    return self._nml['system']['Hubbard_beta']


@Hubbard_beta.setter
def Hubbard_beta(self, Hubbard_beta):
    self._nml['system']['Hubbard_beta'] = [float(h) for h in list(Hubbard_beta)]


@property
def Hubbard_J(self):
    return self._nml['system']['Hubbard_J']


@Hubbard_J.setter
def Hubbard_J(self, Hubbard_J):
    # TODO type not understand leaving this open
    self._nml['system']['Hubbard_J'] = Hubbard_J


@property
def starting_ns_eigenvalue(self):
    return self._nml['system']['starting_ns_eigenvalue']


@starting_ns_eigenvalue.setter
def starting_ns_eigenvalue(self, starting_ns_eigenvalue):
    # TODO type not understand leaving this open
    self._nml['system']['starting_ns_eigenvalue'] = starting_ns_eigenvalue


@property
def U_projection_type(self):
    return self._nml['system']['U_projection_type']


@U_projection_type.setter
def U_projection_type(self, U_projection_type):
    if str(U_projection_type).lower() not in \
            ['atomic', 'ortho-atomic', 'norm-atomic', 'file', 'pseudo']:
        raise ValueError('Invalid value for U_projection_type')
    self._nml['system']['U_projection_type'] = str(U_projection_type).lower()


@property
def edir(self):
    return self._nml['system']['edir']


@edir.setter
def edir(self, edir):
    self._nml['system']['edir'] = int(edir)


@property
def emaxpos(self):
    return self._nml['system']['emaxpos']


@emaxpos.setter
def emaxpos(self, emaxpos):
    self._nml['system']['emaxpos'] = float(emaxpos)


@property
def eopreg(self):
    return self._nml['system']['eopreg']


@eopreg.setter
def eopreg(self, eopreg):
    self._nml['system']['eopreg'] = float(eopreg)


@property
def eamp(self):
    return self._nml['system']['eamp']


@eamp.setter
def eamp(self, eamp):
    self._nml['system']['eamp'] = float(eamp)


@property
def angle1(self):
    return self._nml['system']['angle1']


@angle1.setter
def angle1(self, angle1):
    self._nml['system']['angle1'] = [float(a) for a in list(angle1)]


@property
def angle2(self):
    return self._nml['system']['angle2']


@angle2.setter
def angle2(self, angle2):
    self._nml['system']['angle2'] = [float(a) for a in list(angle2)]


@property
def constrained_magnetization(self):
    return self._nml['system']['constrained_magnetization']


@constrained_magnetization.setter
def constrained_magnetization(self, constrained_magnetization):
    if str(constrained_magnetization).lower() not in \
            ['none', 'total', 'atomic', 'total direction', 'atomic direction']:
        raise ValueError('Invalid value for constrained_magnetization')
    self._nml['system']['constrained_magnetization'] = str(constrained_magnetization).lower()


@property
def fixed_magnetization(self):
    return self._nml['system']['fixed_magnetization']


@fixed_magnetization.setter
def fixed_magnetization(self, fixed_magnetization):
    self._nml['system']['fixed_magnetization'] = \
        [float(fm) for fm in list(fixed_magnetization)[0:3]]


# changed lambda to lambda_qe
@property
def lambda_qe(self):
    return self._nml['system']['lambda']


@lambda_qe.setter
def lambda_qe(self, lambda_qe):
    self._nml['system']['lambda'] = float(lambda_qe)


@property
def report(self):
    return self._nml['system']['report']


@report.setter
def report(self, report):
    self._nml['system']['report'] = int(report)


@property
def lspinorb(self):
    return self._nml['system']['lspinorb']


@lspinorb.setter
def lspinorb(self, lspinorb):
    self._nml['system']['lspinorb'] = bool(lspinorb)


@property
def assume_isolated(self):
    return self._nml['system']['assume_isolated']


@assume_isolated.setter
def assume_isolated(self, assume_isolated):
    if str(assume_isolated).lower() not in \
            ['none', 'makov-payne', 'm-p', 'mp', 'martyna-tuckerman', 'm-t', 'mt', 'esm']:
        raise ValueError('Invalid value for assume_isolated')
    self._nml['system']['assume_isolated'] = str(assume_isolated).lower()


@property
def esm_bc(self):
    return self._nml['system']['esm_bc']


@esm_bc.setter
def esm_bc(self, esm_bc):
    if str(esm_bc).lower() not in \
            ['pbc', 'bc1', 'bc2', 'bc3']:
        raise ValueError('Invalid value for esm_bc')
    self._nml['system']['esm_bc'] = str(esm_bc).lower()


@property
def esm_w(self):
    return self._nml['system']['esm_w']


@esm_w.setter
def esm_w(self, esm_w):
    self._nml['system']['esm_w'] = float(esm_w)


@property
def esm_efield(self):
    return self._nml['system']['esm_efield']


@esm_efield.setter
def esm_efield(self, esm_efield):
    self._nml['system']['esm_efield'] = float(esm_efield)


@property
def esm_nfit(self):
    return self._nml['system']['esm_nfit']


@esm_nfit.setter
def esm_nfit(self, esm_nfit):
    self._nml['system']['esm_nfit'] = int(esm_nfit)


@property
def fcp_mu(self):
    return self._nml['system']['fcp_mu']


@fcp_mu.setter
def fcp_mu(self, fcp_mu):
    self._nml['system']['fcp_mu'] = float(fcp_mu)


@property
def vdw_corr(self):
    return self._nml['system']['vdw_corr']


@vdw_corr.setter
def vdw_corr(self, vdw_corr):
    if str(vdw_corr).lower() not in \
            ['grimme-d2', 'dft-d', 'ts', 'ts-vdw', 'tkatchenko-scheffler', 'xdm']:
        raise ValueError('Invalid value for vdw_corr')
    self._nml['system']['vdw_corr'] = str(vdw_corr).lower()


# Obsolete
# @property
# def london(self):
#     return self._nml['system']['london']


# @london.setter
# def london(self, london):
#     self._nml['system']['london'] = london


@property
def london_s6(self):
    return self._nml['system']['london_s6']


@london_s6.setter
def london_s6(self, london_s6):
    self._nml['system']['london_s6'] = float(london_s6)


@property
def london_c6(self):
    return self._nml['system']['london_c6']


@london_c6.setter
def london_c6(self, london_c6):
    self._nml['system']['london_c6'] = float(london_c6)


@property
def london_rvdw(self):
    return self._nml['system']['london_rvdw']


@london_rvdw.setter
def london_rvdw(self, london_rvdw):
    self._nml['system']['london_rvdw'] = float(london_rvdw)


@property
def london_rcut(self):
    return self._nml['system']['london_rcut']


@london_rcut.setter
def london_rcut(self, london_rcut):
    self._nml['system']['london_rcut'] = float(london_rcut)


@property
def ts_vdw_econv_thr(self):
    return self._nml['system']['ts_vdw_econv_thr']


@ts_vdw_econv_thr.setter
def ts_vdw_econv_thr(self, ts_vdw_econv_thr):
    self._nml['system']['ts_vdw_econv_thr'] = float(ts_vdw_econv_thr)


@property
def ts_vdw_isolated(self):
    return self._nml['system']['ts_vdw_isolated']


@ts_vdw_isolated.setter
def ts_vdw_isolated(self, ts_vdw_isolated):
    self._nml['system']['ts_vdw_isolated'] = bool(ts_vdw_isolated)


# Obsolete
# @property
# def xdm(self):
#     return self._nml['system']['xdm']


# @xdm.setter
# def xdm(self, xdm):
#     self._nml['system']['xdm'] = xdm


@property
def xdm_a1(self):
    return self._nml['system']['xdm_a1']


@xdm_a1.setter
def xdm_a1(self, xdm_a1):
    self._nml['system']['xdm_a1'] = float(xdm_a1)


@property
def xdm_a2(self):
    return self._nml['system']['xdm_a2']


@xdm_a2.setter
def xdm_a2(self, xdm_a2):
    self._nml['system']['xdm_a2'] = float(xdm_a2)


@property
def space_group(self):
    return self._nml['system']['space_group']


@space_group.setter
def space_group(self, space_group):
    self._nml['system']['space_group'] = int(space_group)


@property
def uniqueb(self):
    return self._nml['system']['uniqueb']


@uniqueb.setter
def uniqueb(self, uniqueb):
    self._nml['system']['uniqueb'] = bool(uniqueb)


@property
def origin_choice(self):
    return self._nml['system']['origin_choice']


@origin_choice.setter
def origin_choice(self, origin_choice):
    self._nml['system']['origin_choice'] = int(origin_choice)


@property
def rhombohedral(self):
    return self._nml['system']['rhombohedral']


@rhombohedral.setter
def rhombohedral(self, rhombohedral):
    self._nml['system']['rhombohedral'] = bool(rhombohedral)


@property
def zmon(self):
    return self._nml['system']['zmon']


@zmon.setter
def zmon(self, zmon):
    self._nml['system']['zmon'] = float(zmon)


@property
def realxz(self):
    return self._nml['system']['realxz']


@realxz.setter
def realxz(self, realxz):
    self._nml['system']['realxz'] = bool(realxz)


@property
def block(self):
    return self._nml['system']['block']


@block.setter
def block(self, block):
    self._nml['system']['block'] = bool(block)


@property
def block_1(self):
    return self._nml['system']['block_1']


@block_1.setter
def block_1(self, block_1):
    self._nml['system']['block_1'] = float(block_1)


@property
def block_2(self):
    return self._nml['system']['block_2']


@block_2.setter
def block_2(self, block_2):
    self._nml['system']['block_2'] = float(block_2)


@property
def block_height(self):
    return self._nml['system']['block_height']


@block_height.setter
def block_height(self, block_height):
    self._nml['system']['block_height'] = float(block_height)


'''
######################################################################################
################################## &electrons ########################################
######################################################################################
'''


@property
def electron_maxstep(self):
    return self._nml['electrons']['electron_maxstep']


@electron_maxstep.setter
def electron_maxstep(self, electron_maxstep):
    self._nml['electrons']['electron_maxstep'] = int(electron_maxstep)


@property
def scf_must_converge(self):
    return self._nml['electrons']['scf_must_converge']


@scf_must_converge.setter
def scf_must_converge(self, scf_must_converge):
    self._nml['electrons']['scf_must_converge'] = bool(scf_must_converge)


@property
def conv_thr(self):
    return self._nml['electrons']['conv_thr']


@conv_thr.setter
def conv_thr(self, conv_thr):
    self._nml['electrons']['conv_thr'] = float(conv_thr)


@property
def adaptive_thr(self):
    return self._nml['electrons']['adaptive_thr']


@adaptive_thr.setter
def adaptive_thr(self, adaptive_thr):
    self._nml['electrons']['adaptive_thr'] = bool(adaptive_thr)


@property
def conv_thr_init(self):
    return self._nml['electrons']['conv_thr_init']


@conv_thr_init.setter
def conv_thr_init(self, conv_thr_init):
    self._nml['electrons']['conv_thr_init'] = float(conv_thr_init)


@property
def conv_thr_multi(self):
    return self._nml['electrons']['conv_thr_multi']


@conv_thr_multi.setter
def conv_thr_multi(self, conv_thr_multi):
    self._nml['electrons']['conv_thr_multi'] = float(conv_thr_multi)


@property
def mixing_mode(self):
    return self._nml['electrons']['mixing_mode']


@mixing_mode.setter
def mixing_mode(self, mixing_mode):
    if str(mixing_mode).lower() == 'plain':
        self._nml['system']['mixing_mode'] = str(mixing_mode).lower()
    elif str(mixing_mode) in ['local-TF', 'TF']:
        # don't lower, these are not lower-cased
        self._nml['system']['mixing_mode'] = str(mixing_mode)
    else:
        raise ValueError('Invalid value for mixing_mode')


@property
def mixing_beta(self):
    return self._nml['electrons']['mixing_beta']


@mixing_beta.setter
def mixing_beta(self, mixing_beta):
    self._nml['electrons']['mixing_beta'] = float(mixing_beta)


@property
def mixing_ndim(self):
    return self._nml['electrons']['mixing_ndim']


@mixing_ndim.setter
def mixing_ndim(self, mixing_ndim):
    self._nml['electrons']['mixing_ndim'] = int(mixing_ndim)


@property
def mixing_fixed_ns(self):
    return self._nml['electrons']['mixing_fixed_ns']


@mixing_fixed_ns.setter
def mixing_fixed_ns(self, mixing_fixed_ns):
    self._nml['electrons']['mixing_fixed_ns'] = int(mixing_fixed_ns)


@property
def diagonalization(self):
    return self._nml['electrons']['diagonalization']


@diagonalization.setter
def diagonalization(self, diagonalization):
    if str(diagonalization).lower() not in ['david', 'cg', 'david-serial', 'cg-serial']:
        raise ValueError('Invalid value for diagonalization')
    self._nml['system']['diagonalization'] = str(diagonalization).lower()


# Obsolete
# @property
# def ortho_para(self):
#     return self._nml['electrons']['ortho_para']


# @ortho_para.setter
# def ortho_para(self, ortho_para):
#     self._nml['electrons']['ortho_para'] = int(ortho_para)


@property
def diago_thr_init(self):
    return self._nml['electrons']['diago_thr_init']


@diago_thr_init.setter
def diago_thr_init(self, diago_thr_init):
    self._nml['electrons']['diago_thr_init'] = float(diago_thr_init)


@property
def diago_cg_maxiter(self):
    return self._nml['electrons']['diago_cg_maxiter']


@diago_cg_maxiter.setter
def diago_cg_maxiter(self, diago_cg_maxiter):
    self._nml['electrons']['diago_cg_maxiter'] = int(diago_cg_maxiter)


@property
def diago_david_ndim(self):
    return self._nml['electrons']['diago_david_ndim']


@diago_david_ndim.setter
def diago_david_ndim(self, diago_david_ndim):
    self._nml['electrons']['diago_david_ndim'] = int(diago_david_ndim)


@property
def diago_full_acc(self):
    return self._nml['electrons']['diago_full_acc']


@diago_full_acc.setter
def diago_full_acc(self, diago_full_acc):
    self._nml['electrons']['diago_full_acc'] = bool(diago_full_acc)


@property
def efield(self):
    return self._nml['electrons']['efield']


@efield.setter
def efield(self, efield):
    self._nml['electrons']['efield'] = float(efield)


@property
def efield_cart(self):
    return self._nml['electrons']['efield_cart']


@efield_cart.setter
def efield_cart(self, efield_cart):
    self._nml['electrons']['efield_cart'] = [float(ec) for ec in list(efield_cart)[0:3]]


@property
def efield_phase(self):
    return self._nml['electrons']['efield_phase']


@efield_phase.setter
def efield_phase(self, efield_phase):
    if str(efield_phase).lower() not in ['read', 'write', 'none']:
        raise ValueError('Invalid value for efield_phase')
    self._nml['system']['efield_phase'] = str(efield_phase).lower()


@property
def startingpot(self):
    return self._nml['electrons']['startingpot']


@startingpot.setter
def startingpot(self, startingpot):
    if str(startingpot).lower() not in ['atomic', 'file']:
        raise ValueError('Invalid value for startingpot')
    self._nml['system']['startingpot'] = str(startingpot).lower()


@property
def startingwfc(self):
    return self._nml['electrons']['startingwfc']


@startingwfc.setter
def startingwfc(self, startingwfc):
    if str(startingwfc).lower() not in ['atomic', 'atomic+random', 'random', 'file']:
        raise ValueError('Invalid value for startingwfc')
    self._nml['system']['startingwfc'] = str(startingwfc).lower()


@property
def tqr(self):
    return self._nml['electrons']['tqr']


@tqr.setter
def tqr(self, tqr):
    self._nml['electrons']['tqr'] = bool(tqr)


'''
######################################################################################
##################################### &ions ##########################################
######################################################################################
'''


@property
def ion_dynamics(self):
    return self._nml['ions']['ion_dynamics']


@ion_dynamics.setter
def ion_dynamics(self, ion_dynamics):
    if str(ion_dynamics).lower() not in ['bfgs', 'damp', 'verlet', 'langevin',
                                         'langevin-smc', 'beeman']:
        raise ValueError('Invalid value for ion_dynamics')
    self._nml['system']['ion_dynamics'] = str(ion_dynamics).lower()


@property
def ion_positions(self):
    return self._nml['ions']['ion_positions']


@ion_positions.setter
def ion_positions(self, ion_positions):
    if str(ion_positions).lower() not in ['default', 'from_input']:
        raise ValueError('Invalid value for ion_positions')
    self._nml['system']['ion_positions'] = str(ion_positions).lower()


@property
def pot_extrapolation(self):
    return self._nml['ions']['pot_extrapolation']


@pot_extrapolation.setter
def pot_extrapolation(self, pot_extrapolation):
    if str(pot_extrapolation).lower() not in ['none', 'atomic', 'first_order', 'second_order']:
        raise ValueError('Invalid value for pot_extrapolation')
    self._nml['system']['pot_extrapolation'] = str(pot_extrapolation).lower()


@property
def wfc_extrapolation(self):
    return self._nml['ions']['wfc_extrapolation']


@wfc_extrapolation.setter
def wfc_extrapolation(self, wfc_extrapolation):
    if str(wfc_extrapolation).lower() not in ['none', 'first_order', 'second_order']:
        raise ValueError('Invalid value for wfc_extrapolation')
    self._nml['system']['wfc_extrapolation'] = str(wfc_extrapolation).lower()


@property
def remove_rigid_rot(self):
    return self._nml['ions']['remove_rigid_rot']


@remove_rigid_rot.setter
def remove_rigid_rot(self, remove_rigid_rot):
    self._nml['ions']['remove_rigid_rot'] = bool(remove_rigid_rot)


@property
def ion_temperature(self):
    return self._nml['ions']['ion_temperature']


@ion_temperature.setter
def ion_temperature(self, ion_temperature):
    if str(ion_temperature).lower() in ['rescaling', 'rescaling-v', 'berendsen',
                                        'andersen', 'initial', 'not_controlled']:
        self._nml['system']['ion_temperature'] = str(ion_temperature).lower()
    elif str(ion_temperature) in ['rescaling-T', 'reduce-T']:
        # don't lower, these are not lower-cased
        self._nml['system']['ion_temperature'] = str(ion_temperature)
    else:
        raise ValueError('Invalid value for ion_temperature')


@property
def tempw(self):
    return self._nml['ions']['tempw']


@tempw.setter
def tempw(self, tempw):
    self._nml['ions']['tempw'] = float(tempw)


@property
def tolp(self):
    return self._nml['ions']['tolp']


@tolp.setter
def tolp(self, tolp):
    self._nml['ions']['tolp'] = float(tolp)


@property
def delta_t(self):
    return self._nml['ions']['delta_t']


@delta_t.setter
def delta_t(self, delta_t):
    self._nml['ions']['delta_t'] = float(delta_t)


@property
def nraise(self):
    return self._nml['ions']['nraise']


@nraise.setter
def nraise(self, nraise):
    self._nml['ions']['nraise'] = int(nraise)


@property
def refold_pos(self):
    return self._nml['ions']['refold_pos']


@refold_pos.setter
def refold_pos(self, refold_pos):
    self._nml['ions']['refold_pos'] = bool(refold_pos)


@property
def upscale(self):
    return self._nml['ions']['upscale']


@upscale.setter
def upscale(self, upscale):
    self._nml['ions']['upscale'] = float(upscale)


@property
def bfgs_ndim(self):
    return self._nml['ions']['bfgs_ndim']


@bfgs_ndim.setter
def bfgs_ndim(self, bfgs_ndim):
    self._nml['ions']['bfgs_ndim'] = int(bfgs_ndim)


@property
def trust_radius_max(self):
    return self._nml['ions']['trust_radius_max']


@trust_radius_max.setter
def trust_radius_max(self, trust_radius_max):
    self._nml['ions']['trust_radius_max'] = float(trust_radius_max)


@property
def trust_radius_min(self):
    return self._nml['ions']['trust_radius_min']


@trust_radius_min.setter
def trust_radius_min(self, trust_radius_min):
    self._nml['ions']['trust_radius_min'] = float(trust_radius_min)


@property
def trust_radius_ini(self):
    return self._nml['ions']['trust_radius_ini']


@trust_radius_ini.setter
def trust_radius_ini(self, trust_radius_ini):
    self._nml['ions']['trust_radius_ini'] = float(trust_radius_ini)


@property
def w_1(self):
    return self._nml['ions']['w_1']


@w_1.setter
def w_1(self, w_1):
    self._nml['ions']['w_1'] = float(w_1)


@property
def w_2(self):
    return self._nml['ions']['w_2']


@w_2.setter
def w_2(self, w_2):
    self._nml['ions']['w_2'] = float(w_2)


'''
######################################################################################
##################################### &cell ##########################################
######################################################################################
'''


@property
def cell_dynamics(self):
    return self._nml['cell']['cell_dynamics']


@cell_dynamics.setter
def cell_dynamics(self, cell_dynamics):
    if str(cell_dynamics).lower() not in ['none', 'sd', 'damp-pr', 'damp-w', 'bfgs', 'pr', 'w']:
        raise ValueError('Invalid value for cell_dynamics')
    self._nml['system']['cell_dynamics'] = str(cell_dynamics).lower()


@property
def press(self):
    return self._nml['cell']['press']


@press.setter
def press(self, press):
    self._nml['cell']['press'] = float(press)


@property
def wmass(self):
    return self._nml['cell']['wmass']


@wmass.setter
def wmass(self, wmass):
    self._nml['cell']['wmass'] = float(wmass)


@property
def cell_factor(self):
    return self._nml['cell']['cell_factor']


@cell_factor.setter
def cell_factor(self, cell_factor):
    self._nml['cell']['cell_factor'] = float(cell_factor)


@property
def press_conv_thr(self):
    return self._nml['cell']['press_conv_thr']


@press_conv_thr.setter
def press_conv_thr(self, press_conv_thr):
    self._nml['cell']['press_conv_thr'] = float(press_conv_thr)


@property
def cell_dofree(self):
    return self._nml['cell']['cell_dofree']


@cell_dofree.setter
def cell_dofree(self, cell_dofree):
    if str(cell_dofree).lower() not in ['all', 'x', 'y', 'z',
                                        'xy', 'xz', 'yz', 'xyz', 'shape', 'volume', '2Dxy', '2Dshape']:
        raise ValueError('Invalid value for cell_dofree')
    self._nml['system']['cell_dofree'] = str(cell_dofree).lower()
