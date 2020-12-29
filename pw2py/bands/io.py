import numpy as np
from lxml import etree
import os.path

from ..qesave import qesave
from .._common.constants import ha_to_ev


def _read_qeout_data(f, prec: float = 9):
    '''
    designed read block of data in qe output file e.g. bands and occupations
    '''
    data = []
    while True:
        fline = f.readline()
        if len(fline.strip()) == 0:
            break
        # remove two spaces at beginning, new line character at end
        fline = fline[2:-1]
        # manually split fline
        fline_split = []
        for i in range(len(fline) // prec):
            # grab next prec characters in fline
            value = fline[prec*i:prec*(i+1)]
            try:
                float(value)
                fline_split.append(value)
            except ValueError:
                fline_split.append(np.nan)
        data += fline_split
    # return np.ndarray(data, dtype=dtype)
    return data


def _read_qeout_bands(filename: str):
    '''
    read nspin, occupations, fermi, etc. from qesave
    '''
    nspin = False
    occ_given = False
    fermi_given = False
    vbm_given = False
    cbm_given = False
    with open(filename) as f:
        for line in f:
            if line.startswith('     number of electrons'):
                if 'up' not in line:
                    nelec = [int(float(line.split('=')[-1]))]
                else:
                    # read up/down electrons
                    csplit = line.split(':')
                    nelec = [
                        int(float(csplit[1].split(',')[0])),
                        int(float(csplit[-1].split(')')[0]))
                    ]
            elif line.startswith('     number of k points='):
                nk = int(line.split('=')[-1])
                # skip a line
                f.readline()
                # read weights
                weight = [f.readline().split('=')[-1] for _ in range(nk)]
            elif 'End of self-consistent calculation' in line:
                # preparation to read kpt, eig, occ
                kpt = []
                read_eig = []
                read_occ = []
            elif 'SPIN UP' in line:
                # preparation for reading spin up (read_eig points to eig_up)
                nspin = True
                eig_up = []
                read_eig = eig_up
                occ_up = []
                read_occ = occ_up
            elif 'SPIN DOWN' in line:
                # preparation for reading spin down (read_eig points to eig_dn)
                eig_dn = []
                read_eig = eig_dn
                occ_dn = []
                read_occ = occ_dn
            elif line.startswith('          k ='):
                # read the kpt
                kptstr = line.split('=')[-1].split('(')[0][:-1]
                k = np.array(
                    [kptstr[7*i:7*(i+1)] for i in range(3)], dtype=float
                )
                kpt.append(k)
                # skip a line
                f.readline()
                # read the eigenvalues
                read_eig.append(_read_qeout_data(f))
            elif line.startswith('     occupation numbers '):
                occ_given = True
                read_occ.append(_read_qeout_data(f))
            elif line.startswith('     highest occupied level (ev): '):
                vbm_given = True
                vbm = float(line.split()[-1])
            elif line.startswith('     highest occupied, lowest unoccupied level (ev): '):
                vbm_given, cbm_given = True, True
                data = line[53:]
                vbm, cbm = [float(data[10*i:10*(i+1)]) for i in range(2)]
            elif line.startswith('     the Fermi energy is '):
                fermi_given = True
                fermi = float(line.split()[4])
            elif line.startswith('     the spin up/dw Fermi energies are '):
                fermi_given = True
                data = line[39:]
                fermi = np.array(
                    [float(data[10*i:10*(i+1)]) for i in range(2)]
                )
                if np.isclose(fermi[0], fermi[1], atol=1e-4):
                    fermi = fermi[0]
    # ------------------------------------------------------------
    # collect values read into dictionary 'coll'
    coll = {}
    coll['nspin'] = 2 if nspin else 1
    coll['kpt'] = np.array(kpt, dtype=float).reshape(
        coll['nspin'], len(kpt) // coll['nspin'], 3
    )
    coll['kpt'] = coll['kpt'][0]  # in the end just up
    coll['weight'] = np.array(weight, dtype=float)
    if fermi_given:
        coll['fermi'] = fermi
    if vbm_given:
        coll['vbm'] = vbm
    if cbm_given:
        coll['cbm'] = cbm
    if nspin:
        coll['eig'] = np.array([eig_up, eig_dn], dtype=float)
        if occ_given:
            coll['occ'] = np.array([occ_up, occ_dn], dtype=float)
    else:
        coll['eig'] = np.array([read_eig], dtype=float)
        if occ_given:
            coll['occ'] = np.array([read_occ], dtype=float)
    if not occ_given:
        if (coll['nspin'] == 1) and (nelec[0]//2 == coll['eig'].shape[2]):
            # number of electrons is the same as the number of bands, all states are occupied
            coll['occ'] = np.ones(coll['eig'].shape)
        elif (coll['nspin'] == 2):
            # for the case where len(nelec) != nspin, there is an assumption here that tot_mag = 0
            # if tot_mag != 0 then nelec should have been given for up and down and hence len(nelec) != nspin
            # otherwise occupations should have been given. So this should rarely be needed... maybe never
            nelec = nelec if len(nelec) == coll['nspin'] \
                else [nelec[0]//2, nelec[0]//2]
            if (nelec[0] == coll['eig'][0].shape[1]) and (nelec[1] == coll['eig'][1].shape[1]):
                coll['occ'] = np.ones(coll['eig'].shape)
    return coll


def _read_qesave_bands_61(save: qesave):
    '''
    read bands from xml file in qesave version 6.1
    '''
    # ------------------------------------------------------------
    # read nspin, nkpt, nbnd
    bands_child = save.root.find('BAND_STRUCTURE_INFO')
    nspin = int(bands_child.find('NUMBER_OF_SPIN_COMPONENTS').text)
    nk = int(bands_child.find('NUMBER_OF_K-POINTS').text)
    # TODO read fermi
    # ------------------------------------------------------------
    # read kpt, weight
    bz_child = save.root.find('BRILLOUIN_ZONE')
    kpt, weight = [], []
    for ik in range(nk):
        key = f'K-POINT.{ik+1}'
        kpt_child = bz_child.find(key)
        kpt.append(kpt_child.attrib['XYZ'].split())
        weight.append(kpt_child.attrib['WEIGHT'])
    # ------------------------------------------------------------
    # read eig, occ
    eig, occ = [], []
    for ispin in range(nspin):
        eig.append([])
        occ.append([])
        basename = 'eigenval{}.xml'.format(ispin + 1) if nspin == 2 else \
            'eigenval.xml'
        for ik in range(nk):
            kdir = f'K{ik+1:05d}'
            eigfile = os.path.join(save.savefolder, kdir, basename)
            eigroot = etree.parse(eigfile).getroot()
            for dat, key in ((eig, 'EIGENVALUES'), (occ, 'OCCUPATIONS')):
                dat[-1].append(eigroot.find(key).text.split())
    # ------------------------------------------------------------
    # collect values read into dictionary 'coll'
    coll = {}
    coll['kpt'] = np.array(kpt, dtype=float).reshape(nk, 3)
    coll['weight'] = np.array(weight, dtype=float)
    coll['eig'] = np.array(eig, dtype=float) * ha_to_ev
    coll['occ'] = np.array(occ, dtype=float)
    return coll


def _read_qesave_bands_66(save: qesave):
    '''
    read bands from xml file in qesave version 6.6
    '''
    band_child = save.root.find('output').find('band_structure')
    lsda = band_child.find('lsda').text == 'true'
    nspin = 2 if lsda else 1
    bnd_str = 'nbnd_up' if lsda else 'nbnd'
    nbnd = int(band_child.find(bnd_str).text)
    if lsda:
        nbnd_dw = int(band_child.find('nbnd_dw').text)
        assert nbnd == nbnd_dw, "Cannot read bands when up/down nbnd differ"
    # nk = int(band_child.find('nks').text)
    eig, occ, kpt, weight = [], [], [], []
    for ks_child in band_child.findall('ks_energies'):
        # ------------------------------------------------------------
        # read kpt and weight
        kpt_child = ks_child.find('k_point')
        kpt.append(kpt_child.text.split())
        weight.append(kpt_child.attrib['weight'])
        # ------------------------------------------------------------
        # read eigenvalues and occupations
        for dat, key in ((eig, 'eigenvalues'), (occ, 'occupations')):
            dat_child = ks_child.find(key)
            dat.append(dat_child.text.split())
    nk = len(kpt)
    # ------------------------------------------------------------
    # convert to arrays
    kpt = np.array(kpt, dtype=np.float64)
    weight = np.array(weight, dtype=np.float64)
    eig = np.array(eig, dtype=np.float64).reshape(nk, nspin, nbnd) * ha_to_ev
    occ = np.array(occ, dtype=np.float64).reshape(nk, nspin, nbnd)
    # ------------------------------------------------------------
    # reorder eigenvalues and occupations (spin, kpt, band)
    eig = np.swapaxes(eig, 0, 1)
    occ = np.swapaxes(occ, 0, 1)
    # ------------------------------------------------------------
    # collect values read into dictionary 'coll'
    coll = {}
    coll['kpt'] = kpt
    coll['weight'] = weight
    coll['eig'] = eig
    coll['occ'] = occ
    return coll


@classmethod
def from_save(cls, prefix: str, path: str = '.', version: str = None):
    '''
    read bands from qe savefolder
    '''
    save = qesave(prefix=prefix, path=path)
    version = save.version if version is None else version
    if version == '6.1':
        return cls(**_read_qesave_bands_61(save))
    else:
        return cls(**_read_qesave_bands_66(save))


@classmethod
def from_file(cls, filename: str):
    '''
    read bands from file, creating Bands object

    params
    ---
        filename (str)

    returns
    ---
        Bands

    Note
    ---
        eig will contain np.nan if eigenvalue overflowed, i.e. exceeded
        9 characters (less than -999.9999) or (greater than 9999.9999)
    '''
    return cls(**_read_qeout_bands(filename))  # pylint: disable=E1102
