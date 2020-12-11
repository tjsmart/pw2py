import numpy as np

from ..qeout.methods import _read_qeout_data


def _read_bands(filename: str):
    nspin = False
    occ_given = False
    fermi_given = False
    vbm_given = False
    cbm_given = False
    with open(filename) as f:
        for line in f:
            if 'End of self-consistent calculation' in line:
                kpt = []
                read_eig = []
                read_occ = []
            elif 'SPIN UP' in line:
                nspin = True
                eig_up = []
                read_eig = eig_up
                occ_up = []
                read_occ = occ_up
            elif 'SPIN DOWN' in line:
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

    # collect values read into dictionary 'collection'
    collection = {}
    collection['nspin'] = 2 if nspin else 1
    collection['kpt'] = np.array(kpt, dtype=float).reshape(
        collection['nspin'], len(kpt) // collection['nspin'], 3
    )
    if fermi_given:
        collection['fermi'] = fermi
    if vbm_given:
        collection['vbm'] = vbm
    if cbm_given:
        collection['cbm'] = cbm
    if nspin:
        collection['eig'] = np.array([eig_up, eig_dn], dtype=float)
        if occ_given:
            collection['occ'] = np.array([occ_up, occ_dn], dtype=float)
    else:
        collection['eig'] = np.array([read_eig], dtype=float)
        if occ_given:
            collection['occ'] = np.array([read_occ], dtype=float)
    return collection


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
    return cls(**_read_bands(filename))  # pylint: disable=E1102
