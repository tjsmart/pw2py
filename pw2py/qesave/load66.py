import numpy as np
import mmap
import os.path
import glob
from lxml import etree
from scipy.constants import physical_constants

Ha2eV = physical_constants['Hartree energy in eV'][0]


def read_wfc_file(filename: str):
    '''
    read wfc file (dat format, not hdf5)
    '''
    # see Modules/io_base.f90 from QE src for more details
    with open(filename, 'rb') as f:
        with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as buffer:
            # ------------------------------------------------------------
            # read general information
            buffer.read(4)
            ik = np.frombuffer(buffer.read(4), dtype=np.int32)[0]
            kvec = np.frombuffer(buffer.read(24), dtype=np.float64)
            ispin = np.frombuffer(buffer.read(4), dtype=np.int32)[0]
            gamma_only = np.frombuffer(buffer.read(4), dtype=np.bool)[0]
            scalef = np.frombuffer(buffer.read(8), dtype=np.float64)[0]
            info = {'ik': ik, 'kvec': kvec, 'ispin': ispin,  # noqa
                    'gamma_only': gamma_only, 'scalef': scalef}
            buffer.read(8)  # newline
            # ------------------------------------------------------------
            # read npw, npol, nbnd
            ngw, igwx, npol, nbnd = np.frombuffer(buffer.read(16), dtype=np.int32)  # noqa
            buffer.read(8)  # newline
            # ------------------------------------------------------------
            # read reciprocal lattice vectors
            b = np.frombuffer(buffer.read(72), dtype=np.float64).reshape(3, 3)  # noqa
            buffer.read(8)  # newline
            # ------------------------------------------------------------
            # read gvecs
            gksize = 12 * igwx  # 3 * int32 = 12 bytes
            gk = np.frombuffer(buffer.read(gksize), dtype=np.int32)
            gk = gk.reshape(gk.size//3, 3)
            buffer.read(8)  # newline
            # ------------------------------------------------------------
            # read wfc
            if gamma_only:
                dtype = np.float64
                dsize = 8  # float64 = 8 bytes
            else:
                dtype = np.complex128
                dsize = 16  # complex128 = 16 bytes
            wfc = []
            wtmpsize = dsize * igwx * npol
            for i in range(nbnd):
                wtmp = np.frombuffer(buffer.read(wtmpsize), dtype=dtype)
                wfc.append(wtmp)
                buffer.read(8)  # newline
    #
    return gk, wfc


def read_wavefunction(path: str):
    '''
    read evc from qe save folder (also returns gk)

    in qe6.6 evc and gk are stored togehter so this returns gk as well

    input
    ----
        path - path to .save folder from qe calculation (str)

    returns
    ----
        tuple: (gk, evc)
            gk - g-vectors for each k-point, ndim = 3, indexed by [kpt, pw, (x, y, z)] (array)
            evc - all wavefunctions psi(G), ndim = 4, indexed by [kpt, spin, band, pw] (array)
    '''
    # ------------------------------------------------------------
    # check for file and determine nspin
    file1 = os.path.join(path, 'wfc1.dat')
    file2 = os.path.join(path, 'wfcup1.dat')
    if os.path.exists(file1):
        nspin = 1
    elif os.path.exists(file2):
        nspin = 2
    else:
        raise FileNotFoundError(
            f"Unable to locate wfc*.dat file under: {path}")
    # ------------------------------------------------------------
    # check for number of k-points
    nk = 0
    pre_length = 3 if nspin == 1 else 5
    for filename in glob.glob(os.path.join(path, 'wfc*.dat')):
        f = filename.split('/')[-1]
        ik = int(f[pre_length:-4])
        if ik > nk:
            nk = ik
    # ------------------------------------------------------------
    # loop over ispin and ik reading each dat file
    evc, gk = [], []
    for ik in range(nk):
        evc.append([])
        for ispin in range(nspin):
            if ispin == 0:
                sstr = 'up' if nspin == 2 else ''
            elif ispin == 1:
                sstr = 'dw'
            filename = os.path.join(path, f'wfc{sstr}{ik+1}.dat')
            gk_k, evc_k = read_wfc_file(filename)
            evc[ik].append(evc_k)
            if ispin == 0:
                gk.append(gk_k)
    return gk, evc


def _check_gamma_only(path: str) -> bool:
    xmlfile = os.path.join(path, 'data-file-schema.xml')
    root = etree.parse(xmlfile).getroot()
    return root.find('input').find('basis').find('gamma_only').text == 'true'


def read_eigenvalues(path: str, units='Ha'):
    '''
    Read eigenvalues from save folder in qe-6.6

    default units is Ha (Hartree)

    return
    ---
        tuple (eigenvalues, occupations)
            each array are indexed over: ik, ispin, ib
    '''
    assert units.lower() in ['ha', 'ev', 'ry']
    xmlfile = os.path.join(path, 'data-file-schema.xml')
    root = etree.parse(xmlfile).getroot()
    band_child = root.find('output').find('band_structure')
    lsda = band_child.find('lsda').text == 'true'
    nspin = 2 if lsda else 1
    bnd_str = 'nbnd_up' if lsda else 'nbnd'
    nbnd = int(band_child.find(bnd_str).text)
    # nk = int(band_child.find('nks').text)
    eig, occ = [], []
    for ks_child in band_child.findall('ks_energies'):
        keig = np.fromstring(ks_child.find('eigenvalues').text, sep=' ')
        kocc = np.fromstring(ks_child.find('occupations').text, sep=' ')
        eig.append(keig.reshape(nspin, nbnd))
        occ.append(kocc.reshape(nspin, nbnd))
    eig = np.array(eig, dtype=np.float64)
    if units.lower() == 'ev':
        eig *= Ha2eV
    elif units.lower() == 'ry':
        eig *= 2
    occ = np.array(occ, dtype=np.float64)
    return eig, occ


def read_kvecs(path: str):
    '''
    Read k vectors from save folder in qe-6.6

    return
    ---
        list of k-vectors. Each k-vector is an ndarray of shape (3,)
    '''
    xmlfile = os.path.join(path, 'data-file-schema.xml')
    root = etree.parse(xmlfile).getroot()
    band_child = root.find('output').find('band_structure')
    return [
        np.fromstring(ks_child.find('k_point').text, sep=' ', dtype=np.float64)
        for ks_child in band_child.findall('ks_energies')
    ]


def read_rhog(path: str):
    '''
    read rhog from save folder in qe-6.6 (dat format, not hdf5)

    charge density in G-space

    if nspin == 2, spin density is also read
    '''
    # see Modules/io_base.f90 from QE src for more details
    rhofile = os.path.join(path, 'charge-density.dat')
    with open(rhofile, 'rb') as f:
        with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as buffer:
            # ------------------------------------------------------------
            # read various info
            buffer.read(4)
            gamma_only = np.frombuffer(buffer.read(4), dtype=np.bool)[0]
            ngm_g = np.frombuffer(buffer.read(4), dtype=np.int32)[0]
            nspin = np.frombuffer(buffer.read(4), dtype=np.int32)[0]
            # ------------------------------------------------------------
            # read reciprocal lattice vectors
            buffer.read(8)
            b = np.frombuffer(buffer.read(72), dtype=np.float64).reshape(3, 3)  # noqa
            # ------------------------------------------------------------
            # read g vectors
            buffer.read(8)
            gsize = 12 * ngm_g  # 3 int32's = 12 bytes
            g = np.frombuffer(buffer.read(gsize), dtype=np.int32)
            g = g.reshape(g.size//3, 3)
            # ------------------------------------------------------------
            # read wfc
            buffer.read(8)
            rho = []
            if gamma_only:
                rtmpsize = 8 * ngm_g
                dtype = np.float64
            else:
                rtmpsize = 16 * ngm_g
                dtype = np.complex128
            for ispin in range(nspin):
                rtmp = np.frombuffer(buffer.read(rtmpsize), dtype)
                rho.append(rtmp)
                buffer.read(8)
    return g, rho
