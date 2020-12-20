import numpy as np
import mmap
import os.path
import glob
from lxml import etree


def read_wfc_file(filename: str, isreal: bool = False):
    '''
    read wfc filename (dat format, not hdf5)

    WARNING!!! This should only be used for reading kp=(0, 0, 0), others are not tested

    WARNING!!! gamma_only is not implemented (should read real wfc instead of complex)

    WARNING!!! npol != 1 is not implemented (e.g. noncollinear case)
    '''
    if isreal:
        dtype = np.float64
        dsize = 8  # float64 = 8 bytes
    else:
        dtype = np.complex128
        dsize = 16  # complex128 = 16 bytes

    with open(filename, 'rb') as f:
        with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as buffer:
            # ik_ = np.frombuffer(buffer.read(4), dtype=np.int32)[0]
            # xk = np.frombuffer(buffer.read(24), dtype=np.float64)
            # ispin = np.frombuffer(buffer.read(4), dtype=np.int32)[0]
            # gamma_only = np.frombuffer(buffer.read(8), dtype=np.int64)[0]
            # scalef = np.frombuffer(buffer.read(8), dtype=np.float64)[0]
            # buffer.read(8)  # newline
            # ---- uncomment above to read them but it may be wrong,
            # ---- currently they are skipped
            buffer.read(56)
            # ------------------------------------------------------------
            ngw, igwx_, npol, nbnd_ = \
                np.frombuffer(buffer.read(16), dtype=np.int32)
            # assert npol == 1, "npol != 1, not implemented"
            buffer.read(8)  # newline
            # ------------------------------------------------------------
            # b = np.frombuffer(buffer.read(72), dtype=np.float64).reshape(3, 3)
            # buffer.read(8)  # newline
            # ---- uncomment above to read bvecs
            buffer.read(80)
            # ------------------------------------------------------------
            # read gvecs
            gksize = 12 * igwx_  # 3 * int32 = 12 bytes
            gk = np.frombuffer(buffer.read(gksize), dtype=np.int32)
            gk = gk.reshape(gk.size//3, 3)
            buffer.read(8)  # newline
            # ------------------------------------------------------------
            # read wfc
            wfc = []
            wtmpsize = dsize * igwx_  # complex128 = 16 bytes
            for i in range(nbnd_):
                wtmp = np.frombuffer(buffer.read(wtmpsize), dtype=dtype)  # noqa
                wfc.append(wtmp)
                buffer.read(8)  # newline

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
    # check if calculation is gamma only
    isreal = _check_gamma_only(path)
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
            gk_k, evc_k = read_wfc_file(filename, isreal=isreal)
            evc[ik].append(evc_k)
            if ispin == 1:
                gk.append(gk_k)
    return gk, evc


def _check_gamma_only(path: str) -> bool:
    xmlfile = os.path.join(path, 'data-file-schema.xml')
    root = etree.parse(xmlfile).getroot()
    return root.find('input').find('basis').find('gamma_only').text == 'true'


# def determine_nspin(self):
#
#
# def read_eigenvalues_66(self, dtype=np.float64):
#     root = etree.parse(os.path.join(self.xmlfile)).getroot()
#     spin_child = root.find('spin')
#
#     band_child = root.find('band_structure')
#     nbnd = band_child.find('nbnd').text
#
#
# def xml_root(self):
#     return = etree.parse(self.xmlfile).getroot()
#
#
#     with qesave.xml_root() as root:
#         ....
#
