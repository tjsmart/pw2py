import numpy as np
import mmap
from savefolder import Savefolder
from typing import List, Union
import os.path


def wfc_filename_66(ik: int, ispin: Union[int, None]) -> str:
    '''
    Returns the name of the wfc file
    '''
    if ispin is None:
        sspin = ''
    elif ispin == 0:
        sspin = 'up'
    elif ispin == 1:
        sspin = 'dw'
    return f'wfc{sspin}{ik+1}.dat'


def read_kvec_66(savefolder: Savefolder, ik: int) -> List[float]:
    '''
    return the ik'th k-vector
    '''
    childs = savefolder.root.find('output').find(
        'band_structure').findall('ks_energies')
    return np.fromstring(childs[ik].find('k_point').text, sep=' ', dtype=np.float64)


def read_gamma_only_66(savefolder: Savefolder) -> bool:
    '''
    checks xml file for gamma_only flag
    '''
    return savefolder.root.find('input').find('basis').find('gamma_only').text == 'true'


def read_wfc_66(savefolder: Savefolder,
                ik: int,
                ib: int,
                ispin: Union[int, None] = None,
                dscale: Union[float, None] = None) -> dict:
    '''
    Read wavefunction from  qe-6.6 (dat format, not hdf5)

    WARNING!!! This should only be used for reading kp=(0, 0, 0), others are not tested

    WARNING!!! gamma_only is not implemented (should read real wfc instead of complex)

    WARNING!!! npol != 1 is not implemented (e.g. noncollinear case)
    '''
    kvec = read_kvec_66(savefolder, ik)
    isreal = read_gamma_only_66(savefolder)
    if isreal:
        raise NotImplementedError("Sorry! Gamma only is not yet implemented")

    filename = os.path.join(
        savefolder.savefolder, wfc_filename_66(ik, ispin=ispin))

    with open(filename, 'rb') as f:
        with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as buffer:
            # ------------------------------------------------------------
            seekinfo = 56
            buffer.seek(seekinfo)
            ngw, igwx_, npol, nbnd_ = np.frombuffer(
                buffer.read(16), dtype=np.int32)
            # ------------------------------------------------------------
            seekgk = seekinfo + 104
            buffer.seek(seekgk)
            gksize = 12 * igwx_  # 3*int32 = 12 bytes
            gk = np.frombuffer(buffer.read(gksize), dtype=np.int32)
            gk = gk.reshape(gk.size//3, 3)
            # ------------------------------------------------------------
            # TODO -- implement reading real wfc
            seek = seekgk + gksize + 8
            wfcsize = 16 * igwx_  # complex128 = 16 bytes
            seek += (wfcsize + 8) * ib
            buffer.seek(seek)
            evc = np.frombuffer(buffer.read(wfcsize), dtype=np.complex128)

    return dict(npw=igwx_,
                gk=gk,
                evc=evc,
                kvec=kvec,
                isreal=isreal,
                dscale=dscale
                )
