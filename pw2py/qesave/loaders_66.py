import numpy as np
import mmap


def read_wfc_file(filename: str):
    '''
    read wfc filename (dat format, not hdf5)

    WARNING!!! This should only be used for reading kp=(0, 0, 0), others are not tested

    WARNING!!! gamma_only is not implemented (should read real wfc instead of complex)

    WARNING!!! npol != 1 is not implemented (e.g. noncollinear case)
    '''
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
            wtmpsize = 16 * igwx_  # complex128 = 16 bytes
            for i in range(nbnd_):
                wtmp = np.frombuffer(buffer.read(wtmpsize), dtype=np.complex128)  # noqa
                wfc.append(wtmp)
                buffer.read(8)  # newline
            wfc = np.array(wfc, dtype=np.complex128)

    return gk, wfc
