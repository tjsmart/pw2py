import numpy as np
import scipy.fft
import os.path

from .savefolder import Savefolder
# from .loaders import read_wfc_61, read_wfc_66
from .loaders import read_wfc_66
from ..functions._writers import write_xsf  # noqa


'''
TODO
[X] should have bands object which defines eigs and occ
[X] should define cell object which defines par and reciprocal lattice
        or add reciprocal lattice to qeinp/atomgeo
        a subclass for par may be more ideal
[X] other attributes:
        evc : np.ndarray with dimensions of [ispin, ik, iband, G]
        gk : np.ndarray with dimensions of [ispin, ik, iband, G]
        fft : List[int, int, int]
        isreal : bool (e.g. is real or not)
[X] methods:
        fft R -> G

[X] issue:
        what if we only want to read a portion of the wavefunction e.g. a single band?
        what to do if it is not spin polarized?

Maybe we should define a class for a single wavefunction, then a class for a full 'qesave'

[X] several attributes, such as scaling, should have a setter (if it is changed ffts and fftp should be updated)
[X] other attributes should be read only e.g. evc, gk, etc...


'''


class Wavefunction:
    def __init__(self, **kw):
        '''
        initialize Wavefunction from keyword arguments

        input
        ---
            kw: dict, with keys:
                npw (opt.)      : number of plane waves
                gk              : int array, asserted to be of shape (npw, 3), if npw is given
                evc             : float or complex array, should be same length as gk's first dimension (e.g. npw)
                kvec (opt.)     : float array, length 3 (assumed to be the gamma if not given!!)
                isreal (opt.)   : must be provided as true if the evc is real, otherwise it will be converted to complex
                dscale (opt.)   : scaling for dense grid (integer, default is 2)

        returns
        ---
            Wavefunction object

        '''
        # ------------------------------------------------------------
        # set npw (number of plane waves), if given
        if 'npw' in kw:
            assert isinstance(kw['npw'], (int, np.integer)), \
                "npw should be an integer"
            self.npw = kw['npw']
        # ------------------------------------------------------------
        # set gk vectors (g-vectors at the k-point)
        self.gk = np.array(kw['gk'], dtype=int)
        if 'npw' in kw:
            if self.gk.shape != (self.npw, 3):
                self.gk = self.gk.reshape(self.npw, 3)
        else:
            assert (len(self.gk.shape) == 2) and (self.gk.shape[1] == 3), \
                f"wrong shape for gk, should be (npw, 3): {self.gk.shape}"
            self.npw = self.gk.shape[0]
        # ------------------------------------------------------------
        # set kvec
        if 'kvec' in kw:
            self.kvec = np.array(kw['kvec'], dtype=np.float64)
            assert self.kvec.shape == (3,), \
                f"wrong shape for kvec: {self.kvec.shape}"
        # ------------------------------------------------------------
        # set isreal
        if 'isreal' in kw:
            self.isreal = bool(kw['isreal'])
        else:
            self.isreal = False
        # ------------------------------------------------------------
        # set eigenvector (the wavefunction in G space)
        evc_dtype = np.float64 if self.isreal else np.complex128
        self.evc = np.array(kw['evc'], dtype=evc_dtype)
        assert self.evc.shape == (self.npw,), \
            f"wrong shape for evc: {self.evc.shape}"
        # ------------------------------------------------------------
        # set dscale (scaling factor for dense grid)
        if ('dscale' not in kw) or (kw['dscale'] is None):
            self.dscale = 2
        else:
            self.dscale = float(kw['dscale'])
        # ------------------------------------------------------------
        # set ffts and fftd, smooth and dense grid sizes
        self.ffts = 1 + 2 * np.amax(self.gk, axis=0)
        self.fftd = self.ffts * self.dscale

    @classmethod
    def from_folder(cls,
                    folder,
                    ik,
                    ib,
                    ispin=None,
                    dscale=None,
                    version=None):
        '''
        read wavefunction from a folder

        if calculation is not spin polarized then ispin will be ignored
        but otherwise it should provided as a keyword argument

        input
        ---
            folder          : path to *.save folder in QE, or can be folder containing *.save if it is unique
            ik              : int array, asserted to be of shape (npw, 3), if npw is given
            ib              : float or complex array, should be same length as gk's first dimension (e.g. npw)
            ispin (opt.)    : spin of the wfc (either 1, 2, or None)
            dscale (opt.)   : scaling for dense grid (integer, default is 2)
            version (opt.)  : version of QE folder, will be raises error if different than implemented version

        returns
        ---
            Wavefunction object
        '''
        # ------------------------------------------------------------
        # create savefolder object
        if not os.path.exists(folder):
            raise FileNotFoundError(f"Folder does not exist: {folder}")
        if folder.endswith('.save'):
            prefix = folder[:len('.save')]
            path = os.path.relpath(folder, '..')
        else:
            prefix = None
            path = folder
        savefolder = Savefolder(prefix=prefix, path=path)
        # ------------------------------------------------------------
        # check version of savefolder
        if version is None:
            version = savefolder.version
        assert version in ['6.1', '6.6'], (
            f"Detected version: {version}."
            " However, only 6.1 and 6.6 are implemented and verified."
            " Pass a version number of 6.1 or 6.6 to force reading"
        )
        # ------------------------------------------------------------
        # call the loader
        if version == '6.1':
            # TODO
            raise NotImplementedError("6.1 not yet implemented")
            # return cls(**read_wfc_61(savefolder,
            #                   ik,
            #                   ib,
            #                   ispin=None,
            #                   dscale=None))
        elif version == '6.6':
            return cls(**read_wfc_66(savefolder,
                                     ik,
                                     ib,
                                     ispin=ispin,
                                     dscale=None))

    def reshape_evc3D(self):
        '''
        reshape evc from 1D to 3D
        '''
        evc3D = np.zeros(self.fftd, dtype=self.evc.dtype)
        evc3D[-self.gk[:, 0], -self.gk[:, 1], -self.gk[:, 2]] = self.evc
        return evc3D

    def fft_evcR(self) -> None:
        '''
        perform fourier transform of wave form G -> R

        sets attribute self.evcR which has dimensions of fftd
        '''
        if self.isreal:
            raise NotImplementedError("FFT of real -> complex not implemented")
        self.evcR = self.reshape_evc3D()
        self.evcR = scipy.fft.fftn(self.evcR, overwrite_x=True)
        self.evcR /= (self.evcR.size)**0.5

    def set_evc(self, evc):
        '''
        change the value of evc
        '''
        raise NotImplementedError("TODO")

    def plot_wfc_xsf(self, filename, geo, lsign=False):
        '''
        plot 3D evcR in xsf file

        if evcR is not defined, calls self.fft_evcR() first
        '''
        if not hasattr(self, 'evcR'):
            self.fft_evcR()
        if lsign:
            grid = np.real(np.sign(self.evcR)) * np.abs(self.evcR)**2
        else:
            grid = np.abs(self.evcR)**2
        write_xsf(filename, geo, grid=grid)
