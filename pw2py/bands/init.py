import numpy as np
import warnings


def __init__(self, shift: float = None, **kw):
    self._kpt = kw['kpt']
    self._eig = kw['eig']
    self._weight = kw['weight']
    if 'fermi' in kw:
        self._fermi = kw['fermi']
    if 'occ' in kw:
        self._occ = kw['occ']
    if 'vbm' in kw:
        self._vbm = kw['vbm']
    if 'cbm' in kw:
        self._cbm = kw['cbm']
        self._gap = self._cbm - self._vbm
    # -------- more setup ----------
    # this class assumes, occ, fermi, or vbm must be given
    ass = any([hasattr(self, p) for p in ['_occ', '_fermi', '_vbm']])
    assert ass, "No vbm, occupation or fermi level found in filename"
    # check if metallic, otherwise calc band edges
    self._check_is_metallic()
    # check if occupations exist
    if not self.is_metallic and not hasattr(self, '_occ'):
        _calc_occ(self)
    # calc bandedges, if not provided and system is not metallic
    if not self.is_metallic and not hasattr(self, '_vbm'):
        self._calc_bandedge()
    # optional shift of bands (useful for vacuum or core alignment)
    if shift is not None:
        self.shift_bands(shift)


def _check_is_metallic(self):
    if hasattr(self, '_vbm'):
        self._is_metallic = False
    elif hasattr(self, '_occ'):
        self._is_metallic = np.any(((self.occ != 0) & (self.occ != 1)))
    elif hasattr(self, '_fermi'):
        if (self.nk == 1):
            self._is_metallic = False
        else:
            self._is_metallic = False
            for ispin in range(self.nspin):
                if len(self.fermi) == 2:
                    is_occupied = (self.eig[ispin] < self.fermi[ispin])
                else:
                    is_occupied = (self.eig[ispin] < self.fermi)
                for ik in range(1, self.nk):
                    self._is_metallic = np.any(np.logical_xor(
                        is_occupied[0], is_occupied[ik]))
                    if self._is_metallic:
                        break


def _calc_occ(self):
    assert hasattr(self, '_fermi'), "Cannot compute occ without fermi"
    self._occ = np.ones(self.eig.shape)
    for ispin in range(self.nspin):
        if len(self.fermi) == 2:
            is_occupied = (self.eig[ispin] < self.fermi[ispin])
        else:
            is_occupied = (self.eig[ispin] < self.fermi)
        self._occ[ispin, ~is_occupied] = 0


def _calc_bandedge(self):
    if hasattr(self, '_occ'):
        is_vb = np.isclose(1, self.occ, atol=1e-3)
    # elif hasattr(self, '_fermi'):
    #     is_vb = self.eig < self.fermi
    self._vbm = self.eig[is_vb].max()
    if not np.all(is_vb):
        self._cbm = self.eig[~is_vb].min()
        self._gap = self.cbm - self.vbm
