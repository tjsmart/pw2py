import numpy as np


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
            is_occupied = (self.eig < self.fermi)
            for ispin in range(self.nspin):
                for ik in range(1, self.nk):
                    self._is_metallic = np.any(np.logical_xor(
                        is_occupied[ispin, 0], is_occupied[ispin, ik]))
                    if self._is_metallic:
                        break


def _calc_bandedge(self):
    if hasattr(self, '_fermi'):
        is_vb = self.eig < self.fermi
    elif hasattr(self, '_occ'):
        is_vb = np.isclose(1, self.occ, atol=1e-3)
    self._vbm = self.eig[is_vb].max()
    if not np.all(is_vb):
        self._cbm = self.eig[~is_vb].min()
        self._gap = self.cbm - self.vbm
