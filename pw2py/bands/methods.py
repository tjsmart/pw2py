def shift_bands(self, shift: float):
    '''
    constant shift of bands by float
    updates eig, fermi, vbm, and cbm (if applicapable)
    '''
    self._eig += shift
    if hasattr(self, 'fermi'):
        self._fermi += shift
    if not self.is_metallic:
        self._vbm += shift
        self._cbm += shift
