@property
def eig(self):
    ''' return self._eig '''
    return self._eig


@property
def kpt(self):
    ''' return self._kpt '''
    return self._kpt


@property
def weight(self):
    ''' return self._weight '''
    return self._weight


@property
def fermi(self):
    ''' return self._fermi or None '''
    return self._fermi if hasattr(self, '_fermi') else None


@property
def occ(self):
    ''' return self._occ or None '''
    return self._occ if hasattr(self, '_occ') else None


@property
def is_metallic(self):
    ''' return self._is_metallic '''
    return self._is_metallic


@property
def vbm(self):
    ''' return self._vbm '''
    return self._vbm


@property
def cbm(self):
    ''' return self._cbm '''
    return self._cbm


@property
def gap(self):
    ''' return self._gap '''
    return self._gap


@property
def nspin(self):
    ''' return self.eig.shape[0] '''
    return self.eig.shape[0]


@property
def nk(self):
    ''' return self.eig.shape[1] '''
    return self.eig.shape[1]


@property
def nbnd(self):
    ''' return self.eig.shape[2] '''
    return self.eig.shape[2]
