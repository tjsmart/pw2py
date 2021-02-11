import numpy as np
import scipy.fftpack

from ._writers import write_xsf


def reshape_wfc3D(wfc, gvec, scaling=2, shape=None):
    '''
    reshape wfc to 3D array defined over grid of gvec

    input
    ---
        wfc (np.ndarray)
            - wavefunction of shape = (len(gvec),)
        gvec (np.ndarray)
            - g-vectors of shape = (len(gvec), 3)
        scaling (int or float) (optional)
            - scaling factor to extend the size of the grid (padded with zeros), default=2
        shape (tuple of int) (optional)
            - shape of wfc3D, overrides scaling

    returns
    ---
        wfc3D (np.ndarray)
            - if scaling is used then shape = (np.amax(gvec, axis=0) * 2 + 1) * scaling
            - if shape is provided then shape = shape
    '''
    if shape is None:
        shape = (np.amax(gvec, axis=0) * 2 + 1) * scaling
        shape = tuple([int(s) for s in shape])
    else:
        shape = tuple(shape)
    # initialzie wfc3D with zeros (note indices not covered will default to zero if they are outside of grid gvec)
    wfc3D = np.zeros(shape, dtype=wfc.dtype)
    for ig in range(wfc.size):
        # not sure why, but minus sign is needed here based on testing
        wfc3D[tuple(-gvec[ig])] = wfc[ig]

    return wfc3D


def calc_wfc3D_squared_real(wfc, gvec, scaling=2, shape=None, lsign=False):
    '''
    preform fft of wfc(G) to real space and return 3D wfc(R)

    # TODO not inplace option ...
    '''
    # reshape wfc to be 3-dim array defined over gvec
    wfc3D = reshape_wfc3D(wfc, gvec, scaling=scaling, shape=shape)
    # fourier transform wfc to real space
    wfc3D = scipy.fftpack.fftn(wfc3D, overwrite_x=True) / (wfc3D.size)**0.5
    # calculate |wfc3D|^2 = conj(wfc3D) * wfc3D
    if lsign:
        wfc3D = np.real(np.sign(wfc3D)) * np.abs(wfc3D)**2
    else:
        wfc3D = np.abs(wfc3D)**2

    return wfc3D


def plot_wfc_xsf(filename, geo, wfc, gvec, scaling=2, shape=None, lsign=False):
    '''
    preform fft of wfc(G) to real space and plot 3D wfc(R) in xsf file
    '''
    # calculate fft
    wfc3D = calc_wfc3D_squared_real(
        wfc, gvec, scaling=scaling, shape=shape, lsign=lsign)
    # write xsf file
    write_xsf(filename, geo, grid=wfc3D)


def plot_wfc_averaged(filename, wfc, gvec, scaling=2, shape=None, lsign=False, free_axis=2):
    '''
    preform fft of wfc to real space and average wfc to a 1D function

    if free_axis = 0:
        wfc1D(x), is averaged over y, z
    elif free_axis = 1:
        wfc1D(y), is averaged over x, z
    elif free_axis = 2:
        wfc1D(z), is averaged over x, y
    '''
    # check free_axis
    assert free_axis in [0, 1, 2], "free_axis must be 0, 1, or 2"
    # calculate fft
    wfc3D = calc_wfc3D_squared_real(
        wfc, gvec, scaling=scaling, shape=shape, lsign=lsign)
    # axis to average over
    avg_axis = tuple([i for i in range(3) if i != free_axis])
    # calculate average
    wfc3D = np.average(wfc3D, axis=avg_axis)
    # dump to file
    np.savetxt(filename, wfc3D)


def plot_rho_xsf(filename, geo, rhog, gvec, scaling=1, shape=None):
    '''
    preform fft of rho(G) to real space and plot 3D rho(R) in xsf file
    '''
    # reshape wfc to be 3-dim array defined over gvec
    rho = reshape_wfc3D(rhog, gvec, scaling=scaling, shape=shape)
    # fourier transform wfc to real space
    rho = np.real(scipy.fftpack.fftn(rho, overwrite_x=True))
    # write xsf file
    write_xsf(filename, geo, grid=rho)
