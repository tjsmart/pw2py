#!/usr/bin/env python

import os
import numpy as np
from lxml import etree
import mmap
import glob


def determine_dtype(attrib):
    '''
    attrib : dictionary from xml
    return : dtype based on kind and type
    '''
    dtype = None
    # determine dtype
    if attrib['type'] == 'complex':
        if int(attrib['kind']) == 8:
            dtype = np.complex128
        elif int(attrib['kind']) == 4:
            dtype = np.complex64
    elif attrib['type'] == 'real':
        if int(attrib['kind']) == 8:
            dtype = np.float64
        elif int(attrib['kind']) == 4:
            dtype = np.float32
    elif attrib['type'] == 'integer':
        if int(attrib['kind']) == 4:
            dtype = np.int32
        elif int(attrib['kind']) == 2:
            dtype = np.float16
        elif int(attrib['kind']) == 1:
            dtype = np.float8
    # if dtype wasn't set raise error
    if dtype is None:
        raise ValueError('Unexpected kind: {}, for type {}'.format(attrib['kind'], attrib['type']))
    # otherwise return dtype
    return dtype


def parse_dtype_and_size(buffer: mmap.mmap, pattern: str):
    '''
    parse data type of evc and byte size
    '''
    # find and seek pattern
    ix = buffer.find(pattern.encode())
    if ix == -1:
        raise ValueError("'{}' does not exist in file?".format(pattern))
    buffer.seek(ix)

    # parse evc info line as xml
    info_line = buffer.readline().replace(b'>', b'/>')
    tree = etree.fromstring(info_line)
    dtype = determine_dtype(tree.attrib)

    # read 3 ints of size 4 bytes
    head = np.frombuffer(buffer.read(3 * 4), dtype=np.int32)
    size = head[1] - 4

    # return dtype and size
    return dtype, size


def read_info_fft_grid(buffer: mmap.mmap):
    '''
    parse fft_grid size from info line
    '''
    ix = buffer.find(b'    <INFO ')
    buffer.seek(ix)
    tree = etree.fromstring(buffer.readline())
    fft_grid = np.array(tree.attrib.values(), dtype=np.int32)

    return fft_grid


def read_data_xml(path):
    '''
    parse key values from data_xml file
    '''
    data_xml = os.path.join(path, 'data-file.xml')
    root = etree.parse(data_xml).getroot()

    '''
    read cell
    '''
    cell_child = root.find('CELL')

    # read cell parameters
    par = np.empty((3, 3), dtype=np.float64)
    for i in range(3):
        text = cell_child.find('DIRECT_LATTICE_VECTORS').find('a{}'.format(i + 1)).text
        par[i] = np.fromstring(text, sep=' ', dtype=np.float64)

    # read reciprocal lattice vectors
    rec = np.empty((3, 3), dtype=np.float64)
    for i in range(3):
        text = cell_child.find('RECIPROCAL_LATTICE_VECTORS').find('b{}'.format(i + 1)).text
        rec[i] = np.fromstring(text, sep=' ', dtype=np.float64)

    '''
    read fft
    '''
    # read fft grid dimensions
    values = root.find('PLANE_WAVES').find('FFT_GRID').attrib.values()
    fft_grid = np.array(values, dtype=np.int32)

    '''
    read nspin, nkpt, nbnd, nbgkvec
    '''
    bands_child = root.find('BAND_STRUCTURE_INFO')
    # read nspin
    nspin = int(bands_child.find('NUMBER_OF_SPIN_COMPONENTS').text)

    # read nkpt
    nkpt = int(bands_child.find('NUMBER_OF_K-POINTS').text)

    # read nbnd
    nbnd = int(bands_child.find('NUMBER_OF_BANDS').text)

    # read ngkvec
    ngkvec = []
    for i in range(nkpt):
        ngkvec.append(int(root.find('EIGENVECTORS').find('K-POINT.{}'.format(i+1)).find('NUMBER_OF_GK-VECTORS').text))

    # return parsed data
    return par, rec, fft_grid, nspin, nkpt, nbnd, ngkvec


def read_eigenvalues(path, dtype=np.float64):
    '''
    read eigenvalues from qe save folder
    '''
    eig = []
    for ikpt, kdir in enumerate(glob.glob(os.path.join(path, 'K*'))):
        eig.append([])
        for ispin in range(2):
            xml_file = os.path.join(kdir, 'eigenval{}.xml'.format(ispin + 1))
            root = etree.parse(xml_file).getroot()
            eig_child = root.find('EIGENVALUES')
            eig[ikpt].append(np.fromstring(eig_child.text, sep=' ', dtype=dtype))

    return np.array(eig)


def read_gvectors(path):
    with open(os.path.join(path, 'gvectors.dat')) as f:
        # buffer file
        with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as buffer:
            dtype, size = parse_dtype_and_size(buffer, '<g ')
            gvecs = np.frombuffer(buffer.read(size), dtype=dtype)
            gvecs = gvecs.reshape((int(gvecs.size / 3), 3))
            return gvecs


def read_charge_density(path, filename='charge-density.dat'):
    with open(os.path.join(path, filename)) as f:
        # buffer file
        with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as buffer:
            fft_grid = read_info_fft_grid(buffer)
            rho = []
            for iz in range(fft_grid[2]):
                dtype, size = parse_dtype_and_size(buffer, '<z.{} '.format(iz + 1))
                rho_part = np.frombuffer(buffer.read(size), dtype=dtype)
                rho.append(rho_part)
            rho = np.array(rho).reshape(fft_grid)
            return rho


def read_spin_polarization(path, filename='spin-polarization.dat'):
    return read_charge_density(path, filename=filename)


def read_wavefunction(path):
    '''
    read evc from qe save folder
    '''
    def parse_nbnd():
        '''
        parse the info line for the value of nbnd
            "<INFO ngw=. igwx=. gamma_only=. nbnd=. ik=. ..."
        '''
        ix = buffer.find(b'  <INFO ')
        buffer.seek(ix)
        info_tree = etree.fromstring(buffer.readline())
        nbnd = int(info_tree.attrib['nbnd'])
        return nbnd

    evc = []
    for ikpt, kdir in enumerate(glob.glob(os.path.join(path, 'K*'))):
        # loop through kpts in folders such as 'K00001'
        evc.append([])
        for ispin in range(2):
            # loop through spin, opening files evc1.dat and evc2.dat
            evc[ikpt].append([])
            evc_file = os.path.join(kdir, 'evc{}.dat'.format(ispin + 1))
            # open evc_file
            with open(evc_file, 'rb') as f:
                # buffer file since it is large
                with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as buffer:
                    nbnd = parse_nbnd()
                    for iband in range(nbnd):
                        # loop through bands
                        pattern = '<evc.{} '.format(iband + 1)
                        # parse dtype and size (in bytes)
                        dtype, size = parse_dtype_and_size(buffer, pattern)
                        # read evc data for this band
                        evc_band = np.frombuffer(buffer.read(size), dtype=dtype)
                        evc[ikpt][ispin].append(evc_band)

    return np.array(evc)


def read_gkvectors(path):
    '''
    read array of g-vectors for each k-point
    '''
    gkvectors = []
    for kdir in glob.glob(os.path.join(path, 'K*')):
        with open(os.path.join(kdir, 'gkvectors.dat'), 'rb') as f:
            # buffer file since it is large
            with mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ) as buffer:
                # parse dtype and size (in bytes)
                dtype, size = parse_dtype_and_size(buffer, '  <GRID ')
                # read g-vectors at this kpt
                gvecs_kpt = np.frombuffer(buffer.read(size), dtype=dtype)
                # reshape and append to full list of gkvectors
                gkvectors.append(gvecs_kpt.reshape((int(gvecs_kpt.size / 3), 3)))

    return np.array(gkvectors)


def try_reading_all():
    import time

    path = '/Users/tjsmart/Programs/Ping-Group/NonRad/Examples/local/lin-gs/ratio-0.0000/BN.save'

    absolute_begin = time.time()

    begin = time.time()
    print('reading par, rec, fft_grid, nspin, nkpt, nbnd, ngkvec')
    par, rec, fft_grid, nspin, nkpt, nbnd, ngkvec = read_data_xml(path)
    end = time.time()
    print('... done in {}'.format(end - begin))
    print(par)
    print(rec)
    print(fft_grid)
    print(nspin)
    print(nkpt)
    print(nbnd)
    print(ngkvec)

    begin = time.time()
    print('\nreading eig')
    eig = read_eigenvalues(path)
    end = time.time()
    print('... done in {}'.format(end - begin))
    print(eig.size)
    print(eig.shape)
    print(eig.dtype)
    print(eig[0, 0, 0])

    begin = time.time()
    print('\nreading gvecs')
    gvecs = read_gvectors(path)
    end = time.time()
    print('... done in {}'.format(end - begin))
    print(gvecs.size)
    print(gvecs.shape)
    print(gvecs.dtype)
    print(gvecs[0])

    begin = time.time()
    print('\nreading charge density')
    rho = read_charge_density(path)
    end = time.time()
    print('... done in {}'.format(end - begin))
    print(rho.size)
    print(rho.shape)
    print(rho.dtype)
    print(rho[0, 0, 0])

    begin = time.time()
    print('\nreading spin polarization')
    sigma = read_spin_polarization(path)
    end = time.time()
    print('... done in {}'.format(end - begin))
    print(sigma.size)
    print(sigma.shape)
    print(sigma.dtype)
    print(sigma[0, 0, 0])

    begin = time.time()
    print('\nreading wavefunction')
    evc = read_wavefunction(path)
    end = time.time()
    print('... done in {}'.format(end - begin))
    print(evc.ndim)
    print(evc.shape)
    print(evc.size)
    print(evc.dtype)
    print(evc[0, 0, 0, 0])

    begin = time.time()
    print('\nreading kpt gvectors')
    gkvectors = read_gkvectors(path)
    end = time.time()
    print('... done in {}'.format(end - begin))
    print(gkvectors.ndim)
    print(gkvectors.shape)
    print(gkvectors.size)
    print(gkvectors.dtype)
    print(gkvectors[0, 0])

    absolute_end = time.time()

    print('\nDone reading everything in {}'.format(absolute_end - absolute_begin))


if __name__ == "__main__":
    try_reading_all()
