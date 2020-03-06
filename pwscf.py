#!/usr/bin/env python3

import numpy as np
import re
import f90nml

import tjs_resource as tjs

# bohr to angstrom
bohr_to_angstrom = 0.529177210903
# filter for stripping non-alphabetic characters from string
nonalpha = re.compile('[^a-zA-Z]')


def read_atomic_positions(lines, nat, only_pos=False):
    '''
    iterate through lines of atomic positons returning ion, pos, if_pos
    '''
    pos, ion, if_pos = [], [], []
    for _ in range(nat):
        nl = next(lines).split('!')[0].split()
        pos.append(np.array(nl[1:4], dtype=np.float64))
        if not only_pos:
            ion.append(nl[0])
            try:
                if_pos.append(np.array([nl[4], nl[5], nl[6]], dtype=int)) # explicit so indexError will be thrown
            except IndexError:
                if_pos.append(np.array([1,1,1], dtype=int))
    if only_pos:
        return lines, np.array(pos)
    else:
        return lines, np.array(pos), ion, np.array(if_pos)


def ibrav_to_par(system_nml, units="angstrom"):
    '''
    calculate cell parameters from bravais lattice using qe system namelist as input

    note the default cell parameters are in angstrom
    '''
    # check value of ibrav
    if int(system_nml['ibrav']) == 0:
        return None
    
    elif int(system_nml['ibrav']) == 1:
        # v1 = a(1,0,0),  v2 = a(0,1,0),  v3 = a(0,0,1)
        par = np.array( [[1,0,0],[0,1,0],[0,0,1]] , dtype=np.float64)
    elif int(system_nml['ibrav']) == 2:
        # v1 = (a/2)(-1,0,1),  v2 = (a/2)(0,1,1), v3 = (a/2)(-1,1,0)
        par = np.array( [[-1,0,1],[0,1,1],[-1,1,0]] , dtype=np.float64) * 0.5
    elif int(system_nml['ibrav']) == 3:
        # v1 = (a/2)(1,1,1),  v2 = (a/2)(-1,1,1),  v3 = (a/2)(-1,-1,1)
        par = np.array( [[1,1,1],[-1,1,1],[-1,-1,1]] , dtype=np.float64) * 0.5
    elif int(system_nml['ibrav']) == 4:
        # v1 = a(1,0,0),  v2 = a(-1/2,sqrt(3)/2,0),  v3 = a(0,0,c/a)
        par = np.array( [[1,0,0],[-1/2,np.sqrt(3)/2,0],[0,0,1]] , dtype=np.float64)
    elif int(system_nml['ibrav']) == 6:
        # v1 = a(1,0,0),  v2 = a(0,1,0),  v3 = a(0,0,c/a)
        par = np.array( [[1,0,0],[0,1,0],[0,0,1]] , dtype=np.float64)
    elif int(system_nml['ibrav']) == 8:
        # v1 = (a,0,0),  v2 = (0,b,0), v3 = (0,0,c)
        par = np.array( [[1,0,0],[0,1,0],[0,0,1]] , dtype=np.float64)
    else:
        tjs.die("Unsupported value for ibrav '{}'".format(system_nml['ibrav']))

    # rescale lattice parameters
    if 'a' in system_nml.keys():
        par *= np.float64(system_nml['a'])
        if int(system_nml['ibrav']) in [8]:
            par[1,1] = np.float64(system_nml['b'])
        if int(system_nml['ibrav']) in [4, 6, 8]:
            par[2,2] = np.float64(system_nml['c'])
    elif 'celldm' in system_nml.keys():
        par *= np.float64(system_nml['celldm'][0]) * bohr_to_angstrom
        if int(system_nml['ibrav']) in [8]:
            par[1,1] *= np.float64(system_nml['celldm'][1])
        if int(system_nml['ibrav']) in [4, 6, 8]:
            par[2,2] *= np.float64(system_nml['celldm'][2])
    
    # check units
    if units == "bohr":
        par /= bohr_to_angstrom
    
    return par


def convert_par(par, in_units, out_units="angstrom", alat=None, alat_units="angstrom"):
    '''
    Convert cell parameters from 'in_units' to 'out_units'
        ['alat', 'angstrom', 'bohr']
    '''
    # convert alat
    if alat_units == "bohr":
        alat *= alat * bohr_to_angstrom
    elif alat_units != "angstrom":
        die("Invalid value of alat_units: '{}".format(alat_units))
    # define unit dictionary for unit conversion
    conversion = {'alat' : alat, 'angstrom' : 1, 'bohr' : bohr_to_angstrom}
    # check input
    if in_units == out_units:
        return par
    elif not in_units in conversion:
        die("Invalid value of in_units: '{}'".format(in_units))
    elif not out_units in conversion:
        die("Invalid value of out_units: '{}'".format(out_units))
    elif (in_units == "alat" or out_units == "alat") and alat == None:
        die("Requested alat conversion but value for alat was not provided.")
    # calculate prefactor
    return par / np.float64( conversion[out_units.lower()] ) / np.float64( conversion[in_units.lower()] )


def convert_pos(pos, in_units, out_units="angstrom", alat=None, alat_units="angstrom", par=None):
    '''
    Convert atomic positions from 'in_units' to 'out_units'
        ['alat', 'bohr', 'angstrom', 'crystal']
    '''
    # calculate prefactor
    if in_units == out_units:
        return pos
    elif in_units != 'crystal' and out_units != 'crystal':
        return convert_par(pos, in_units=in_units, out_units=out_units, alat=alat, alat_units=alat_units)
    elif par is None:
        die("Requested crystal conversion but cell parameters (par) were not provided.")
    elif in_units == 'crystal':
        return pos.dot(par)
    elif out_units == 'crystal':
        return pos.dot(np.linalg.inv(par))


class qegeo:
    '''
    class for atomic geometry included cell parameters, atomic positions, species names, and if_pos
    '''

    def calcNat(self):
        lengths = [ len(t) for t in [self.ion, self.pos, self.if_pos] if t is not None]
        if len(lengths) == 0:
            self.nat = None
        else:
            if all( [ lengths[0] == l for l in lengths ] ):
                self.nat = int(lengths[0])
            else:
                self.nat = None


    # def __init__(self, par=None, ion=None, pos=None, if_pos=None, par_units=None, pos_units=None):
    def __init__(self, par=None, ion=None, pos=None, if_pos=None, nat=None):
        self.par        = np.array(par, dtype=np.float) if par is not None else None
        self.ion        = list(ion) if ion is not None else None
        self.pos        = np.array(pos, dtype=np.float) if pos is not None else None
        self.if_pos     = np.array(if_pos, dtype=np.int) if if_pos is not None else None
        # self.par        = par
        # self.ion        = ion
        # self.pos        = pos
        # self.if_pos     = if_pos
        if nat is None:
            self.calcNat()

    
    def from_file(filename):
        '''
        generate qegeo from file (automatically detects if file is input or output)
        '''
        # store lines of file to iterator
        with open(filename) as f:
            lines = iter(f.readlines())

        nml = f90nml.read(filename)
        if len(nml) == 0:
            # then file is output
            for line in lines:
                if "lattice parameter (alat)" in line:
                    alat = np.float64(line.split()[-2]) * bohr_to_angstrom
                elif "number of atoms/cell" in line:
                    nat = int(line.split()[-1])
                elif "crystal axes:" in line:
                    # cell parameters in angstrom
                    par = np.array([ next(lines).split()[3:6] for _ in range(3) ], dtype=np.float64) * alat
                elif "ATOMIC_POSITIONS" in line:
                    pos_units = nonalpha.sub('', line.split('ATOMIC_POSITIONS')[1]).lower()
                    lines, pos, ion, if_pos = read_atomic_positions(lines, nat)
        else:
            # then file is input
            alat = None
            for line in lines:
                if 'a' in nml['system']:
                    alat = nml['system']['a']
                elif 'celldm' in nml['system']:
                    alat = nml['system']['celldm'][0] * bohr_to_angstrom
                nat = int(nml['system']['nat'])
                line = line.split('!')[0]   # trim away comments
                if 'CELL_PARAMETERS' in line:
                    par_units = nonalpha.sub('', line.split('CELL_PARAMETERS')[1]).lower()
                    par = np.array( [np.fromstring(next(lines).split('!')[0], sep=' ') for _ in range(3)] , dtype=np.float64 )
                    par = convert_par(par, in_units=par_units, alat=alat)
                elif 'ATOMIC_POSITIONS' in line:
                    pos_units = nonalpha.sub('', line.split('ATOMIC_POSITIONS')[1]).lower()
                    lines, pos, ion, if_pos = read_atomic_positions(lines, nat)

            if int(nml['system']['ibrav']) != 0:
                par = ibrav_to_par(nml)
        
        # convert atomic positions to angstrom
        convert_pos(pos, pos_units, alat=alat, par=par)

        return qegeo(par=par, ion=ion, pos=pos, if_pos=if_pos, nat=nat)



class qeinp(qegeo):
    '''
    class for quantum espresso input file
        qedict = dictionary of namelists (nml), ATOMIC_SPECIES, CELL_PARAMETERS, ATOMIC_POSITIONS, and K_POINTS
        par = CELL_PARAMETERS data
        ion = ATOMIC_POSITIONS ions data
        pos = ATOMIC_POSITIONS position data
        if_pos = ATOMIC_POSITIONS if_pos data
        kpt = K_POINTS data
    '''

    def __init__(self, qedict, par, ion, pos, if_pos, kpt):
        '''
        initialize qeinp object
        '''
        # if any(None in locals().values()):
        #     tjs.die("Cannot initialize qeinp with 'None'\n{}".format(locals))

        super().__init__(par=par, ion=ion, pos=pos, if_pos=if_pos)
        self.qedict = qedict
        self.kpt = kpt
    

    def __str__(self):
        '''
        convert qeinp to string
        '''
        out = str(self.qedict['nml']) + "\n"
        out += "ATOMIC_SPECIES\n"
        for spec in self.qedict['ATOMIC_SPECIES']:
            out += "    {}  {}  {}\n".format(spec[0], spec[1], spec[2])
        out += "CELL_PARAMETERS {}\n".format(self.qedict['CELL_PARAMETERS'])
        for par in self.par:
            out += "    {:16.9f}  {:16.9f}  {:16.9f}\n".format(par[0], par[1], par[2])
        out += "ATOMIC_POSITIONS {}\n".format(self.qedict['ATOMIC_POSITIONS'])
        for ion, pos, if_pos in zip(self.ion, self.pos, self.if_pos):
            out += "    {}    {:16.9f}  {:16.9f}  {:16.9f}".format(ion, pos[0], pos[1], pos[2])
            if all(if_pos == np.array([1,1,1], dtype=int)):
                out += "\n"
            else:
                out += "    {}  {}  {}\n".format(if_pos[0], if_pos[1], if_pos[2])
        out += "K_POINTS {}\n".format(self.qedict['K_POINTS'])
        if self.qedict['K_POINTS'] == "automatic":
            out += "    {}  {}  {}  {}  {}  {}\n".format(
                self.kpt[0,0], self.kpt[0,1], self.kpt[0,2], self.kpt[1,0], self.kpt[1,1], self.kpt[1,2]
            )

        return out


    def from_file(filename):
        '''
        create qeinp object from file
        '''
        # intialize qedict object
        qedict = {'nml' : f90nml.read(filename)}

        # store lines of file to iterator
        with open(filename) as f:
            lines = iter(f.readlines())

        # loop through lines for ['CELL_PARAMETERS', 'ATOMIC_SPECIES', 'ATOMIC_POSITIONS', 'K_POINTS']
        for line in lines:
            line = line.split('!')[0]   # trim away comments
            if 'CELL_PARAMETERS' in line:
                try:
                    qedict['CELL_PARAMETERS'] = nonalpha.sub('', line.split('CELL_PARAMETERS')[1]).lower()
                except:
                    tjs.die("Invalid CELL_PARAMETERS line: \n{}".format(line))
                # read cell parameters
                par = np.array( [np.fromstring(next(lines).split('!')[0], sep=' ') for _ in range(3)] , dtype=np.float64 )

            elif 'ATOMIC_SPECIES' in line:
                qedict['ATOMIC_SPECIES'] = [next(lines).split('!')[0].split() for _ in range(int(qedict['nml']['system']['ntyp']))]

            elif 'ATOMIC_POSITIONS' in line:
                try:
                    qedict['ATOMIC_POSITIONS'] = nonalpha.sub('', line.split('ATOMIC_POSITIONS')[1]).lower()
                except:
                    tjs.die("Invalid ATOMIC_POSITIONS line: \n{}".format(line))
                # read ion, pos, if_pos
                ion, pos, if_pos = [], [], []
                for _ in range(int(qedict['nml']['system']['nat'])):
                    spl = next(lines).split('!')[0].split()
                    ion.append(spl[0])
                    pos.append(np.array(spl[1:4], dtype=np.float64))
                    try:
                        if_pos.append(np.array([spl[4], spl[5], spl[6]], dtype=int)) # explicit so indexError will be thrown
                    except IndexError:
                        if_pos.append(np.array([1,1,1], dtype=int))
                        
                pos = np.array(pos)
                if_pos = np.array(if_pos)

            elif 'K_POINTS' in line:
                try:
                    qedict['K_POINTS'] = nonalpha.sub('', line.split('K_POINTS')[1]).lower()
                except:
                    tjs.die("Invalid K_POINTS line: \n{}".format(line))
                # read k-points
                if qedict['K_POINTS'] == 'automatic':
                    kpt = np.fromstring(next(lines).split('!')[0], sep=' ', dtype=int).reshape((2,3))
                elif qedict['K_POINTS'] not in [ 'automatic', 'gamma']:
                    tjs.die("K_POINTS option '{}' not supported, please use 'automatic' or 'gamma'")
        
        if int(qedict['nml']['system']['ibrav']) != 0:
            par = ibrav_to_par(qedict['nml'])

        return qeinp(qedict, par, ion, pos, if_pos, kpt)

    def write_file(self, filename):
        '''
        create qeinput object from file
        '''
        with open(filename, 'w') as f:
            f.write(str(self))

        return None


class qeout():
    '''
    class for quantum espresso output file
        Input objects:
            nat         = number of atoms
            ntyp        = number of atomic types
            # TODO calc        = calculation type [scf, nscf, relax, vc-relax]
            # TODO is_exx      = hybrid calculation flag
            # TODO kpt         = kpt data
            # TODO celldm      = lattice parameter

        Convergence (conv dictionary):
            E           = list of total energy
            nsteps      = number of steps for scf convergence
            nsteps_exx  = number of steps for hybrid convergence
            !           = list of total energy (after scf convergence, only if calc is hybrid and/or relax)
            tot_forc    = list of total force
            max_forc    = list of maximum force acting on an atom
            !!          = list of hybrid total energy (if calc is hybrid)
            tot_mag     = list of total magnetization
            abs_mag     = list of absolute magnetization
            time        = length of time for scf steps
            time_exx    = length of time for exx steps

        Geometry data:
            par         = CELL_PARAMETERS data (if calculation is vc-relax)
            ion         = ATOMIC_POSITIONS ions data (if calculation is relax/vc-relax)
            list_pos    = ATOMIC_POSITIONS position data (if calculation is relax/vc-relax)
            if_pos      = ATOMIC_POSITIONS if_pos data (if present)
            pos_units   = units of position data

        Other:
            # TODO eig         = eigenvalue data
            # TODO occ         = occupation data
    '''

    def __init__(self, nat, ntyp, conv, par, ion, list_pos, if_pos, pos_units):
        '''
        initialize qeout object
        '''
        # if any(None in locals().values()):
        #     tjs.die("Cannot initialize qeinp with 'None'\n{}".format(locals))
        self.nat = nat
        self.ntyp = ntyp
        self.conv = conv
        self.par = par
        self.ion = ion
        self.list_pos = list_pos
        self.if_pos = if_pos
        self.pos_units = pos_units


    def from_file(filename):
        '''
        create qeout object from file
        '''
        # store lines of file to iterator
        with open(filename) as f:
            lines = iter(f.readlines())

        is_exx = False
        conv = {}
        for key in ['E', '!', '!!', 'dE', 'nsteps', 'nsteps_exx', 'tot_forc', 'max_forc', 
                    'mag_tot', 'mag_abs', 'time', 'time_exx', 'tot_mag', 'abs_mag']:
            conv[key] = []
        list_pos, if_pos, ion = [], [], []

        is_first = True
        for line in lines:
            if "lattice parameter (alat)" in line:
                alat = np.float64(line.split()[-2]) * bohr_to_angstrom
            elif "number of atoms/cell" in line:
                nat = int(line.split()[-1])
            elif "number of atomic types" in line:
                ntyp = int(line.split()[-1])
            elif "EXX-fraction" in line:
                is_exx = True
            elif "crystal axes:" in line:
                par = np.array([ next(lines)[0].split()[3:6] for _ in range(3) ], dtype=np.float64) * alat
            elif "total cpu time spent up to now is" in line:
                conv['time'].append(np.float64(line.split()[-2]))
            elif "total energy              =" in line:
                if "!!" in line:
                    conv['!!'].append(np.float64(line.split()[-2]))
                    conv['nsteps_exx'].append( sum(conv['nsteps']) - sum(conv['nsteps_exx']) )
                    conv['time_exx'].append( conv['time'][-1] )
                elif "!" in line:
                    conv['!'].append(np.float64(line.split()[-2]))
                    next(lines)
                    conv['dE'].append(np.float64(next(lines).split()[-2]))
                    next(lines)
                    conv['tot_mag'].append(np.float64(next(lines).split()[-3]))
                    conv['abs_mag'].append(np.float64(next(lines).split()[-3]))

                else:
                    conv['E'].append(np.float64(line.split()[-2]))
            elif "has" in line:
                conv['nsteps'].append(int(line.split()[-2]))
            elif "Forces acting on atoms" in line:
                max_forc = -1.0
                next(lines)
                for _ in range(nat):
                    this_forc = np.linalg.norm(np.fromstring(next(lines).split("force =")[1], sep=' ', dtype=np.float64))
                    if this_forc > max_forc:
                        max_forc = this_forc
                conv['max_forc'].append(max_forc)
            elif "Total force" in line:
                conv['tot_forc'].append(np.float64(line.split()[3]))
            elif "ATOMIC_POSITIONS" in line:
                if is_first:
                    pos_units = nonalpha.sub('', line.split('ATOMIC_POSITIONS')[1]).lower()
                pos = []
                for _ in range(nat):
                    nl = next(lines).split()
                    if is_first:
                        ion.append(nl[0])
                        try:
                            if_pos.append(np.array([nl[4],nl[5],nl[6]], dtype=int))
                        except:
                            if_pos.append(np.array([1,1,1], dtype=int))
                    pos.append(np.array(nl[1:4], dtype=np.float64))
                if is_first:
                    is_first = False
                list_pos.append(np.array(pos))

        return qeout(nat, ntyp, conv, par, ion, list_pos, if_pos, pos_units)
