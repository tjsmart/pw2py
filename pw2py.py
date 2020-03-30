#!/usr/bin/env python3

import copy
import f90nml
import numpy as np
import pandas as pd
import re

import tjs_resource as tjs

# bohr to angstrom
bohr_to_angstrom = 0.529177210903
# filter for stripping non-alphabetic characters from string
nonalpha = re.compile('[^a-zA-Z]')


def _resolve_continuation_lines(filename):
    '''
    return lines of file without continuation lines (used for jdftx)
    '''
    lines = []
    with open(filename) as f:
        for line in f:
            line = line.rstrip('\n')
            while line.endswith('\\'):
                line = line[:-1] + next(f).rstrip('\n')
            lines.append(line)
    return lines


def _read_atomic_positions(lines, nat, no_if_pos=True):
    '''
    iterate through lines of atomic positons returning ion, pos, if_pos (QE format)
    '''
    pos, ion, if_pos = [], [], []
    for _ in range(nat):
        nl = next(lines).split('!')[0].split()
        ion.append(nl[0])
        pos.append(np.array(nl[1:4], dtype=np.float64))
        if not no_if_pos:
            try:
                if_pos.append(np.array([nl[4], nl[5], nl[6]], dtype=int)) # explicit so indexError will be thrown
            except IndexError:
                if_pos.append(np.array([1,1,1], dtype=int))
    if no_if_pos:
        return lines, np.array(pos), ion
    else:
        return lines, np.array(pos), ion, np.array(if_pos)


def _ibrav_to_par(system_nml, units="angstrom"):
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
        tjs.die("Invalid value of alat_units: '{}".format(alat_units))
    # define unit dictionary for unit conversion
    conversion = {'alat' : alat, 'angstrom' : 1, 'bohr' : bohr_to_angstrom}
    # check input
    if in_units == out_units:
        return par
    elif not in_units in conversion:
        tjs.die("Invalid value of in_units: '{}'".format(in_units))
    elif not out_units in conversion:
        tjs.die("Invalid value of out_units: '{}'".format(out_units))
    elif (in_units == "alat" or out_units == "alat") and alat == None:
        tjs.die("Requested alat conversion but value for alat was not provided.")
    # calculate prefactor
    return par / np.float64( conversion[out_units.lower()] ) * np.float64( conversion[in_units.lower()] )


def convert_pos(pos, in_units, out_units="angstrom", alat=None, alat_units="angstrom", par=None, par_units=None):
    '''
    Convert atomic positions from 'in_units' to 'out_units'
        ['alat', 'bohr', 'angstrom', 'crystal']
    '''
    if in_units == out_units:
        return pos
    elif in_units != 'crystal' and out_units != 'crystal':
        return convert_par(pos, in_units=in_units, out_units=out_units, alat=alat, alat_units=alat_units)
    elif par is None or par_units is None:
        tjs.die("Requested crystal conversion but cell parameters (par) were not provided.")

    # cases with crystal conversion
    prefactor = 1.0
    if (par_units == 'angstrom' and out_units == 'bohr') or (par_units == 'bohr' and in_units == 'angstrom'):
        prefactor /= bohr_to_angstrom
    elif (par_units == 'bohr' and out_units == 'angstrom') or (par_units == 'angstrom' and in_units == 'bohr'):
        prefactor *= bohr_to_angstrom
    
    if in_units == 'crystal':
        return pos.dot(par) * prefactor
    elif out_units == 'crystal':
        return pos.dot(np.linalg.inv(par)) * prefactor


def _readEigs(lines, nbnd):
    '''
    (private) parse lines for eigenvalues
    '''
    next(lines)
    eigs = []
    for _ in range(nbnd // 8 + 1):
        eigs.append( np.fromstring(next(lines), sep=' ', dtype=np.float64) )

    return np.hstack(eigs), lines


class atomgeo:
    '''
    class for atomic geometry

    Suggested to make object from file:
    -----
    geo = atomgeo.from_file(filename, ftype="auto")
    
    Attributes:
    -----

    par (np.array, float, shape = (3,3))
        - cell parameters
    
    par_units (str)
        - units of par, acceptable values include "angstrom", "bohr", "alat"

    nat (int)
        - number of atoms
    
    ion (list, str)
        - list of atom/ion names (QE format)
    
    pos (np.array, float, shape = (nat, 3))
        - atomic positions
    
    pos_units (str)
        - units of pos, acceptable values include "angstrom", "bohr", "crystal", "alat"
    '''


    def __init__(self, par=None, ion=None, pos=None, nat=None, par_units=None, pos_units=None):
    # def __init__(self, par=None, ion=None, pos=None, if_pos=None, nat=None):
        self.par        = np.array(par, dtype=np.float) if par is not None else None
        self.ion        = list(ion) if ion is not None else None
        self.pos        = np.array(pos, dtype=np.float) if pos is not None else None
        self.par_units  = par_units
        self.pos_units  = pos_units
        self.nat        = int(nat)


    def __str__(self):
        '''
        convert atomgeo to str (QE format)
        '''
        # TODO need to be able to handle cases with ibrav != 0
        ntyp = len(set(self.ion))
        out = "&control\n/\n&system\n    ibrav = 0\n    ntyp = {}\n    nat = {}\n/\n&electrons\n/\n"\
            .format(ntyp, str(self.nat))
        out += "CELL_PARAMETERS {}\n".format(self.par_units)
        for par in self.par:
            out += "    {:16.9f}  {:16.9f}  {:16.9f}\n".format(par[0], par[1], par[2])
        out += "ATOMIC_POSITIONS {}\n".format(self.pos_units)
        for ion, pos in zip(self.ion, self.pos):
            out += "    {:5s}  {:16.9f}  {:16.9f}  {:16.9f}\n".format(ion, pos[0], pos[1], pos[2])
        return out

    
    # def replace_atom(self, a, b):
    #     '''
    #     replace atom a with atom b
    #     '''
    #     pass
        
    
    # def add_atom(self, a):
    #     '''
    #     add atom a
    #     '''
    #     replace(self, a, None)
    #     pass
        
    
    # def remove_atom(self, a):
    #     '''
    #     remove atom a
    #     '''
    #     replace(self, None, a)
    #     pass


    def change_units_par(self, units="angstrom", inplace=True):
        '''
        Convert cell parameters to 'units'

        units (str)
            - acceptable values: ['alat', 'angstrom', 'bohr']
        
        updates values of: par, par_units, (qedict['CELL_PARAMETERS'])
        '''
        out = self if inplace else copy.deepcopy(self)

        out.par = convert_par(out.par, out.par_units, out_units=units)
        out.par_units = units
        # TODO bottom part should only belong to qeinp
        try:
            out.qedict["CELL_PARAMETERS"] = units
        except:
            None


    def change_units_pos(self, units="angstrom", inplace=True):
        '''
        Convert atomic positions to 'units'

        units (str)
            - acceptable values: ['alat', 'bohr', 'angstrom', 'crystal']
        
        updates values of: pos, pos_units, (qedict['ATOMIC_POSITIONS'])
        '''
        out = self if inplace else copy.deepcopy(self)

        out.pos = convert_pos(out.pos, out.pos_units, out_units=units, par=out.par, par_units=out.par_units)
        out.pos_units = units
        # TODO bottom part should only belong to qeinp
        try:
            out.qedict["ATOMIC_POSITIONS"] = units
        except:
            None
        return out


    def from_file(filename, ftype='auto'):
        '''
        generate atomgeo object from file

        input
        ----
        filename (str)
            - name/path to file to be read
        
        ftype (str)
            - specify filetype otherwise ftype='auto' will use filename to detect file type
            - acceptable values: ['auto', 'vasp', 'qeinp', 'qeout', 'jdftx', 'xsf', 'xyz']
        
        output
        ----
        atomgeo (atomgeo object)
        '''
        if ftype.lower() not in ['auto', 'vasp', 'qeinp', 'qeout', 'jdftx', 'xsf', 'xyz']:
            tjs.warn("ftype value not recognized '{}'".format(ftype))

        # store lines of file to iterator
        with open(filename) as f:
            lines = iter(f.readlines())

        if any([filename.lower().endswith(vasp.lower()) for vasp in ["OUTCAR", "CONTCAR", "vasp"]]) \
            or ftype.lower() == "vasp":
            # then file is vasp
            next(lines) # skip first line
            # read cell par
            alat = np.float64(next(lines))
            par = np.array([ next(lines).split()[0:3] for _ in range(3) ], dtype=np.float64) * alat
            par_units = "angstrom"
            # read ions
            _ion = next(lines).split()
            _count = np.fromstring(next(lines), sep=' ', dtype=int)
            ion = []
            for i, c in zip(_ion, _count):
                ion += [i] * c
            nat = len(ion)
            # read pos type
            pos_units = "crystal" if next(lines).strip().lower() == "direct" else "angstrom"
            pos = np.array([ next(lines).split()[0:3] for _ in range(nat) ], dtype=np.float64)

        elif filename.lower().endswith("xyz") or ftype.lower() == "xyz":
            # then file is xyz format
            tjs.die("xyz format not implemented")

        elif filename.lower().endswith("xsf") or ftype.lower() == "xsf":
            # then file is xsf format
            tjs.die("xsf format not implemented")

        elif filename.lower().endswith("jdftx") or filename.lower().endswith("pos") \
            or ftype.lower() == "jdftx":
            # reread lines where all continuation lines have been resolved (i.e. lines that end in '\')
            lines = iter(_resolve_continuation_lines(filename))
            ion = []
            pos = []
            for line in lines:
                line = line.split('#')[0].strip()
                if line.startswith('lattice '):
                    par = np.fromstring(line.strip('lattice'), dtype=float, sep=' ').reshape((3,3)).T
                elif line.startswith('ion '):
                    ion.append(line.split()[1])
                    pos.append(np.array(line.split()[2:5], dtype=float))

            pos = np.array(pos)
            nat = len(ion)
            par_units = "bohr"
            tjs.warn("pos_units not set, defaulting to crystal")
            pos_units = "crystal"

        else:
            # then file is QE
            nml = f90nml.read(filename)
            if len(nml) == 0 or ftype.lower() == "qeout":
                # then file is output
                for line in lines:
                    if "lattice parameter (alat)" in line:
                        alat = np.float64(line.split()[-2]) * bohr_to_angstrom
                    elif "number of atoms/cell" in line:
                        nat = int(line.split()[-1])
                    elif "crystal axes:" in line:
                        # cell parameters in angstrom
                        par_units = "angstrom"
                        par = np.array([ next(lines).split()[3:6] for _ in range(3) ], dtype=np.float64) * alat
                    elif "ATOMIC_POSITIONS" in line:
                        pos_units = nonalpha.sub('', line.split('ATOMIC_POSITIONS')[1]).lower()
                        lines, pos, ion = _read_atomic_positions(lines, nat, no_if_pos=True)
            else:
                # ftype.lower() == "qeinp", no need to actually check
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
                        lines, pos, ion = _read_atomic_positions(lines, nat, no_if_pos=True)

                if int(nml['system']['ibrav']) != 0:
                    # if lattice is specified by ibrav then build par from ibrav
                    par = _ibrav_to_par(nml['system'])
                    par_units = "angstrom"
            
            # # convert atomic positions to angstrom
            # convert_pos(pos, pos_units, alat=alat, par=par)
        return atomgeo(par=par, ion=ion, pos=pos, nat=nat, pos_units=pos_units, par_units=par_units)
    

    def sort_ions(self, inplace=False):
        '''
        sort ions into categories and reorder corresponding positions
        '''
        out = self if inplace else copy.deepcopy(self)

        # save ion and pos to dataframe
        columns = ['ion'] + ['pos{}'.format(i) for i in range(3)]
        df = pd.DataFrame(columns=columns)
        df.ion = self.ion
        for i in range(3):
            df['pos{}'.format(i)] = self.pos[:,i]

        # sort dataframe by ion column
        df.sort_values(by=['ion'], inplace=True)

        # update ion and pos
        out.ion = list(df['ion'])
        out.pos = np.array(df.loc[:,'pos0':'pos2'])
        return out


    def write_file(self, filename):
        '''
        write atomgeo to file (QE format)
        '''
        with open(filename, 'w') as f:
            f.write(str(self))

        return None


    def write_vasp(self, filename, alat=1.0, inplace=False):
        '''
        write atomgeo to poscar/vasp file type
        '''
        out = self.sort_ions(inplace=inplace)

        with open(filename, 'w') as f:
            # header
            f.write("Generated from PW2PY python module\n")
            # alat
            f.write("{}\n".format(float(alat)))
            # par
            for par in out.par:
                f.write( "    {:16.9f}  {:16.9f}  {:16.9f}\n".format(par[0], par[1], par[2]) )
            # ions and count (convert to dict then write)
            ionAndCount = {}
            for i in out.ion:
                if i not in ionAndCount:
                    ionAndCount[i] = 1
                else:
                    ionAndCount[i] += 1
            f.write( "     {}\n".format("     ".join(ionAndCount.keys())) )             # ions
            f.write( "     {}\n".format("     ".join(map(str, ionAndCount.values()))) ) # count
            # pos type
            if out.pos_units == "angstrom":
                f.write("{}\n".format("cartesian"))
            elif out.pos_units == "crystal":
                f.write("{}\n".format("direct"))
            # pos
            for pos in out.pos:
                f.write( "    {:16.9f}  {:16.9f}  {:16.9f}\n".format(pos[0], pos[1], pos[2]) )
        
        return None


    def write_jdftx(self, filename, inplace=False, if_pos=0, write_latt=True):
        '''
        write qegeo to jdftx
        '''
        out = self if inplace else copy.deepcopy(self)

        # change par and pos to bohr
        out.change_units_par(units="bohr", inplace=True)

        with open(filename, 'w') as f:
            f.write("# Generated from PW2PY python module\n")
            if write_latt:
                # write lattice
                f.write("lattice \\\n")
                for i, par in enumerate(out.par.T):
                    if i == 0 or i == 1:
                        f.write( "    {:16.9f}  {:16.9f}  {:16.9f} \\\n".format(par[0], par[1], par[2]) )
                    elif i == 2:
                        f.write( "    {:16.9f}  {:16.9f}  {:16.9f}\n".format(par[0], par[1], par[2]) )
            # write ion and pos
            for ion, pos in zip(out.ion, out.pos):
                f.write("ion  {:5s}  {:16.9f}  {:16.9f}  {:16.9f} {}\n".format(ion, pos[0], pos[1], pos[2], if_pos))


    def write_xyz(self, filename, inplace=False):
        '''
        write qegeo to xyz
        '''
        out = self.change_units_pos(units="angstrom", inplace=inplace)

        with open(filename, 'w') as f:
            # nat
            f.write("{}\n".format(out.nat))
            # description
            f.write("Generated from PW2PY python module\n")
            for ion, pos in zip(out.ion, out.pos):
                f.write("    {:5s}  {:16.9f}  {:16.9f}  {:16.9f}\n".format(ion, pos[0], pos[1], pos[2]))



    # def write_xsf():
    #     '''
    #     write qegeo to xsf
    #     '''
    #     pass


class qegeo(atomgeo):
    '''
    same as atomgeo with if_pos
    '''


    def __init__(self, par=None, ion=None, pos=None, if_pos=None, nat=None, par_units=None, pos_units=None):
        super().__init__(par=par, ion=ion, pos=pos, nat=nat, par_units=par_units, pos_units=pos_units)
        self.if_pos     = np.array(if_pos, dtype=np.int) if if_pos is not None else None


    def __str__(self, par_units="angstrom", pos_units="crystal"):
        '''
        convert qegeo to str (QE format)
        '''
        # TODO need to be able to handle cases with ibrav != 0
        ntyp = len(set(self.ion))
        out = "&control\n/\n&system\n    ibrav = 0\n    ntyp = {}\n    nat = {}\n/\n&electrons\n/\n"\
            .format(ntyp, str(self.nat))
        out += "CELL_PARAMETERS {}\n".format(self.par_units)
        for par in self.par:
            out += "    {:16.9f}  {:16.9f}  {:16.9f}\n".format(par[0], par[1], par[2])
        out += "ATOMIC_POSITIONS {}\n".format(self.pos_units)
        for ion, pos, if_pos in zip(self.ion, self.pos, self.if_pos):
            out += "    {:5s}  {:16.9f}  {:16.9f}  {:16.9f}\n".format(ion, pos[0], pos[1], pos[2])
            if np.array_equal(if_pos, np.array([1,1,1], dtype=int)):
                out += "\n"
            else:
                try:
                    out += "    {}  {}  {}\n".format(if_pos[0], if_pos[1], if_pos[2])
                except IndexError:
                    None

        return out


    def load_if_pos(self, filename, if_pos_all=[1,1,1]):
        '''
        load if_pos data from filename, default to if_pos_all in many cases
        '''

        if any([filename.lower().endswith(vasp.lower()) for vasp in ["OUTCAR", "CONTCAR", "vasp"]]):
            self.if_pos = np.array([ if_pos_all for _ in range(nat) ], dtype=int)

        elif filename.lower().endswith("xyz"):
            self.if_pos = np.array([ if_pos_all for _ in range(nat) ], dtype=int)

        elif filename.lower().endswith("xsf"):
            self.if_pos = np.array([ if_pos_all for _ in range(nat) ], dtype=int)

        elif filename.lower().endswith("jdftx") or filename.lower().endswith("pos"):
            #TODO read if_pos from jdftx file
            tjs.die("jdftx format not implemented")

        else:
            # store lines of file to iterator
            with open(filename) as f:
                lines = iter(f.readlines())
            # then file is QE
            nml = f90nml.read(filename)
            if len(nml) == 0:
                # then file is output
                for line in lines:
                    if "ATOMIC_POSITIONS" in line:
                        _, _, _, self.if_pos = _read_atomic_positions(lines, nat, no_if_pos=False)
            else:
                for line in lines:
                    line = line.split('!')[0]   # trim away comments
                    if 'ATOMIC_POSITIONS' in line:
                        _, _, _, self.if_pos = _read_atomic_positions(lines, nat, no_if_pos=False)

        return self.if_pos
    

    # TODO implement sort_ions with if_pos
    # def sort_ions(self, inplace=False):
    #     '''
    #     sort ions into categories and reorder corresponding positions
    #     '''
    #     super().sort_ions(self)
    #     self.if_pos = np.array(df.loc[:,'if_pos0':'if_pos2'], dtype=int)



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

    def __init__(self, qedict, par, ion, pos, if_pos, kpt, nat=None, par_units=None, pos_units=None):
        '''
        initialize qeinp object
        '''
        # if any(None in locals().values()):
        #     tjs.die("Cannot initialize qeinp with 'None'\n{}".format(locals))

        if par_units is None:
            # TODO what to do here when ibrav != 0!!! -- should 
            try:
                par_units = qedict['CELL_PARAMETERS']
            except:
                par_units = None
        else:
            par_units = par_units
        
        if pos_units is None:
            # TODO what to do here when ibrav != 0!!! -- should 
            try:
                pos_units = qedict['ATOMIC_POSITIONS']
            except:
                pos_units = None
        else:
            pos_units = pos_units
        

        super().__init__(par=par, ion=ion, pos=pos, if_pos=if_pos, nat=nat, par_units=par_units, pos_units=pos_units)
        self.qedict = qedict
        self.kpt = kpt
    
    
    # TODO
    # def convert_ibrav(self, newBrav):
    #     '''
    #     convert ibrav of self to newBrav (intended for ibrav != 0 to 0 or back)
    #     '''
    #     pass
    

    def __str__(self):
        '''
        convert qeinp to string
        '''
        out = str(self.qedict['nml']) + "\n"
        out += "ATOMIC_SPECIES\n"
        for spec in self.qedict['ATOMIC_SPECIES']:
            out += "    {}  {}  {}\n".format(spec[0], spec[1], spec[2])
        # TO-DO need to be able to handle cases with ibrav != 0
        # if "CELL_PARAMETERS" in self.qedict.keys():
        if self.qedict['nml']['system']['ibrav'] == 0:
            out += "CELL_PARAMETERS {}\n".format(self.qedict['CELL_PARAMETERS'])
            for par in self.par:
                out += "    {:16.9f}  {:16.9f}  {:16.9f}\n".format(par[0], par[1], par[2])
        out += "ATOMIC_POSITIONS {}\n".format(self.qedict['ATOMIC_POSITIONS'])
        for ion, pos, if_pos in zip(self.ion, self.pos, self.if_pos):
            out += "    {}    {:16.9f}  {:16.9f}  {:16.9f}".format(ion, pos[0], pos[1], pos[2])
            if np.array_equal(if_pos, np.array([1,1,1], dtype=int)):
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
                elif qedict['K_POINTS'] == 'gamma':
                    kpt = None
                else:
                    tjs.die("K_POINTS option '{}' not supported, please use 'automatic' or 'gamma'")
        
        if int(qedict['nml']['system']['ibrav']) != 0:
            par = _ibrav_to_par(qedict['nml']['system'])
            qedict['CELL_PARAMETERS'] = "angstrom"

        nat = qedict['nml']['system']['nat']

        return qeinp(qedict=qedict, par=par, ion=ion, pos=pos, if_pos=if_pos, kpt=kpt, nat=nat)


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

    def __init__(self, nat, ntyp, conv, par, ion, list_pos, if_pos, pos_units, list_eigs):
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
        self.list_eigs = list_eigs


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
                    'mag_tot', 'mag_abs', 'time', 'time_exx', 'tot_mag', 'abs_mag', 'Fermi']:
            conv[key] = []
        list_pos, if_pos, ion, list_eigs = [], [], [], []

        is_first = True
        for line in lines:
            if "lattice parameter (alat)" in line:
                alat = np.float64(line.split()[-2]) * bohr_to_angstrom
            elif "number of atoms/cell" in line:
                nat = int(line.split()[-1])
            elif "number of electrons" in line:
                nelec = int(float(line.split()[-1]))
            elif "number of Kohn-Sham states" in line:
                nbnd = int(line.split()[-1])
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
                    [ next(lines) for _ in range(10) ]
                    conv['tot_mag'].append(np.float64(next(lines).split()[-3]))
                    conv['abs_mag'].append(np.float64(next(lines).split()[-3]))

                else:
                    conv['E'].append(np.float64(line.split()[-2]))
            elif "has" in line:
                conv['nsteps'].append(int(line.split()[-2]))
            elif "SPIN UP" in line:
                [ next(lines) for _ in range(2) ]
                if  "k = 0.0000 0.0000 0.0000" in line:
                    #TODO need to other kpoints
                    eigs, lines = _readEigs(lines, nbnd)
                    list_eigs.append(eigs)
                    

            elif "SPIN DOWN" in line:
                [ next(lines) for _ in range(2) ]
                if  "k = 0.0000 0.0000 0.0000" in line:
                    eigs, lines = _readEigs(lines, nbnd)
                    list_eigs.append(eigs)
    
            elif "k = 0.0000 0.0000 0.0000" in line:
                eigs, lines = _readEigs(lines, nbnd)
                list_eigs.append(eigs)

            elif "Fermi" in line:
                conv['Fermi'].append(np.float64(line.split()[-2]))
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

        return qeout(nat, ntyp, conv, par, ion, list_pos, if_pos, pos_units, list_eigs)
    

    def calcEigs(self):
        '''
        calculate VBM, CBM, and gap (only works with smearing)
        '''
        try:
            fermi = self.conv['Fermi'][-1]
        except:
            # TODO -- implement no smearing case
            tjs.warn("In calcEigs: Cannot calculate eigs for systems w/o smearing")

        eigs = self.list_eigs[-1]
        vbm = max(np.where(eigs < fermi, eigs, -1E10))
        cbm = min(np.where(eigs > fermi, eigs, 1E10))

        return vbm, cbm, cbm - vbm
