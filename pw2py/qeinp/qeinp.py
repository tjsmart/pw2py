
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

        # TODO this will be fixed when par_units/pos_units move to be an attribute
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
            if np.array_equal(if_pos, np.array([1, 1, 1], dtype=int)):
                out += "\n"
            else:
                out += "    {}  {}  {}\n".format(if_pos[0], if_pos[1], if_pos[2])
        out += "K_POINTS {}\n".format(self.qedict['K_POINTS'])
        if self.qedict['K_POINTS'] == "automatic":
            out += "    {}  {}  {}  {}  {}  {}\n".format(
                self.kpt[0, 0], self.kpt[0, 1], self.kpt[0, 2], self.kpt[1, 0], self.kpt[1, 1], self.kpt[1, 2]
            )

        return out

    @classmethod
    def from_file(cls, filename):
        '''
        create qeinp object from file
        '''
        # intialize qedict object
        qedict = {'nml': f90nml.read(filename)}

        # store lines of file to iterator
        with open(filename) as f:
            lines = iter(f.readlines())

        # loop through lines for ['CELL_PARAMETERS', 'ATOMIC_SPECIES', 'ATOMIC_POSITIONS', 'K_POINTS']
        for line in lines:
            line = line.split('!')[0]   # trim away comments
            if 'CELL_PARAMETERS' in line:
                qedict['CELL_PARAMETERS'] = nonalpha.sub('', line.split('CELL_PARAMETERS')[1]).lower()
                # read cell parameters
                par = np.array([np.fromstring(next(lines).split('!')[0], sep=' ') for _ in range(3)], dtype=np.float64)

            elif 'ATOMIC_SPECIES' in line:
                qedict['ATOMIC_SPECIES'] = [next(lines).split('!')[0].split()
                                            for _ in range(int(qedict['nml']['system']['ntyp']))]

            elif 'ATOMIC_POSITIONS' in line:
                qedict['ATOMIC_POSITIONS'] = nonalpha.sub('', line.split('ATOMIC_POSITIONS')[1]).lower()
                # read ion, pos, if_pos
                ion, pos, if_pos = [], [], []
                for _ in range(int(qedict['nml']['system']['nat'])):
                    spl = next(lines).split('!')[0].split()
                    ion.append(spl[0])
                    pos.append(np.array(spl[1:4], dtype=np.float64))
                    try:
                        # explicit so indexError will be thrown
                        if_pos.append(np.array([spl[4], spl[5], spl[6]], dtype=int))
                    except IndexError:
                        if_pos.append(np.array([1, 1, 1], dtype=int))

                pos = np.array(pos)
                if_pos = np.array(if_pos)

            elif 'K_POINTS' in line:
                qedict['K_POINTS'] = nonalpha.sub('', line.split('K_POINTS')[1]).lower()
                # read k-points
                if qedict['K_POINTS'] == 'automatic':
                    kpt = np.fromstring(next(lines).split('!')[0], sep=' ', dtype=int).reshape((2, 3))
                elif qedict['K_POINTS'] == 'gamma':
                    kpt = None
                else:
                    # TODO add support for crystal_b etc.
                    raise ValueError("K_POINTS option '{}' not supported, please use 'automatic' or 'gamma'")

        if int(qedict['nml']['system']['ibrav']) != 0:
            par = _ibrav_to_par(qedict['nml']['system'])
            qedict['CELL_PARAMETERS'] = "angstrom"

        nat = qedict['nml']['system']['nat']

        return cls(qedict=qedict, par=par, ion=ion, pos=pos, if_pos=if_pos, kpt=kpt, nat=nat)

    def write_file(self, filename):
        '''
        create qeinput object from file
        '''
        with open(filename, 'w') as f:
            f.write(str(self))

        return None
