
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

    @staticmethod
    def final_energy(filename, conv_level='automatic', units='ev'):
        '''
        grab the final energy from filename and return convergence level

        filename (str):
            - quantum espresso output file

        conv_level (str):
            - 'total energy'    : return last instance of total energy (could be from inner SCF loop)
            - '!'               : return last instance of total energy with SCF loop convergence
            - '!!'              : return last instance of total energy from EXX loop convergence
            - 'Final'           : return final energy (from relax/vc-relax)
            - 'automatic'       : return highest convergence level achieved

        returns
        ---
            tuple(energy, conv_level)
                - energy(float) = energy in Ry
                - conv_level(str) = level of convergence achieved
        '''
        # build dictionary of energy instances
        energy_instances = {}
        instances = ['Final', '!!', '!', 'total energy']
        conv_level = conv_level.lower()
        # check given conv_level is a valid option
        if conv_level not in instances + ['automatic']:
            return None, "invalid conv_level: '{}'".format(conv_level)
        # load file to memory, storing in lines
        with open(filename) as f:
            lines = f.readlines()
        # iterate backwards through file
        for line in reversed(lines):
            # skip lines without 'energy' in them
            if 'energy' not in line:
                continue
            # otherwise iterate through possible instances and extract energy when found
            for instance in instances:
                if instance in line and instance not in energy_instances:
                    try:
                        energy_instances[instance] = float(line.split('=')[1].split('R')[0])
                    except IndexError:
                        # false line, i.e. no total energy was given
                        pass
            # if all instances have been found then break loop
            if instances == list(energy_instances.keys()):
                break

        # now return the appropriate energy, based on user conv_level
        if conv_level == 'automatic':
            # for automatic iterate through instances and return first available instance
            for instance in instances:
                if instance in energy_instances:
                    if units == 'eV':
                        return energy_instances[instance] * ry_to_ev, instance
                    else:
                        return energy_instances[instance] * ry_to_ev, instance

        elif conv_level not in energy_instances:
            # user specified conv_level was not found in the file
            return None, "conv_level '{}' not achieved".format(conv_level)

        else:
            # return user specified conv_level
            if units == 'eV':
                return energy_instances[conv_level] * ry_to_ev, conv_level
            else:
                return energy_instances[conv_level], conv_level

        # this should not happen
        return None, 'Unexpected error in final_energy'

    @classmethod
    def from_file(cls, filename):
        '''
        create qeout object from file
        '''
        # store lines of file to iterator
        with open(filename) as f:
            lines = iter(f.readlines())

        is_mag = False
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
            elif "Starting magnetic structure" in line:
                is_mag = True
            elif "crystal axes:" in line:
                par = np.array([next(lines)[0].split()[3:6] for _ in range(3)], dtype=np.float64) * alat
            elif "total cpu time spent up to now is" in line:
                conv['time'].append(np.float64(line.split()[-2]))
            elif "total energy              =" in line:
                if "!!" in line:
                    conv['!!'].append(np.float64(line.split()[-2]))
                    conv['nsteps_exx'].append(sum(conv['nsteps']) - sum(conv['nsteps_exx']))
                    conv['time_exx'].append(conv['time'][-1])
                elif "!" in line:
                    conv['!'].append(np.float64(line.split()[-2]))
                    next(lines)
                    conv['dE'].append(np.float64(next(lines).split()[-2]))
                    if is_mag:
                        # skip throught till total magnetization
                        while "total magnetization" not in line:
                            line = next(lines)
                        conv['tot_mag'].append(np.float64(line.split()[-3]))
                        conv['abs_mag'].append(np.float64(next(lines).split()[-3]))

                else:
                    conv['E'].append(np.float64(line.split()[-2]))

            elif "has" in line:
                conv['nsteps'].append(int(line.split()[-2]))
            elif "SPIN UP" in line:
                [next(lines) for _ in range(2)]
                if "k = 0.0000 0.0000 0.0000" in line:
                    # TODO need to other kpoints
                    eigs, lines = _readEigs(lines, nbnd)
                    list_eigs.append(eigs)

            elif "SPIN DOWN" in line:
                [next(lines) for _ in range(2)]
                if "k = 0.0000 0.0000 0.0000" in line:
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
                    this_forc = np.linalg.norm(np.fromstring(
                        next(lines).split("force =")[1], sep=' ', dtype=np.float64))
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
                            if_pos.append(np.array([nl[4], nl[5], nl[6]], dtype=int))
                        except IndexError:
                            if_pos.append(np.ones(3, dtype=int))
                    pos.append(np.array(nl[1:4], dtype=np.float64))
                if is_first:
                    is_first = False
                list_pos.append(np.array(pos))

        return cls(nat, ntyp, conv, par, ion, list_pos, if_pos, pos_units, list_eigs)

    def calcEigs(self):
        '''
        calculate VBM, CBM, and gap (only works with smearing)
        '''
        try:
            fermi = self.conv['Fermi'][-1]
        except KeyError:
            # TODO -- implement no smearing case
            tjs.warn("In calcEigs: Cannot calculate eigs for systems w/o smearing")

        eigs = self.list_eigs[-1]
        vbm = max(np.where(eigs < fermi, eigs, -1E10))
        cbm = min(np.where(eigs > fermi, eigs, 1E10))

        return vbm, cbm, cbm - vbm
