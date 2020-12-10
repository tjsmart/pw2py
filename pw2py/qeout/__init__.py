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

        Methods:
            from_file       : create qeout object from_file -- intended usage
            calcEigs        : calculate eigenvalues vbm/cbm/gap (only works for gamma w/ smearing TODO)
            final_energy    : return final instance of energy (can specifiy convergence level)
            # TODO repr/str
    '''

    from .init import __init__

    from .io import from_file

    from .methods import calcEigs, final_energy, read_bands
