from .. import atomgeo


class qegeo(atomgeo):
    '''
    same as atomgeo with if_pos and ibrav capabilities
    '''

    from ._init import __init__
    from ._str import __str__
    from ._properties import A, B, C, celldm, ibrav, if_pos

    # TODO CONTINUE HERE!!!! Start by fixing str method and clearning up below mehtods:

    # def load_if_pos(self, filename, if_pos_all=[1, 1, 1]):
    #     '''
    #     load if_pos data from filename, default to if_pos_all in many cases
    #     '''

    #     # TODO this likely needs to be rewritten or checked

    #     if any([filename.lower().endswith(vasp.lower()) for vasp in ["OUTCAR", "CONTCAR", "vasp"]]):
    #         self.if_pos = np.array([if_pos_all for _ in range(self.nat)], dtype=int)

    #     elif filename.lower().endswith("xyz"):
    #         self.if_pos = np.array([if_pos_all for _ in range(self.nat)], dtype=int)

    #     elif filename.lower().endswith("xsf"):
    #         self.if_pos = np.array([if_pos_all for _ in range(self.nat)], dtype=int)

    #     elif filename.lower().endswith("jdftx") or filename.lower().endswith("pos"):
    #         # TODO read if_pos from jdftx file
    #         raise ValueError('jdftx format not implemented')

    #     else:
    #         # store lines of file to iterator
    #         with open(filename) as f:
    #             lines = iter(f.readlines())
    #         # then file is QE
    #         nml = f90nml.read(filename)
    #         if len(nml) == 0:
    #             # then file is output
    #             for line in lines:
    #                 if "ATOMIC_POSITIONS" in line:
    #                     _, _, _, self.if_pos = _read_atomic_positions(lines, nat, no_if_pos=False)
    #         else:
    #             for line in lines:
    #                 line = line.split('!')[0]   # trim away comments
    #                 if 'ATOMIC_POSITIONS' in line:
    #                     _, _, _, self.if_pos = _read_atomic_positions(lines, nat, no_if_pos=False)

    #     return self.if_pos

    # TODO implement sort_ions with if_pos
    # def sort_ions(self, inplace=False):
    #     '''
    #     sort ions into categories and reorder corresponding positions
    #     '''
    #     super().sort_ions(self)
    #     self.if_pos = np.array(df.loc[:,'if_pos0':'if_pos2'], dtype=int)
