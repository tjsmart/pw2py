import os.path

from ..functions import final_energy
from ..bands import bands
from ..atomgeo import atomgeo


@classmethod
def from_file(cls, filename: str):
    '''
    create qeout object from file
    '''
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File not found: {filename}")
    # dictionary to collect results
    coll = {}
    # read bands from file
    coll['bands'] = bands.from_file(filename)
    # read geometry from file
    coll['geo'] = atomgeo.from_file(filename)
    # final energy
    coll['final_energy'] = final_energy(filename)
    # TODO read convergences and more (see below)
    return cls(**coll)

# def read_conv(filename: str) -> dict:
#     '''
#     Read various aspects of convergence from file
#     '''
#     conv = {}
#     count = {}
#
#     with open(filename, 'r') as f:
#         for line in f:
#             if "lattice parameter (alat)" in line:
#                 alat = np.float64(line.split()[-2]) * bohr_to_angstrom
#             elif "number of atoms/cell" in line:
#                 nat = int(line.split()[-1])
#             elif "number of electrons" in line:
#                 nelec = int(float(line.split()[4]))    # noqa: F841
#             elif "number of Kohn-Sham states" in line:
#                 nbnd = int(line.split()[-1])
#             elif "number of atomic types" in line:
#                 ntyp = int(line.split()[-1])
#             elif "EXX-fraction" in line:
#                 is_exx = True   # noqa: F841
#             elif "Starting magnetic structure" in line:
#                 is_mag = True
#             elif "crystal axes:" in line:
#                 par = np.array([next(lines)[0].split()[3:6]
#                                 for _ in range(3)], dtype=np.float64) * alat
#             elif "total cpu time spent up to now is" in line:
#                 conv['time'].append(np.float64(line.split()[-2]))
#             elif "total energy              =" in line:
#                 if "!!" in line:
#                     conv['!!'].append(np.float64(line.split()[-2]))
#                     conv['nsteps_exx'].append(
#                         sum(conv['nsteps']) - sum(conv['nsteps_exx']))
#                     conv['time_exx'].append(conv['time'][-1])
#                 elif "!" in line:
#                     conv['!'].append(np.float64(line.split()[-2]))
#                     next(lines)
#                     conv['dE'].append(np.float64(next(lines).split()[-2]))
#                     if is_mag:
#                         # skip throught till total magnetization
#                         while "total magnetization" not in line:
#                             line = next(lines)
#                         conv['tot_mag'].append(np.float64(line.split()[-3]))
#                         conv['abs_mag'].append(np.float64(next(lines).split()[-3]))
#
#                 else:
#                     conv['E'].append(np.float64(line.split()[-2]))
#
#             elif "has" in line:
#                 _nsteps = line.split()[-2]
#                 try:
#                     conv['nsteps'].append(int(_nsteps))
#                 except ValueError:
#                     continue
#             elif "Forces acting on atoms" in line:
#                 max_forc = -1.0
#                 next(lines)
#                 for _ in range(nat):
#                     this_forc = np.linalg.norm(np.fromstring(
#                         next(lines).split("force =")[1], sep=' ', dtype=np.float64))
#                     if this_forc > max_forc:
#                         max_forc = this_forc
#                 conv['max_forc'].append(max_forc)
#             elif "Total force" in line:
#                 conv['tot_forc'].append(np.float64(line.split()[3]))
#             elif "ATOMIC_POSITIONS" in line:
#                 if is_first:
#                     pos_units = line.split('ATOMIC_POSITIONS')[1]
#                     pos_units = ''.join(filter(str.isalpha, pos_units)).lower()
#                 pos = []
#                 for _ in range(nat):
#                     nl = next(lines).split()
#                     if is_first:
#                         ion.append(nl[0])
#                         try:
#                             if_pos.append(
#                                 np.array([nl[4], nl[5], nl[6]], dtype=int))
#                         except IndexError:
#                             if_pos.append(np.ones(3, dtype=int))
#                     pos.append(np.array(nl[1:4], dtype=np.float64))
#                 if is_first:
#                     is_first = False
#                 list_pos.append(np.array(pos))
