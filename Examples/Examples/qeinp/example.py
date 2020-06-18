#!/usr/bin/env python3

import pw2py as pw

# create an 'qienp' object from a qe file
# note that if since the provided calculation is a relax I will specifiy the prefix of the file
#   this way both relax.in and relax.out will be read to create the qeinp object
filename = 'relax'
inp = pw.qeinp.from_file(filename)

# everything that we could do with atomgeo can do with qeinp

# print some things about the geometry
print(inp.ion)          # ion names
print(inp.pos)          # xyz positions
print(inp.par)          # cell parameters

# but now we have information of the qe namelist!
print(inp.nml)  # namelists include control, system, electrons, etc.

# let's make an nscf calculation
inp.nml.calculation = 'nscf'
inp.nml.nbnd = 1000
inp.kpt = [[4, 4, 4], [0, 0, 0]]

# and write the file
inp.write_file('nscf.in')


# let's refresh things
inp = pw.qeinp.from_file(filename)

# we can also see all the cards of the input file
print(inp.card)  # cards include ATOMIC_POSITIONS, ATOMIC_SPECIES, K_POINTS, etc.
print(inp.card.ATOMIC_SPECIES)

# if we replace all instances of Sn with Zr
inp.replace_ion('Sn', 'Zr')
# the ATOMIC_SPECIES will updated automatically!
print(inp.card.ATOMIC_SPECIES)

# print our full input file
print(inp)

# write the file
inp.write_file('zr.in')
