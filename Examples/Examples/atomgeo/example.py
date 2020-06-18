#!/usr/bin/env python3

import pw2py as pw

# create an 'atomgeo' object from a qe file
filename = 'relax.out'
geo = pw.atomgeo.from_file(filename)

# print some things about the geometry
print(geo.ion)          # ion names
print(geo.pos_units)    # units for the positions (crystal, angstrom, and bohr are supported)
print(geo.pos)          # xyz positions
print(geo.par_units)    # units for the cell parameters (ibrav supported!)
print(geo.par)          # cell parameters


# some cool things you can do
# convert positions to angstrom
geo.pos_units = 'angstrom'
print(geo.pos)
# replace all instances of Sn with Zr
geo.replace_ion('Sn', 'Zr')
print(geo.ion)
# return mass of each element
print(geo.mass())
# calculate 20 nearest neighbors of Sn120 (careful python is base 0 so I use 119)
df = geo.nearest_neighbor(119, return_type='df', N=20, include_site=True)
df.to_csv('neighbors.csv')  # write to csv file


# write some new file formats (will use extension but you can also use the ftype keyword)
geo.write_file('relax.vasp')
geo.write_file('relax.xsf')
