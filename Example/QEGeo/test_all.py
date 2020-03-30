#!/usr/bin/env python3

import os
import pw2py as pw

test_files = [
    'demo.jdftx',
    'test.vasp',
    'vc.in',
    'vc.out',
    'vc.vasp',
    'vc_0.in',
    'vc_cart.vasp'
]

sep = "-------------------------------\n"
print("reading all files in this directory")
print(sep)
for filename in os.listdir('.'):
    if not filename in test_files:
        continue
    print("reading {}\n".format(filename))
    geo = pw.atomgeo.from_file(filename)
    print(geo)
    print("writing {}\n".format("demo." + filename + ".qe"))
    geo.write_file("demo." + filename + ".qe")
    print("writing {}\n".format("demo." + filename + ".vasp"))
    geo.write_vasp("demo." + filename + ".vasp")
    print(sep)

print("Done!")


