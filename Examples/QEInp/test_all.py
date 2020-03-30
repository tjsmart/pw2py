#!/usr/bin/env python3

import os
import pw2py as pw

sep = "-------------------------------\n"
print("reading all input files in this directory")
print(sep)
for filename in os.listdir('.'):
    if not filename.endswith('.in'):
        continue
    print("reading {}\n".format(filename))
    qeinp = pw.qeinp.from_file(filename)
    print(qeinp)
    print("writing {}\n".format(filename + ".out"))
    qeinp.write_file(filename + ".out")
    print(sep)

print("Done!")

