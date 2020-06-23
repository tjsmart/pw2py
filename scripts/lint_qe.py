#!/usr/bin/env python

import pw2py as pw
import sys


inp = pw.qeinp.from_file(sys.argv[1])
inp.write_file('lint.in')

