#!/usr/bin/env python

from setuptools import setup, find_packages


with open("README.md", 'r') as f:
    long_description = f.read()

setup(
    name='pw2py',
    version='0.0.0',
    description='python interface for reading quantum espresso input/output',
    license=None,
    long_description=long_description,
    author='Tyler J. Smart',
    author_email='tjsmart@ucsc.edu',
    url="https://gitlab.com/tjsmart/pw2py",
    packages=find_packages(),
    install_requires=[
        'f90nml',
        'lxml',
        'numpy',
        'pandas',
        'scipy'
    ],
    scripts=[
        'scripts/build_scell',
        'scripts/conv_geo',
        'scripts/cp_geo',
        'scripts/create_defect',
        'scripts/jdftx_conv',
        'scripts/linear_mix',
        'scripts/reorder_atoms',
        'scripts/scell_embed',
        'scripts/shift_xsf',
    ]
)
