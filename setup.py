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
        'scripts/conv_geo.py',
        'scripts/cp_geo.py',
        'scripts/linear_mix.py',
    ]
)
