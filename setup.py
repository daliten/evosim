#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import setuptools

NAME='evosim'
VERSION = '0.0.1'
DESCRIPTION='Continuous-time branching process fast simulation package for clonal growth and dynamics'

with open("README.md", "r", encoding="utf-8") as rmf:
    LONG_DESCRIPTION = rmf.read()
    

AUTHOR =  'Dalit Engelhardt'
URL = 'https://github.com/daliten/evosim'


CLASSIFIERS = [
    "Development Status :: 4 - Beta",
    "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
    "Operating System :: OS Independent",
    "Programming Language :: Python",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3",
    "Programming Language :: Python",
    "Intended Audience :: Science/Research",
    "Topic :: Scientific/Engineering",
    "Topic :: Scientific/Engineering :: Bio-Informatics",
]

LICENSE = 'GPLv3+'

KEYWORDS = 'simulation, branching process, tau-leaping, evolution, SSA, clonal dynamics, birth-death process, Python'

PYTHON_REQ=">=3.8"
with open("install-requirements.txt") as f:
    REQUIREMENTS = f.read().splitlines()

setuptools.setup(
    name=NAME,
    version=VERSION,
    description = DESCRIPTION,
    long_description = LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    keywords=KEYWORDS,
    url = URL,
    author = AUTHOR,
    classifiers = CLASSIFIERS,
    python_requires=PYTHON_REQ,
    install_requires = REQUIREMENTS,
    package_dir={"": "src"},
    packages=setuptools.find_packages(where='src',exclude=['tests']))
