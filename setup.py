# ==============================================================================
# Regrid W5E5 data for HighRes ISIMIP Experiments
# ------------------------------------------------------------------------------
# Copyright: 2022, Johanna Malle
# Licence: xxx
# ==============================================================================

import os
import setuptools
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy as np


with open('requirements.txt') as requirements_file:
    requirements = requirements_file.read().splitlines()


setup(
    name='W5E5_regrid',
    version='0.1',
    #test_suite='tests',
    packages=['w5e5_regrid'],
    install_requires=requirements,
    license='MIT',
    author='Dr. Johanna Malle',
    author_email='johanna.malle@slf.ch',
    description="Regrid W5E5 data for HighRes ISIMIP Experiments"
)
