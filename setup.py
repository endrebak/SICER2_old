import os
import sys

from distutils.core import setup

from setuptools import find_packages, Extension, Command
from Cython.Build import cythonize

extensions = [Extension("SICER2.src.reads_to_bins",
                        ["SICER2/src/reads_to_bins.pyx"], language="c++"),
              Extension("SICER2.src.statistics",
                        ["SICER2/src/statistics.pyx"], language="c++"),
              Extension("SICER2.src.find_islands",
                        ["SICER2/src/find_islands.pyx"], language="c++")]

setup(ext_modules = cythonize(extensions, annotate=True))
