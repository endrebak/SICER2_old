import os
import sys

from distutils.core import setup

from setuptools import find_packages, Extension, Command
from Cython.Build import cythonize

compile_options = ["-Ofast", "-Wall"] #, "-frename-registers", "-funroll-loops"] #
                   # -fprofile-generate
                   #"-fopenmp", "-D_GLIBCXX_PARALLEL"]

extensions = [Extension("SICER2.src.reads_to_bins",
                        ["SICER2/src/reads_to_bins.pyx"], language="c++", extra_compile_args=compile_options),
              Extension("SICER2.src.statistics",
                        ["SICER2/src/statistics.pyx"], language="c++", extra_compile_args=compile_options),
              Extension("SICER2.src.find_islands",
                        ["SICER2/src/find_islands.pyx"], language="c++", extra_compile_args=compile_options)]

setup(ext_modules = cythonize(extensions, annotate=True))
