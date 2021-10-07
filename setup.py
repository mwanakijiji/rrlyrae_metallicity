# We will need:
# C++ (for Carrell's script)
# Python (for the wrapper in the first place)
# GSL ver1 as the libgsl0-dev package (for Robospect)

from setuptools import setup, Extension
from setuptools.command.install import install
import sys
import shlex
import subprocess
from subprocess import Popen,PIPE

long_description = "For determining metallicities of RR Lyraes from low-res spectroscopy see `here <https://github.com/mwanakijiji/rrlyrae_metallicity>`__ for more info"

setup(name="rrlyrae_metallicity",
      version="0.0.3",
      description="For finding FeH from low-res survey spectra of RRab stars",
      long_description=long_description,
      author="Eckhart Spalding, Ron Wilhelm, Nathan De Lee, Kenneth Carrell",
      author_email="espaldin@nd.edu",
      url="https://github.com/mwanakijiji/rrlyrae_metallicity",
      license="MIT",
      packages=['modules'],
      include_package_data=True
      )
