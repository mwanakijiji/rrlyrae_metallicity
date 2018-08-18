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

# check --compile-background flag
# if true, compile Carrell's normalization script

'''
try:
    compile_back = sys.argv.index('--compile-background')
except ValueError:
    _COMPILE_BACKGROUND = False
else:
    del sys.argv[compile_back]
    _COMPILE_BACKGROUND = True
if _COMPILE_BACKGROUND:
    bkgrnd_compile = Popen(["g++","-o","bkgrnd351","bkgrnd.cc"],stdout=PIPE,stderr=PIPE)
'''
    
#bkgrnd_compile = Popen(["g++","-o","bkgrnd351","bkgrnd.cc"],stdout=PIPE,stderr=PIPE)
    
setup(name="RRab metallicity",
      version="1.0.0",
      description="For finding FeH from low-res survey spectra of RRab stars",
      long_description=long_description,
      author="Eckhart Spalding, Ron Wilhelm, Nathan De Lee, Kenneth Carrell",
      author_email="spalding at email dot arizona dot edu",
      url="https://github.com/mwanakijiji/rrlyrae_metallicity",
      license="MIT",
      include_package_data=True
      )
