# We will need:
# C++ (for Carrell's script)
# Python (for the wrapper in the first place)
# GSL ver1 as the libgsl0-dev package (for Robospect)

from setuptools import setup, Extension
from setuptools.command.install import install
import shlex
import subprocess
#from numpy.distutils.core import setup, Extension

long_description = "For determining metallicities of RR Lyraes from low-res spectroscopy see `here <https://github.com/mwanakijiji/rrlyrae_metallicity>`__ for more info"

'''
class CustomInstall(install):
    def run(self):
        command = "g++ -o bkgrnd11 bkgrnd.cc"
        normzn_compile = shlex.split(command)
        process = subprocess.run(normzn_compile) #, shell=True, stderr=subprocess.STDOUT)
        #process = subprocess.Popen(normzn_compile, shell=True, stderr=subprocess.STDOUT)
        #process.wait()
        print('allo?')
        print('_')
        print('_')
        print('_')
        install.run(self)
'''

# NEEDED? 
module1 = Extension("bkgrnd13", sources = ["bkgrnd.cc"]) # if extension were to be converted into a *.so file on Mac OSX

setup(name="RRab metallicity",
      version="1.0.0",
      description="For finding FeH from low-res survey spectra of RRab stars",
      long_description=long_description,
      author="Eckhart Spalding, Ron Wilhelm, Nathan De Lee, Kenneth Carrell",
      author_email="spalding at email dot arizona dot edu",
      url="https://github.com/mwanakijiji/rrlyrae_metallicity",
      license="MIT",
      include_package_data=True,
      ext_modules=[module1]
      )
# cmdclass={'install': CustomInstall}
# install_requires=['gsl==1']
