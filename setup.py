# We will need:
# C++ (for Carrell's script)
# Python (for the wrapper in the first place)
# GSL ver1 as the libgsl0-dev package (for Robospect)

from setuptools import setup, Extension
from setuptools.command.install import install
import shlex
import subprocess
#from numpy.distutils.core import setup, Extension

'''
class CustomInstall(install):
    def run(self):
        command = "g++ -o bkgrnd7 bkgrnd.cc"
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
        
module1 = Extension("bkgrnd", sources = ["bkgrnd.cc"])

setup(name="RRab metallicity",
      version="1.0.0",
      description="For finding FeH from low-res survey spectra of RRab stars",
      author="Eckhart Spalding, Ron Wilhelm, Nathan De Lee, Kenneth Carrell",
      author_email="spalding at email dot arizona dot edu",
      url="https://github.com/mwanakijiji/rrlyrae_metallicity",
      license="MIT",
      include_package_data=True,
      ext_modules=[module1]
      ) # cmdclass={'install': CustomInstall}
