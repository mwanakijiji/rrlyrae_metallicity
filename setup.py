# We will need:
# C++ (for Carrell's script)
# Python (for the wrapper in the first place)
# GSL ver1 as the libgsl0-dev package (for Robospect)
# maybe Fortran (for RW's code)

from setuptools import setup, Extension

module1 = Extension('hello_module', sources = ['wrapper/hello_module.cc'])

setup(name='RRab metallicity wrapper',
      version='1.0',
      description='For finding FeH from low-res survey spectra',
      author='Eckhart Spalding, Ron Wilhelm, Nathan De Lee, Kenneth Carrell',
      author_email='spalding at email dot arizona dot edu',
      url='https://github.com/mwanakijiji/rrlyrae_metallicity',
      license='MIT',
      ext_modules=[module1],
      zip_safe=False)
