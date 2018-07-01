# We will need:
# C++ (for Carrell's script)
# Python (for the wrapper in the first place)
# GSL ver1 as the libgsl0-dev package (for Robospect)
# maybe Fortran (for RW's code)

from setuptools import setup, Extension

module1 = Extension('test_c', sources = ['wrapper/test.cc'])

setup(name='Metallicity wrapper',
      version='1.0',
      description='For finding FeH from low-res survey spectra',
      author='Eckhart Spalding, Ron Wilhelm, Nathan De Lee, Kenneth Carrell',
      author_email='spalding@email.arizona.edu',
      license='MIT',
      ext_modules=[module1],
      zip_safe=False)
