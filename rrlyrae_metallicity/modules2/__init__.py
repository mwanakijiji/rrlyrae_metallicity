#synthetic_in_dir = "empirical_spectra_unnorm"
#synthetic_out_dir = "synthetic_output"

goodPhaseRange = [0.05,0.90]

scrapedEWdataFilename = "scrapedEWdataFilename_test.dat"

smooth_val = int(22)

from setuptools import Distribution
from setuptools.command.install import install

# The class OnlyGetScriptPath() and function get_setuptools_script_dir() are from the setup.py script
# in the Apogee repository by jobovy
# https://github.com/jobovy/apogee/blob/master/setup.py

class OnlyGetScriptPath(install):
    def run(self):
        self.distribution.install_scripts = self.install_scripts

def get_setuptools_script_dir():
    " Get the directory setuptools installs scripts to for current python "
    dist = Distribution({'cmdclass': {'install': OnlyGetScriptPath}})
    dist.dry_run = True  # not sure if necessary
    dist.parse_config_files()
    command = dist.get_command_obj('install')
    command.ensure_finalized()
    command.run()
    return dist.install_scripts
