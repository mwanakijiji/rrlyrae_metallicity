'''
Initialization
'''

import configparser
import os
import multiprocessing
from setuptools import Distribution
from setuptools.command.install import install

# configuration data
config = configparser.ConfigParser() # for parsing values in .init file
# config for reduction to find a, b, c, d
config.read("rrlyrae_metallicity/modules2/config.ini")
# config for applying a, b, c, d
config_apply = configparser.ConfigParser()
config_apply.read("rrlyrae_metallicity/modules2/config_apply.ini")

ncpu = multiprocessing.cpu_count()

# The class OnlyGetScriptPath() and function get_setuptools_script_dir()
# are from the setup.py script in the Apogee repository by jobovy
# https://github.com/jobovy/apogee/blob/master/setup.py

class OnlyGetScriptPath(install):
    def run(self):
        self.distribution.install_scripts = self.install_scripts

def get_setuptools_script_dir():
    '''
    Get the directory setuptools installs scripts to for current python
    '''
    dist = Distribution({'cmdclass': {'install': OnlyGetScriptPath}})
    dist.dry_run = True  # not sure if necessary
    dist.parse_config_files()
    command = dist.get_command_obj('install')
    command.ensure_finalized()
    command.run()
    return dist.install_scripts

def make_dirs(objective="apply_abcd"):
    '''
    Make directories for housing files/info if they don't already exist
    '''

    # make directories for
    # 1. reduction of spectra to find a, b, c, d (objective = "find_abcd"), or
    # 2. to apply the solution (objective = "apply_abcd"; default)
    if (objective == "apply_abcd"):
        config_choice = config_apply
    elif (objective == "find_abcd"):
        config_choice = config

    # loop over all directory paths we will need
    for vals in config_choice["data_dirs"]:
        abs_path_name = str(config_choice["data_dirs"][vals])
        print("Directory exists: " + abs_path_name)

        # if directory does not exist, create it
        if not os.path.exists(abs_path_name):
            os.makedirs(abs_path_name)
            print("Made directory " + abs_path_name)

def phase_regions():
    '''
    Read in the boundary between good and bad phase regions
    '''

    # obtain values as floats
    value1 = config.getfloat("phase", "MIN_GOOD")
    value2 = config.getfloat("phase", "MAX_GOOD")

    return value1, value2
