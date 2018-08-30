import os
from setuptools import setup #, Extension
from setuptools import Distribution
from setuptools.command.install import install
import sys
import shutil
import subprocess
import tempfile

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
    print(dist.install_scripts)
    return 

# shutil.copy('ferre',get_setuptools_script_dir()) ## copy [arg1] to directory [arg2]

if __name__ == "__main__":
    get_setuptools_script_dir()
