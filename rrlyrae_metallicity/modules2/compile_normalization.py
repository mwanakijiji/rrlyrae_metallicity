'''
Compile spectral normalization script
'''

from subprocess import Popen, PIPE, check_call, CalledProcessError
import shutil
from . import * # read in config file, basic functions

def compile_bkgrnd(objective = "apply_abcd"):

    # choose the config file for the objective
    if (objective == "apply_abcd"):
        config_choice = config_apply
    elif (objective == "find_abcd"):
        config_choice = config
    
    _COMPILE_BKGRND = True
    if _COMPILE_BKGRND:
        if True:

            print("--------------------------")
            print("Compiling background normalization script...")
            bkgrnd_compile = Popen(["g++", "-o",
                                    config_choice["data_dirs"]["DIR_BIN"] + "bkgrnd",
                                    config_choice["data_dirs"]["DIR_SRC"] + "bkgrnd.cc"],
                                    stdout=PIPE, stderr=PIPE)
            print(get_setuptools_script_dir())
            # now move the executable to the right path (the Python setup tools path)
            # but what if I'm not using a conda environment? I won't know where to move
            # it, beyond just get_setuptools_script_dir()
            shutil.copy(config_choice["data_dirs"]["DIR_BIN"] + "bkgrnd",
                        get_setuptools_script_dir())

            print("Binary for spectrum normalization saved to")
            print(config_choice["data_dirs"]["DIR_BIN"] + "bkgrnd")
            print("--------------------------")

        try:
            # (why is removal of original file not working without throwing an error??)
            check_call(["rm", "-rf", get_setuptools_script_dir() + "bkgrnd"])
        except CalledProcessError:
            print("Binary bkgrnd not removed from original location")
