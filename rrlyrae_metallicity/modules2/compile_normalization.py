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

            cc_path_abs = config_choice["data_dirs"]["DIR_SRC"] + "/bkgrnd.cc"
            compile_to_temp_path_abs = config_choice["data_dirs"]["DIR_BIN"] + "/bkgrnd"

            logging.info("--------------------------")
            logging.info("Compiling background normalization script...")
            bkgrnd_compile = Popen(["g++", "-o",
                                    compile_to_temp_path_abs,
                                    cc_path_abs],
                                    stdout=PIPE, stderr=PIPE)
            logging.info(get_setuptools_script_dir())
            # now move the executable to the right path (the Python setup tools path)
            # but what if I'm not using a conda environment? I won't know where to move
            # it, beyond just get_setuptools_script_dir()
            print("yada1")
            shutil.copy(compile_to_temp_path_abs,
                        get_setuptools_script_dir())
            print("yada2")

            logging.info("Binary for spectrum normalization saved to")
            logging.info(get_setuptools_script_dir())
            logging.info("--------------------------")

        try:
            # (why is removal of original file not working without throwing an error??)
            check_call(["rm", "-rf", compile_to_temp_path_abs])
        except CalledProcessError:
            print("Binary bkgrnd not removed from original location")
