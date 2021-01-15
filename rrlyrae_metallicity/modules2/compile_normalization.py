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

            cc_file_path_abs = config_choice["data_dirs"]["DIR_SRC"] + "/bkgrnd.cc"
            compile_to_file_path_abs = config_choice["data_dirs"]["DIR_BIN"] + "/bkgrnd"

            logging.info("--------------------------")
            logging.info("Compiling background normalization script...")
            bkgrnd_compile = Popen(["g++", "-o",
                                    compile_to_file_path_abs,
                                    cc_file_path_abs],
                                    stdout=PIPE, stderr=PIPE)

            logging.info("Binary for spectrum normalization saved to")
            logging.info(compile_to_file_path_abs)
            logging.info("--------------------------")
