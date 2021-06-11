'''
Compile spectral normalization script
'''

from subprocess import Popen, PIPE#, check_call, CalledProcessError
import shutil
from . import * # read in config file, basic functions (logging)

def compile_bkgrnd(
                    compiled_bkgrnd_file_path_abs_pass = compiled_bkgrnd_file_path_abs,
                    cc_bkgrnd_file_path_abs_pass = cc_bkgrnd_file_path_abs
                    ):

    _COMPILE_BKGRND = True
    if _COMPILE_BKGRND:
        if True:

            logging.info("--------------------------")
            logging.info("Compiling background normalization script...")
            bkgrnd_compile = Popen(["g++", "-o",
                                    compiled_bkgrnd_file_path_abs_pass,
                                    cc_bkgrnd_file_path_abs_pass],
                                    stdout=PIPE, stderr=PIPE)

            logging.info("Binary for spectrum normalization saved to")
            logging.info(compiled_bkgrnd_file_path_abs_pass)
            logging.info("--------------------------")
