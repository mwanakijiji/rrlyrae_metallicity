from subprocess import Popen,PIPE,check_call,CalledProcessError
import os
import shutil
from modules2 import *

'''
Compile spectral normalization script
'''

def compile_bkgrnd():
    _COMPILE_BKGRND = True
    if _COMPILE_BKGRND:
        if True:
            print("Compiling background normalization script...")
            bkgrnd_compile = Popen(["g++","-o","./rrlyrae_metallicity/bin/bkgrnd","./rrlyrae_metallicity/src/bkgrnd.cc"],stdout=PIPE,stderr=PIPE)
            # now move the executable to the right path... but what if I'm not using a conda environment? I won't know where to move it, beyond just get_setuptools_script_dir()
            shutil.copy("./rrlyrae_metallicity/bin/bkgrnd",get_setuptools_script_dir()) # put copy of executable in Python setup tools path
        try:
            check_call(["rm","-rf","./rrlyrae_metallicity/bin/bkgrnd"]) # (why is removal of original file not working without throwing an error??)
        except CalledProcessError:
            print("Binary bkgrnd not removed from original location")

    
