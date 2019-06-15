import subprocess
from subprocess import call
from subprocess import Popen
import shlex
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os, os.path
from os import listdir
from os.path import isfile, join
import pandas as pd
import sys
from pylab import * 
import glob
from IPython.display import clear_output
from astropy.io import fits


def run_robospect():
### THIS NOT WORKING!!!

    ##############################################################################
    # STEP 3: RUN ROBOSPECT ON ANY SPECTRA AND WRITE OUT EW VALUES AS *.c.dat FILES (applicable to A and B)
    ##############################################################################

    print("Running Robospect")
    
    # for applying to synthetic spectra 

    # accumulate list of filenames of normalized synthetic spectra
    fileNameList = glob.glob("../*_*.c.dat") # (or search whatever other directory the *.c.dat files are in)

    # for-loop to write out *.robolines and *.robospect files
    for p in fileNameList: 
        # default command: (-F: find all lines; -P: sets path of output files)
        # robospect -F -P rs.out example.dat
        args = ['./src/robospect', '-F', '-P', 'rs.out', p]
        q = subprocess.call(args, shell=True)

        # command we have used a lot
        #robospect -L lines.dat -C null --strict_width=16 X_Ari__10.dat


# N.b. Wrapper should catch anything that robospect tries to print to terminal 
# (3 pipes in any process: StdIn, StdOut, and StdError)
