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

class norm_spec:
    
    ##############################################################################
    # STEP 1: NORMALIZE SPECTRA (applicable to A and B)
    ##############################################################################
    
    # this is a superclass that will be inherited by other classes
    
    def __init__(self, input_file):
        self.smoothing = 22 # smoothing applied by Carrell's normalization code
        self.input_file = input_file
        
    def __call__(self):
        
        # compile Carrell's normalization code
        # (see carrell_readme.txt)
        
        # Carrell's C code should already be compiled by setup.py. But here are the manual commands for posterity:
        #normzn_compile1 = shlex.split("g++ -o bkgrnd bkgrnd.cc")
        #normzn_compile2 = subprocess.Popen(normzn_compile1) # run
        
        # run the normalization routine on the data
        normzn_run1 = shlex.split("./bkgrnd --smooth "+str(self.smoothing)+" "+self.input_file) # self.input_file can be 'in.data'
        normzn_run2 = subprocess.Popen(normzn_run1, stdout=subprocess.PIPE, stderr=subprocess.PIPE) # run and capture output
                                  
        # make list of all output files
        dir_name = "test_output/"
        list_output_files = [name for name in os.listdir(dir_name) if os.path.isfile(os.path.join(dir_name, name))]
                                  
        # divide the second column of the output files (empirical) with the third (normalization) 
        header = ['WAVELENGTH', 'FLUX'] # headers of output file that Robospect will use (see Robospect user manual)
        for filenum in range(0,len(list_output_files)):
            df = pd.read_csv('test_output/'+str(list_output_files[filenum]), delim_whitespace=True, header=None)
            df['WAVELENGTH'] = df[0]
            df['FLUX'] = np.divide(df[1],df[2]) # normalize
            df.to_csv('test_normzed_output/output.csv', columns = header, index = False, sep = ' ') # write out file 
            del df


## ## test: run normalization of raw data
do_normzn = norm_spec("input_file") # initialize class instance
# do_normzn.smoothing = 22 # can overload default smoothing
do_normzn() # call it
                                  
