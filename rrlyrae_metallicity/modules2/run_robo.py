'''
Calls Robospect to find EWs of the normalized, noise-churned spectra
'''

import os
import glob
from rrlyrae_metallicity.modules2 import *


def run_robospect(norm_spec_source_dir=config["data_dirs"]["DIR_SYNTH_SPEC_NORM_FINAL"],
                  robo_dir=config["data_dirs"]["DIR_ROBO"]):
    '''
    INPUTS:
    norm_spec_source_dir: directory containing the normalized spectra
    robo_dir: directory of the robospect.py repo

    OUTPUTS:
    (writes files to disk)
    '''

    # for applying to synthetic spectra

    # accumulate list of filenames of normalized synthetic spectra
    ## ## this is the non-testing command
    #file_name_list = glob.glob(norm_spec_source_dir+"*.dat_*")
    ## ## this is the test command
    file_name_list = glob.glob(robo_dir+"tmp/"+"*.smo")

    for p in file_name_list:

        print("Running Robospect on "+p + " \n")

        # define string for output base names
        file_specific_string = p.split(".")[-2].split("/")[-1]

        os.system("python " +
                  robo_dir + "rSpect.py -i 4 " +
                  robo_dir + "tmp/input_spectrum.dat -P" +
                  robo_dir + "tmp/" + file_specific_string +
                  " --line_list " + robo_dir + "tmp/ll " +
                  "-C name null " +
                  "-D name null " +
                  "-N name boxcar " +
                  "-I range 10.0 " +
                  "-F chi_window 20.0 " +
                  "-vvvv")

    print("Done with Robospect")
    print("-------------------")
