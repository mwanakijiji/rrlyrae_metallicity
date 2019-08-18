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
    file_name_list = glob.glob(robo_dir+"tmp/"+"*.dat")

    for p in file_name_list:

        print("Running Robospect on "+p)

        ## NEW COMMAND FOR 1 SPECTRUM
        ## rSpect.py -i 1 ./tmp/input_spectrum.dat -P /tmp/output_base_name --line_list ./tmp/lines.dat
        #args = ["python", robo_dir+"rSpect.py", "-i", "1", p,
        #        "-P", robo_dir+"tmp/output_base_name", "--line-list", robo_dir+"tmp/lines.dat"]
        #q = subprocess.call(args, shell=True)

        os.system("python " + robo_dir + "rSpect.py -i 1 " +
                  robo_dir + "tmp/input_spectrum.dat -P" +
                  robo_dir + "tmp/output_base_name --line_list " +
                  robo_dir + "tmp/lines.dat")

    print("Done with Robospect")
    print("-------------------")
