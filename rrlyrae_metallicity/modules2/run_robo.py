'''
Calls Robospect to find EWs of the normalized, noise-churned spectra
'''

import os
import glob
import multiprocessing
import logging
from rrlyrae_metallicity.modules2 import *


class RunRobo:

    def __init__(self, config_data = config):

        '''
        This just configures IO
        '''

        self.norm_spec_deposit_dir = config_data["data_dirs"]["DIR_ROBO_OUTPUT"]
        self.robo_dir = config_data["data_dirs"]["DIR_ROBO"]

    def __call__(self, file_name):

        '''
        INPUTS:
        file_name: the absolute file name of one file to run Robospect on
        norm_spec_source_dir: directory containing the normalized spectra
        robo_dir: directory of the robospect.py repo

        OUTPUTS:
        (writes files to disk)
        '''

        # for applying to synthetic spectra

        logging.info("Running Robospect on "+ file_name + " \n")

        # define string for output base names
        #file_specific_string = p.split(".")[-2].split("/")[-1]
        ## ## the below is specific to *smo* files
        file_specific_string = file_name.split("/")[-1]

        # example rSpect command:
        ## python rSpect.py -i 4       ['Run four iterations of the full fitting loop,
        ##                               with noise/continuum/detection/initial line
        ##                               fits/non-linear least squares line fits all
        ##                               happening on each iteration.' --C. Waters]
        ## ./tmp/input_spectrum.dat    [run Robospect on this input spectrum]
        ## -P ./tmp/czw.20190823       [deposit results here, with this file stem]
        ## --line_list /tmp/ll         [use this line list]
        ## -C name null                [consider the continuum fixed at 1]
        ## -D name null                ['Do not detect lines not listed in the line
        ##                               list.' --C. Waters]
        ## -N name boxcar              ['Generate noise estimate from boxcar smoothing.'
        ##                               --C. Waters]
        ## -I range 10.0               ['Set initial line estimate range to +/- 10AA
        ##                               when estimating FWHM value.' --C. Waters]
        ## -F chi_window 20.0          ['Fit chi^2 calculation using windows of
        ##                               +/- 20AA from the line center.' --C. Waters]
        ## -vvvv                       ['Enable all debug messages for all components.'
        ##                               --C. Waters]

        logging.info("Robospect cmd:")
        logging.info("python "+self.robo_dir + "bin/rSpect.py -i 4 " +str(file_name) +" -P " + self.norm_spec_deposit_dir + file_specific_string +" --line_list " + self.robo_dir + "tmp/ll" +" -C name null" +" -D name null" +" -N name boxcar" + " -I range 10.0" + " -F chi_window 20.0 " + "-vvvv")
        os.system("python " +
              self.robo_dir + "bin/rSpect.py -i 4 " +
              str(file_name) +
              " -P " + self.norm_spec_deposit_dir + file_specific_string +
              " --line_list " + self.robo_dir + "tmp/ll" +
              " -C name null" +
              " -D name null" +
              " -N name boxcar" +
              " -I range 10.0" +
              " -F chi_window 20.0 " +
              "-vvvv")

        logging.info("Robospect output files written to " + \
            self.norm_spec_deposit_dir + file_specific_string + "*")


def main():

    # accumulate list of filenames of normalized synthetic spectra
    ## ## note that I have put in a specific string to look for
    ## ## in the file name here; this might be a weakness later on

    pool = multiprocessing.Pool(ncpu)

    norm_spec_source_dir = config["data_dirs"]["DIR_SYNTH_SPEC_NORM_FINAL"]

    file_name_list = glob.glob(norm_spec_source_dir+"*.smo*")
    logging.info('Reading in spectra from directory')
    logging.info(norm_spec_source_dir)

    # Check to see if it is empty (if not, there is data from a previous
    # run that will inadvertently be used later)
    write_dir = config["data_dirs"]["DIR_ROBO_OUTPUT"]
    preexisting_file_list = glob.glob(write_dir + "/*", recursive=False)
    print(preexisting_file_list)
    print(len(preexisting_file_list))
    if (len(preexisting_file_list) > 0):
        logging.info("------------------------------")
        logging.info("Directory to receive Robospect output not empty!!")
        logging.info(write_dir)
        logging.info("------------------------------")
        input("Do what you want with those files, then hit [Enter]")

    # run Robospect on normalized spectra in parallel
    run_robospect_instance = RunRobo()
    pool.map(run_robospect_instance, file_name_list)

    # serial (testing only)
    #run_robospect_instance(file_name_list[0])

    logging.info("Done with Robospect")
    logging.info("-------------------")
