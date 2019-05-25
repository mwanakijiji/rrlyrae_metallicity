import sys
import os
from modules2 import * 
from modules2 import \
     compile_normalization, \
     create_spec_realizations, \
     scrape_ew_and_errew, \
     make_high_res_feh_basis, \
     error_propagation_and_mapping
from subprocess import Popen,PIPE
import ipdb

# This is the reduction pipeline for applying the Layden coefficients a, b, c, d
# to low-resolution spectra

def main():

    # Make all the directories
    make_dirs(type = "apply_abcd")
    
    # Compile the C spectral normalization script
    compile_normalization.compile_bkgrnd()

    # Take list of unnormalized empirical science spectra and normalize them
    create_spec_realizations.create_spec_realizations_main()
    
    # run_robospect on normalized spectra
    ## ## IMPLEMENT THE PYTHON VERSION OF ROBOSPECT WHEN ITS OUT
    #run_robo.run_robospect()

    # scrape_ew_from_robo and calculate EWs + err_EW
    #scraper_instance = scrape_ew_and_errew.Scraper() # instantiate EW file scraper
    #scraper_instance() # call instance

    # find equivalent widths
    #find_HK_instance = scrape_ew_and_errew.findHK() # instantiate Balmer, CaIIK finder
    #find_HK_instance() # call instance
    
    # find Fe/H values, sampling from the a, b, c, d posteriors and while
    # incorporating equivalent width errors
    #find_feh.find_feh()

# entry point
if __name__ == '__main__':
    sys.exit(main())
