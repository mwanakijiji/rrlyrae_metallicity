import sys
import os
from modules2 import * 
from modules2 import \
     compile_normalization, \
     scrape_ew_and_errew, \
     make_high_res_feh_basis, \
     error_propagation_and_mapping, \
     find_feh
from subprocess import Popen,PIPE
import ipdb

# This is the reduction pipeline for applying the Layden coefficients a, b, c, d
# to low-resolution spectra

def main():

    # Make all the directories
    make_dirs(type = "apply_abcd")
    
    # Compile the C spectral normalization script
    compile_normalization.compile_bkgrnd() ## will have to update Xcode before this works

    # Take list of unnormalized empirical science spectra and normalize them
    #normalize_simple.normalize_simple() ## will have to debug this once compile_bkgrnd() stuff fixed
    
    # run_robospect on normalized spectra
    ## ## IMPLEMENT THE PYTHON VERSION OF ROBOSPECT WHEN ITS OUT
    #run_robo.run_robospect()

    # scrape_ew_from_robo and calculate EWs + err_EW
    '''
    scraper_instance = scrape_ew_and_errew.Scraper(subdir = config_apply["data_dirs"]["DIR_ROBO_OUTPUT"],
                                                   file_scraped_info = config_apply["file_names"]["SCRAPED_SCIENCE_SPECTRA_FILE_NAME"])
    #scraper_instance() # call instance

    # find Balmer, CaIIK equivalent widths
    find_HK_instance = scrape_ew_and_errew.findHK(source_subdir = config_apply["data_dirs"]["DIR_ROBO_OUTPUT"],
                                                  hk_write_subdir = config_apply["data_dirs"]["DIR_SRC"],
                                                  plot_write_subdir = config_apply["data_dirs"]["DIR_FYI_INFO"])
    find_HK_instance() # call instance
    '''
    # find Fe/H values, sampling from the a, b, c, d posteriors and while
    # incorporating equivalent width errors
    find_feh.find_feh().sample_feh()

# entry point
if __name__ == '__main__':
    sys.exit(main())
