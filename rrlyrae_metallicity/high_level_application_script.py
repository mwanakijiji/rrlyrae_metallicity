'''
This is the reduction pipeline for applying the updated Layden coefficients
[a, b, c, d] to low-resolution spectra
'''

import sys
from modules2 import *
from modules2 import (compile_normalization,
                      create_spec_realizations,
                      scrape_ew_and_errew,
                      error_propagation_and_mapping,
                      find_feh)

def main():

    # set the objective
    # "apply_abcd": apply pre-determined coefficients to spectra to find Fe/H
    # "find_abcd": determine coefficients [a,b,c,d] in the first place
    objective_choice = "apply_abcd"

    # Make all the directories
    make_dirs(objective = objective_choice)

    # Compile the C spectral normalization script
    compile_normalization.compile_bkgrnd()

    # Take list of unnormalized empirical science spectra and normalize them
    # (N.b. only make 1 realization, with 0 noise; note also that application of
    # the calibration uses non-default directories)
    create_spec_realizations.create_spec_realizations_main(num = 1,
                                                            spec_file_type="ascii.no_header",
                                                            noise_level=0,
                                                            input_spec_list_dir = config_apply["data_dirs"]["DIR_SRC"],
                                                            input_list = config_apply["data_dirs"]["DIR_SRC"] + config_apply["file_names"]["LIST_SPEC_APPLY"],
                                                            unnorm_empirical_spectra_dir = config_apply["data_dirs"]["DIR_SCI_SPECTRA"],
                                                            unnorm_noise_churned_spectra_dir = config_apply["data_dirs"]["DIR_SYNTH_SPEC"],
                                                            bkgrnd_output_dir = config_apply["data_dirs"]["DIR_SYNTH_SPEC_NORM"],
                                                            final_dir = config_apply["data_dirs"]["DIR_SYNTH_SPEC_NORM_FINAL"])

    '''
    # run_robospect on normalized spectra
    run_robo.run_robospect()

    # scrape_ew_from_robo and calculate EWs + err_EW
    scraper_instance = scrape_ew_and_errew.Scraper(subdir = config_apply["data_dirs"]["DIR_ROBO_OUTPUT"],
                                                   file_scraped_info = config_apply["file_names"]["SCRAPED_SCIENCE_SPECTRA_FILE_NAME"])
    #scraper_instance() # call instance

    # find Balmer, CaIIK equivalent widths
    find_HK_instance = scrape_ew_and_errew.findHK(source_subdir = config_apply["data_dirs"]["DIR_ROBO_OUTPUT"],
                                                  hk_write_subdir = config_apply["data_dirs"]["DIR_SRC"],
                                                  plot_write_subdir = config_apply["data_dirs"]["DIR_FYI_INFO"])
    find_HK_instance() # call instance

    # find Fe/H values, sampling from the a, b, c, d posteriors and while
    # incorporating equivalent width errors
    find_feh.find_feh().sample_feh()
    '''

# entry point
if __name__ == '__main__':
    sys.exit(main())
