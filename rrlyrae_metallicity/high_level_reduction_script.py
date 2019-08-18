'''
This is the high-level script which runs all the pieces of the pipeline
'''

import sys
from modules2 import *
from modules2 import \
     compile_normalization, \
     create_spec_realizations, \
     scrape_ew_and_errew, \
     make_high_res_feh_basis, \
     ca_correction, \
     consolidate_pre_mcmc, \
     run_emcee, \
     error_propagation_and_mapping

def main():

    # make all the directories
    '''
    make_dirs(objective = "find_abcd") ## find_abcd as opposed to apply_abcd

    # compile the C spectral normalization script
    compile_normalization.compile_bkgrnd()

    # Take list of unnormalized empirical spectra and generate
    # noise-churned spectra
    ## ## just 2 realizations for testing (otherwise, 100)
    create_spec_realizations.create_spec_realizations_main(num = 2)

    # run_robospect on normalized synthetic spectra
    ## ## IMPLEMENT THE PYTHON VERSION OF ROBOSPECT WHEN ITS OUT
    #run_robo.run_robospect()

    # scrape_ew_from_robo and calculate EWs + err_EW
    scraper_instance = scrape_ew_and_errew.Scraper()
    scraper_instance() # call instance

    # find equivalent widths
    find_HK_instance = scrape_ew_and_errew.findHK()
    find_HK_instance() # call instance

    # apply offsets to Fe/H values, etc., to map Fe/H values based
    # on a basis set, and pickle results
    make_high_res_feh_basis.calc_feh_program_stars()

    # bootstrap to obtain mapped Fe/H values with errors
    error_propagation_and_mapping.feh_mapper().do()
    '''
    # graft mapped FeH values onto table of EWs
    consolidate_pre_mcmc.graft_feh()
    # remove data corresponding to bad phase values, wrong type
    consolidate_pre_mcmc.winnow_by_phase_type(remove_rrl_subtype="c")

    # run_emcee with input data_table_winnowed
    emcee_instance = run_emcee.RunEmcee()
    emcee_instance() # call instance


# entry point
if __name__ == '__main__':
    sys.exit(main())
