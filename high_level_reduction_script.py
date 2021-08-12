'''
This is the high-level script which runs all the pieces of the pipeline to
obtain updated Layden coefficients [a, b, c, d]
'''

import sys
from conf import *
from modules import *
from modules import (compile_normalization,
                      create_spec_realizations,
                      run_robo,
                      scrape_ew_and_errew,
                      make_high_res_feh_basis,
                      ca_correction,
                      consolidate_pre_mcmc,
                      run_emcee)

def main():

    # make all the directories

    make_dirs(objective = "find_abcd") ## find_abcd as opposed to apply_abcd

    # compile the C spectral normalization script
    compile_normalization.compile_bkgrnd()
    '''
    # Take list of unnormalized empirical spectra and noise-churned the
    # spectra, normalize them, and write out normalizations
    ## ## just 1 or 2 realizations for testing (default is 100)
    create_spec_realizations.create_spec_realizations_main(num = 1, noise_level="None", spec_file_type="ascii.no_header")

    # run_robospect on normalized synthetic spectra
    run_robo.main()

    # scrape_ew_from_robo and calculate EWs + err_EW

    scraper_instance = scrape_ew_and_errew.Scraper()
    scraper_instance() # call instance
    scrape_ew_and_errew.quality_check()

    # put the good EW data into a table with
    # rows corresponding to files and cols for the lines
    scrape_ew_and_errew.stack_spectra(objective="find_abcd")

    # run_emcee with input data_table_winnowed
    # coeff defs: K = a + bH + cF + dHF + f(H^2) + g(F^2) + h(H^2)F + kH(F^2) + m(H^3) + n(F^3)
    # where K is CaII K EW; H is Balmer EW; F is [Fe/H]
    emcee_instance = run_emcee.RunEmcee()
    #emcee_instance(model = 'abcd') # call instance
    emcee_instance(model = 'abcdfghk')

    run_emcee.corner_plot(coeffs = 'abcdfghk')
    '''

# entry point
if __name__ == '__main__':
    sys.exit(main())
