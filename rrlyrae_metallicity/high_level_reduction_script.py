import sys
import os
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
from subprocess import Popen,PIPE
import ipdb

########################
## MAKE THE SOLUTION IN THE FIRST PLACE; WORRY ABOUT APPLYING IT LATER
########################

def main():

    # Make all the directories
    make_dirs()
    
    # Compile the C spectral normalization script
    #compile_normalization.compile_bkgrnd()

    # Take list of unnormalized empirical spectra and generate noise-churned spectra
    #create_spec_realizations.create_spec_realizations_main()
    
    # run_robospect on normalized synthetic spectra
    ## ## IMPLEMENT THE PYTHON VERSION OF ROBOSPECT WHEN ITS OUT
    #run_robo.run_robospect()

    # scrape_ew_from_robo and calculate EWs + err_EW
    #scraper_instance = scrape_ew_and_errew.Scraper() # instantiate EW file scraper
    #scraper_instance() # call instance

    # find equivalent widths
    #find_HK_instance = scrape_ew_and_errew.findHK() # instantiate Balmer, CaIIK finder
    #find_HK_instance() # call instance
    
    # apply_interstellar_ca_absorption (needed?)
    ## ## ca_correction.ca_corrxn("maps_EW(CaNa)_20150318.fits")

    # apply offsets to Fe/H values, etc., to map Fe/H values based on a basis set, and pickle results
    #make_high_res_feh_basis.calc_FeH_program_stars()
    
    # bootstrap to obtain mapped Fe/H values with errors
    #error_propagation_and_mapping.FeHmapper().do()
    
    # consolidate data to feed into the MCMC: K and Balmer EWs, mapped Fe/Hs, errors; data is also winnowed based on phase
    consolidate_pre_mcmc.graft_feh() # graft mapped FeH values onto table of EWs
    consolidate_pre_mcmc.winnow_by_phase_type(remove_rrl_subtype = "c") # remove data corresponding to bad phase values, wrong type
    
    # run_emcee with input data_table_winnowed
    emcee_instance = run_emcee.RunEmcee()
    emcee_instance() # call instance    


# entry point
if __name__ == '__main__':
    sys.exit(main())
