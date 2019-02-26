import sys
import os
from modules2 import *  # import stuff in init file
from modules2 import \
     compile_normalization, \
     create_spec_realizations, \
     run_robo, \
     scrape_ew_and_errew, \
     make_high_res_feh_basis, \
     ca_correction, \
     graft_phases, \
     run_emcee, \
     error_propagation_and_mapping
from subprocess import Popen,PIPE
import ipdb

########################
## MAKE THE SOLUTION IN THE FIRST PLACE; WORRY ABOUT APPLYING IT LATER
########################

def main():

    '''
    Compile the C spectral normalization script
    '''
    compile_normalization.compile_bkgrnd()
    
    '''
    Take list of unnormalized empirical spectra and generate synthetic spectra
    '''
    
    # [outdir]/norm/ contains bkgrnd output
    # [outdir]/final/ contains normalized spectra
    
    # the input_list of spectra requires
    # col [0]: spectrum filename
    # col [1]: RR type (ab, c)
    # col [2]: phase (0. to 1.)

    ## ## COMMENTED OUT TO SAVE TIME BUG-CHECKING
    #create_spec_realizations.create_spec_realizations_main(input_list="./rrlyrae_metallicity/src/spec_phases_fake.list", outdir="./synthetic_output")
    ## ## END COMMENT
    
    '''
    # run_robospect on normalized synthetic spectra
    run_robo.run_robospect()
    '''

    # scrape_ew_from_robo and calculate EWs + err_EW
    mamluk2 = scrape_ew_and_errew.scraper() # create scraper instance
    print('----')
    ## ## COMMENTED OUT TO SAVE TIME BUG-CHECKING
    #test = mamluk2() # call instance
    ## ## END COMMENT
    scrapedEWdataFilename = mamluk2.get_list() # return the name of the file containing all the EW data from Robospect

    # findHK
    mamluk3 = scrape_ew_and_errew.findHK(scrapedEWdataFilename) # create findHK instance
    mamluk3() # call instance
    reducedHKdataFilename = mamluk3.get_hk_file() # return the name of the file of H and K data points we will use to do the MCMC on
    
    # apply_interstellar_ca_absorption
    ## ## ca_correction.ca_corrxn("maps_EW(CaNa)_20150318.fits")

    # THE FOLLOWING COMMAND IS OLD
    ## make FeH basis from literature (stand-alone part of code)
    #make_high_res_feh_basis.make_basis() ## ## make output of this bit get appended in cols to file corresp to scrapedEWdataFilename

    # obtain high-res metallicities for the program stars by mapping the basis set
    ab_ab_off, ab_c_off, c_c_off, c_ab_off = make_high_res_feh_basis.MetalBasisTypeSpecific(plot_name='name_here').calc_FeH_program_stars()
    error_propagation_and_mapping.FeHmapper().map() # actually do the mapping here

    
    ipdb.set_trace()
    
    # put data into giant table, winnow data based on phase
    data_table_winnowed = graft_phases.winnow(data_table) ## ## implement once we have reliable phases
    ipdb.set_trace()

    
    
    print(reducedHKdataFilename)
    # run_emcee with input data_table_winnowed
    mamluk5 = run_emcee.run_emcee(reducedHKdataFilename)
    ipdb.set_trace()
    mamluk5() # call instance
    ipdb.set_trace()
    mcmcOutputFilename = mamluk5.get_mcmc_output() # return file name of MCMC output

    # yield the four coefficients with errors
    


# entry point
if __name__ == '__main__':
    sys.exit(main())
