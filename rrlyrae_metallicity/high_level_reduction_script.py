import sys
import os
from modules2 import * 
from modules2 import \
     compile_normalization, \
     create_spec_realizations, \
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

    import ipdb; ipdb.set_trace()

    # Make all the directories
    make_dirs()
    
    # Compile the C spectral normalization script
    compile_normalization.compile_bkgrnd()
    import ipdb; ipdb.set_trace()
    
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
    '''
    create_spec_realizations.create_spec_realizations_main()
    ## ## END COMMENT
    
    # run_robospect on normalized synthetic spectra
    ## ## IMPLEMENT THE PYTHON VERSION OF ROBOSPECT WHEN ITS OUT
    run_robo.run_robospect()
    '''

    # scrape_ew_from_robo and calculate EWs + err_EW
    mamluk2 = scrape_ew_and_errew.Scraper() # create scraper instance
    print('----')
    ## ## COMMENTED OUT TO SAVE TIME BUG-CHECKING
    mamluk2() # call instance
    ## ## END COMMENT

    # findHK
    ## ## COMMENTED OUT TO SAVE TIME BUG-CHECKING
    mamluk3 = scrape_ew_and_errew.findHK() # create findHK instance
    mamluk3() # call instance
    ## ## END COMMENT
    
    # apply_interstellar_ca_absorption
    ## ## ca_correction.ca_corrxn("maps_EW(CaNa)_20150318.fits")

    # THE FOLLOWING COMMAND IS OLD
    ## make FeH basis from literature (stand-alone part of code)
    #make_high_res_feh_basis.make_basis() ## ## make output of this bit get appended in cols to file corresp to scrapedEWdataFilename

    # obtain high-res metallicities for the program stars by mapping the basis set
    make_high_res_feh_basis.calc_FeH_program_stars()
    error_propagation_and_mapping.FeHmapper().do() # actually do the mapping here
    
    # put data into giant table, winnow data based on phase
    graft_phases.graft_feh() # graft FeH values onto table 
    graft_phases.winnow() # remove data corresponding to bad phase values

    # run_emcee with input data_table_winnowed
    mamluk5 = run_emcee.RunEmcee()
    mamluk5() # call instance    


# entry point
if __name__ == '__main__':
    sys.exit(main())
