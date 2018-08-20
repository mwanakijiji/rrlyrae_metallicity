import sys
import os
from modules2 import *  # import stuff in init file
from modules2 import create_spec_realizations, run_robo, scrape_ew_and_errew, make_high_res_feh_basis, ca_correction, graft_phases, run_emcee
from subprocess import Popen,PIPE

########################
## MAKE THE SOLUTION IN THE FIRST PLACE; WORRY ABOUT APPLYING IT LATER
########################

def main():

    '''
    Compile spectral normalization script
    '''
    bkgrnd_compile = Popen(["g++","-o","./bin/bkgrnd","./src/bkgrnd.cc"],stdout=PIPE,stderr=PIPE)

    
    '''
    Take list of unnormalized empirical spectra and generate synthetic spectra
    '''
    
    # [outdir]/norm/ contains bkgrnd output
    # [outdir]/final/ contains normalized spectra
    
    # the input_list of spectra requires
    # col [0]: spectrum filename
    # col [1]: RR type (ab, c)
    # col [2]: phase (0. to 1.)
    
    create_spec_realizations.create_spec_realizations_main(input_list="./src/spec_phases.list", outdir="./synthetic_output")

    '''
    # run_robospect on normalized synthetic spectra
    run_robo.run_robospect()
    '''

    # scrape_ew_from_robo and calculate EWs + err_EW
    mamluk2 = scrape_ew_and_errew.scraper() # create scraper instance
    print('----')
    test = mamluk2() # call instance
    scrapedEWdataFilename = mamluk2.get_list()
    

    # findHK
    mamluk3 = scrape_ew_and_errew.findHK(scrapedEWdataFilename) # create findHK instance
    mamluk3() # call instance
    
    
    # apply_interstellar_ca_absorption
    ## ## ca_correction.ca_corrxn("maps_EW(CaNa)_20150318.fits")

    # make FeH basis
    make_high_res_feh_basis.make_basis()

    '''
    # assign phase values to spectra, put remaining data into giant table
    data_table = graft_phases.graft_phases("spec_phases.list") # not made yet

    # put data into giant table, winnow data based on phase
    data_table_winnowed = graft_phases.winnow(data_table)

    # run_emcee with input data_table_winnowed
    mamluk5 = run_emcee.run_emcee()
    mamluk5()

    # yield the four coefficients with errors
    '''


# entry point (to ensure user executes this script explicitly, as opposed to importing it)
if __name__ == '__main__':
    sys.exit(main())
