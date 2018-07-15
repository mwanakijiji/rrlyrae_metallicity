from modules2 import *  # import stuff in init file
from modules2 import create_spec_realizations, run_robo, scrape_ew_and_errew, make_high_res_feh_basis, ca_correction, graft_phases, run_emcee

########################
## MAKE THE SOLUTION IN THE FIRST PLACE; WORRY ABOUT APPLYING IT LATER
########################

# normalize spectra for making the calibration in the first place (no! not first step! first we need to generate synthetic spectra, and THEN normalize)
#mamluk = norm_spec.norm_spec("in.list") # create instance 
#mamluk() # call instance

# take empirical spectra and generate_synthetic_spec
create_spec_realizations.create_spec_realizations_main("spec_names.list", synthetic_out_dir)

# run_robospect on synthetic spectra
run_robo.run_robospect()

# scrape_ew_from_robo and calculate EWs + err_EW
mamluk2 = scrape_ew_and_errew.scraper() # create scraper instance
mamluk2() # call instance

# findHK
mamluk3 = scrape_ew_and_errew.findHK() # create findHK instance
mamluk3() # call instance

# apply_interstellar_ca_absorption
ca_correction.ca_corrxn("maps_EW(CaNa)_20150318.fits")

# read_lit_metallicities
mamluk4 = make_high_res_feh_basis.lit_metallicities()
mamluk4()

# assign phase values to spectra, put remaining data into giant table
data_table = graft_phases.graft_phases("spec_phases.list") # not made yet

# put data into giant table, winnow data based on phase
data_table_winnowed = graft_phases.winnow(data_table)

# run_emcee with input data_table_winnowed
mamluk5 = run_emcee.run_emcee()
mamluk5()

# yield the four coefficients with errors
