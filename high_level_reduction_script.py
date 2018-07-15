from modules2 import *  # import stuff in init file
from modules2 import create_spec_realizations #norm_spec

# normalize spectra for making the calibration in the first place (no! not first step! first we need to generate synthetic spectra, and THEN normalize)
#mamluk = norm_spec.norm_spec("in.list") # create instance 
#mamluk() # call instance

#generate_synthetic_spec
create_spec_realizations.create_spec_realizations_main("spec_names.list", synthetic_out_dir)

#run_robospect

#scrape_ew_from_robo

#apply_interstellar_ca_absorption

#findHK

#read_lit_metallicities

#rescale_lit_metallicities

#run_emcee
