# GENERATE OUR SOLUTION
# NEED:
# 1. FE/H VALUES FROM RESCALED HIGH-RES DATA FROM LITERATURE
# 2. SPECTRA PHASE VALUES
# 3. OUR OWN EW DATA, WITH CORRXNS FOR CA ABSORPTION

norm_spec

generate_synthetic_spec

run_robospect

scrape_ew_from_robo

apply_interstellar_ca_absorption

findHK

read_lit_metallicities

rescale_lit_metallicities

run_emcee
