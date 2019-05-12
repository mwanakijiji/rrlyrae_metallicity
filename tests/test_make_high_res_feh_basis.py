import matplotlib
matplotlib.use('Agg')
from rrlyrae_metallicity.modules2 import *
from rrlyrae_metallicity.modules2 import make_high_res_feh_basis
import numpy as np
import pandas as pd
import pickle

def test_LitMetallicities(test_subdir = config["data_dirs"]["TEST_DIR_LIT_HIGH_RES_FEH"]):
    '''
    Read in fake literature metallicities
    '''

    # instantiate
    test_instance = make_high_res_feh_basis.LitMetallicities(source_dir = test_subdir)

    # obtain some read-in data
    clementini_fake, pancino_fake = test_instance.return_some_raw_data()

    
    ## test __init__() function
    # assign star names column to be index
    clementini_fake = clementini_fake.set_index("name")
    pancino_fake = pancino_fake.set_index("name")

    # test some read-in Fe/H values
    assert clementini_fake.at["RR Cet","log_eps_feI"] == 6.0
    assert pancino_fake.at["TU UMa","feh"] == -6.88
    assert pancino_fake.at["X Ari","feh"] != 4

    
    ## test matchmaker() function
    test_input_table = pd.read_csv(test_subdir + "pancino_2015_abundances.dat")
    test_basis_table = pd.read_csv(test_subdir + "layden_1994_abundances.dat")
    test_single_match = test_instance.matchmaker(input_table = test_input_table,
                                          basis_table = test_basis_table,
                                          basis_dataset_name = "test_layden",
                                          highres_dataset_name = "test_pancino")

    # the only overlapping stars should be X Ari, TU UMa, RR Cet
    test_names = ["X Ari", "TU UMa", "RR Cet"]
    assert len(test_single_match["name_star"]) == 3

    # the names should appear in the cross-matched table
    assert test_single_match["name_star"].isin(test_names).all()

    # the matched high-res Fe/H values should all be in the input table
    assert test_single_match["FeH_highres"].isin(test_input_table["feh"]).all()

    # the matched basis Fe/H values should all be in the input basis table
    assert test_single_match["FeH_basis"].isin(test_basis_table["feh"]).all()

    # check the names
    assert test_single_match["name_basis_dataset"].all() == "test_layden"
    assert test_single_match["name_highres_dataset"].all() == "test_pancino"


    ## test match_highres_w_basis() function: matches with ALL high-res studies
    test_all_matches_abs = test_instance.match_highres_w_basis(star_type = "RRab")
    # all the basis dataset names should be "layden_1994", and all the matching ones
    # should still be "pancino_2015" or "wallertein_2010", where the names here are
    # determined within the function in the actual pipeline
    assert test_all_matches_abs["name_basis_dataset"].all() == "layden_1994"
    # conversion to list to avoid indexing weirdness
    list_datasets = test_all_matches_abs["name_highres_dataset"].to_numpy()
    assert np.any(list_datasets == "pancino_2015")
    assert np.any(list_datasets == "wallerstein_2010")
    
    
def test_return_offsets_and_make_basis(test_subdir = config["data_dirs"]["TEST_DIR_LIT_HIGH_RES_FEH"],
                                       pickle_subdir = config["data_dirs"]["TEST_DIR_PICKLE"]):
    '''
    Tests
    return_offsets() function: offsets to apply to make a common scale
    test make_basis_via_offsets() function: make a metallicity basis by applying
        offsets to multiple high-res studies

    This test function basically goes through the functions called by calc_feh_program_stars()
    '''

    # cross-match all the fake RRab stars again
    test_instance = make_high_res_feh_basis.LitMetallicities(source_dir = test_subdir)
    test_rrab_matches = test_instance.match_highres_w_basis(star_type = "RRab")

    # find offset to be consistent with Chadid+ 2017
    test_rrab_offsets = make_high_res_feh_basis.return_offsets(data_postmatch = test_rrab_matches,
                                                          chadid_offset=True)

    # check offset for fake data is the same (to within N decimals) as that
    # found using independent means (see KH_check_workbook.xlsx)
    assert round(test_rrab_offsets["offset_highres_residuals"].values[0], 3) == -2.425

    test_rrab_basis_w_rrab_offsets = make_high_res_feh_basis.make_basis_via_offsets(df_to_offset = test_rrab_matches,
                                                                  df_offsets = test_rrab_offsets,
                                                                  plot_string = config["data_dirs"]["TEST_DIR_PLOTS"]+"test_rrab_w_rrab_offsets.png", make_plot = False)

    # check found slopes and y-intercepts with those found
    # using independent means (see KH_check_workbook.xlsx)
    assert round(test_rrab_basis_w_rrab_offsets["m_merged_highres"], 3) == 3.453 # slope for NON-shifted lit Fe/H vs. basis set Fe/H
    assert round(test_rrab_basis_w_rrab_offsets["b_merged_highres"], 3) == 2.026 # y-intercept for NON-shifted lit Fe/H vs. basis set Fe/H
    assert round(test_rrab_basis_w_rrab_offsets["m_merged_shifted_resid"], 3) == 3.096 # slope for residuals of shifted lit Fe/H vs. basis set Fe/H
    assert round(test_rrab_basis_w_rrab_offsets["b_merged_shifted_resid"], 3) == 1.789 # y-intercept for residuals of shifted lit Fe/H vs. basis set Fe/H
    assert round(test_rrab_basis_w_rrab_offsets["m_merged_highres_postshift"], 3) == 4.096 # slope for residuals of shifted lit Fe/H vs. basis set Fe/H
    assert round(test_rrab_basis_w_rrab_offsets["b_merged_highres_postshift"], 3) == 1.789 # y-intercept for residuals of shifted lit Fe/H vs. basis set Fe/H
    
    print("huzzah")
    print(test_rrab_basis_w_rrab_offsets)

    # map the fake RRab metallicities
    test_rrab_feh_highres_ab_offsets = np.add(test_rrab_basis_w_rrab_offsets["m_merged_highres"]*test_rrab_matches["FeH_basis"],
                                         test_rrab_basis_w_rrab_offsets["b_merged_highres"])
    
    # write out the plot of the mapping
    test_write_pickle = make_high_res_feh_basis.plot_mapping(input = test_rrab_matches,
                                                             mapped = test_rrab_feh_highres_ab_offsets,
                                                             title_string = "*TEST* RRab Fe/H mapping, w/ ab-based offsets",
                                                             plot_file_name = "rrab_w_ab_offsets_basis.png",
                                                             write_plot_subdir = config["data_dirs"]["TEST_DIR_FYI_INFO"])

    # pickle info
    pickle.dump( [test_rrab_matches, test_rrab_feh_highres_ab_offsets],
                     open( pickle_subdir + config["file_names"]["RRAB_RRAB_OFFSETS"], "wb" ) )
