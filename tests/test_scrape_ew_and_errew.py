import matplotlib
matplotlib.use('Agg')

import sys, os
import configparser
import pandas as pd
import astropy

current_dir = os.path.dirname(__file__)
target_dir = os.path.abspath(os.path.join(current_dir, "../"))
sys.path.insert(0, target_dir)

# import more things with changed system path
from modules import *
from modules import scrape_ew_and_errew
from conf import *
import numpy as np
import glob

# configuration data for reduction
config_red = configparser.ConfigParser() # for parsing values in .init file
# config for reduction to find a, b, c, d
config_red.read(os.path.join(os.path.dirname(__file__), '../conf', 'config_red.ini'))


def test_Scraper():

    '''
    write_dir_test = config_red["data_dirs"]["TEST_DIR_BIN"]
    robo_dir = config_red["data_dirs"]["DIR_ROBO"]
    file_names_test = glob.glob(config_red["data_dirs"]["TEST_DIR_SRC"] + "spec_norm_final/*")
    '''

    # instantiate
    scraper_instance = scrape_ew_and_errew.Scraper(subdir = config_red["data_dirs"]["TEST_DIR_SRC"],
                                                   file_scraped_info = config_red["data_dirs"]["TEST_DIR_BIN"]+"scraper_output/"+config_red["file_names"]["SCRAPED_EW_ALL_DATA"])


    # try a single instance; does it work?
    # note the writing of files is not directly tested here
    function_state = True
    try:
        scraper_instance()
    except Exception as e:
        # e contains printable attributes of exception object
        function_state = False

    assert function_state


def test_quality_check():

    data_out_test = scrape_ew_and_errew.quality_check(
                        read_in_filename = config_red["data_dirs"]["TEST_DIR_BIN"]+"scraper_output/"+config_red["file_names"]["TEST_SCRAPED_EW_ALL_DATA"],
                        write_out_filename = config_red["data_dirs"]["TEST_DIR_BIN"]+"scraper_output/"+config_red["file_names"]["TEST_SCRAPED_EW_DATA_GOOD_ONLY"])

    # lots of checks of data types
    # note this uses .iloc[0] instead of [0], because bad rows with index 0 may
    # have been removed
    assert isinstance(data_out_test["wavel_stated_center"].iloc[0],np.float64)
    assert isinstance(data_out_test["wavel_found_center"].iloc[0],np.float64)
    assert isinstance(data_out_test["gaussianSigma"].iloc[0],np.float64)
    assert isinstance(data_out_test["gaussianAmp"].iloc[0],np.float64)
    assert isinstance(data_out_test["uncertaintyMu"].iloc[0],np.float64)
    assert isinstance(data_out_test["uncertaintySigma"].iloc[0],np.float64)
    assert isinstance(data_out_test["uncertaintyAmp"].iloc[0],np.float64)
    assert isinstance(data_out_test["priorMu"].iloc[0],np.float64)
    assert isinstance(data_out_test["priorSigma"].iloc[0],np.float64)
    assert isinstance(data_out_test["priorAmp"].iloc[0],np.float64)
    assert isinstance(data_out_test["EQW"].iloc[0],np.float64)
    assert isinstance(data_out_test["uncertaintyEQW"].iloc[0],np.float64)
    assert isinstance(data_out_test["chiSqr"].iloc[0],np.float64)
    assert isinstance(data_out_test["flags"].iloc[0],str)
    assert isinstance(data_out_test["blendGroup"].iloc[0],np.int64)
    assert isinstance(data_out_test["line_name"].iloc[0],str)
    assert isinstance(data_out_test["robolines_file_name"].iloc[0],str)
    assert isinstance(data_out_test["realization_spec_file_name"].iloc[0],str)
    assert isinstance(data_out_test["quality"].iloc[0],str)

def test_stack_spectra():

    print("input list")
    print(config_red["data_dirs"]["TEST_DIR_SRC"] + "test_input_file_list.list")
    print("read in file name")

    data_stacked_test = scrape_ew_and_errew.stack_spectra(
                            read_in_filename = config_red["data_dirs"]["TEST_DIR_BIN"]+"scraper_output/"+config_red["file_names"]["TEST_SCRAPED_EW_DATA_GOOD_ONLY"],
                            write_out_filename = config_red["data_dirs"]["TEST_DIR_BIN"]+"scraper_output/"+config_red["file_names"]["TEST_RESTACKED_EW_DATA_GOOD_ONLY"],
                            input_list = config_red["data_dirs"]["TEST_DIR_SRC"] + "test_input_file_list.list")

    print("data_stacked")
    print(data_stacked_test.keys())

    # lots of checks of data types
    # note this uses .iloc[0] instead of [0], because bad rows with index 0 may
    # have been removed

    assert isinstance(data_stacked_test["realization_spec_file_name"].iloc[0],str)
    assert isinstance(data_stacked_test["orig_spec_file_name"].iloc[0],str)
    assert isinstance(data_stacked_test["EW_Hbeta"].iloc[0],np.float64)
    assert isinstance(data_stacked_test["err_EW_Hbeta_from_robo"].iloc[0],np.float64)
    assert isinstance(data_stacked_test["EW_Hdelta"].iloc[0],np.float64)
    assert isinstance(data_stacked_test["err_EW_Hdelta_from_robo"].iloc[0],np.float64)
    assert isinstance(data_stacked_test["EW_Hgamma"].iloc[0],np.float64)
    assert isinstance(data_stacked_test["err_EW_Hgamma_from_robo"].iloc[0],np.float64)
    assert isinstance(data_stacked_test["EW_Heps"].iloc[0],np.float64)
    assert isinstance(data_stacked_test["err_EW_Heps_from_robo"].iloc[0],np.float64)
    assert isinstance(data_stacked_test["EW_CaIIK"].iloc[0],np.float64)
    assert isinstance(data_stacked_test["err_EW_CaIIK_from_robo"].iloc[0],np.float64)


def test_generate_net_balmer():

    # generate the fake data: H_del =

    params_data, params_1to1, data_net_balmer_test = scrape_ew_and_errew.generate_net_balmer(
                                                    read_in_filename = config_red["data_dirs"]["TEST_DIR_BIN"]+"scraper_output/test_stacked_data_pre_net_balmer_calc.csv",
                                                    write_out_filename = config_red["data_dirs"]["TEST_DIR_BIN"]+"scraper_output/test_stacked_data_post_net_balmer_calc.csv")

    # extract parameters of the fit to the data
    m_test = params_data[0]
    err_m_test = params_data[1]
    b_test = params_data[2]
    err_b_test = params_data[3]

    # extract parameters of the fit to the rescaled lines
    m_test_1to1 = params_1to1[0]
    err_m_test_1to1 = params_1to1[1]
    b_test_1to1 = params_1to1[2]
    err_b_test_1to1 = params_1to1[3]

    # was the best-fit line fit correctly?
    assert round(m_test, 3) == 4.100
    assert round(b_test, 3) == 3.700
    # is the best fit to Hdelta and the rescaled Hgamma really a 1-to-1 line?
    assert round(m_test_1to1, 3) == 1.000
    assert round(b_test_1to1, 3) == 0.000
    # check data type of newly-added data
    assert isinstance(data_net_balmer_test["EW_Balmer"].iloc[0],np.float64)
    assert isinstance(data_net_balmer_test["err_EW_Balmer_based_Robo"].iloc[0],np.float64)


def test_generate_addl_ew_errors():
    # placeholder for now, until more decisions about how to calculate EW errors

    #dummy = scrape_ew_and_errew.generate_addl_ew_errors()

    assert 1<2


def add_synthetic_meta_data():
    # placeholder for now, until more decisions about how to calculate EW errors

    #dummy = scrape_ew_and_errew.generate_addl_ew_errors()

    assert 1<2
