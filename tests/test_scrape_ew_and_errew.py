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

    data_out = scrape_ew_and_errew.quality_check(
                        read_in_filename = config_red["data_dirs"]["TEST_DIR_BIN"]+"scraper_output/"+config_red["file_names"]["SCRAPED_EW_ALL_DATA"],
                        write_out_filename = config_red["data_dirs"]["TEST_DIR_BIN"]+"scraper_output/"+config_red["file_names"]["SCRAPED_EW_DATA_GOOD_ONLY"])

    # lots of checks of data types
    # note this uses .iloc[0] instead of [0], because bad rows with index 0 may
    # have been removed
    assert isinstance(data_out["wavel_stated_center"].iloc[0],np.float64)
    assert isinstance(data_out["wavel_found_center"].iloc[0],np.float64)
    assert isinstance(data_out["gaussianSigma"].iloc[0],np.float64)
    assert isinstance(data_out["gaussianAmp"].iloc[0],np.float64)
    assert isinstance(data_out["uncertaintyMu"].iloc[0],np.float64)
    assert isinstance(data_out["uncertaintySigma"].iloc[0],np.float64)
    assert isinstance(data_out["uncertaintyAmp"].iloc[0],np.float64)
    assert isinstance(data_out["priorMu"].iloc[0],np.float64)
    assert isinstance(data_out["priorSigma"].iloc[0],np.float64)
    assert isinstance(data_out["priorAmp"].iloc[0],np.float64)
    assert isinstance(data_out["EQW"].iloc[0],np.float64)
    assert isinstance(data_out["uncertaintyEQW"].iloc[0],np.float64)
    assert isinstance(data_out["chiSqr"].iloc[0],np.float64)
    assert isinstance(data_out["flags"].iloc[0],str)
    assert isinstance(data_out["blendGroup"].iloc[0],np.int64)
    assert isinstance(data_out["line_name"].iloc[0],str)
    assert isinstance(data_out["robolines_file_name"].iloc[0],str)
    assert isinstance(data_out["realization_spec_file_name"].iloc[0],str)
    assert isinstance(data_out["original_spec_file_name"].iloc[0],str)
    assert isinstance(data_out["quality"].iloc[0],str)


def test_stack_spectra():

    assert 1<2
