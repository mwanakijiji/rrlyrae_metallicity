from rrlyrae_metallicity.modules2 import *
from rrlyrae_metallicity.modules2 import scrape_ew_and_errew
import numpy as np
import pandas as pd

def test_Scraper(test_subdir = config["data_dirs"]["TEST_DIR_ROBO_OUTPUT"]):
    '''
    Read in fake *.fits.robolines files and check to make sure the data
    is scraped and collated correctly.
    '''

    '''
    As of 2019 Apr. 7, the fake files are
    
    all_good_01_102.fits.robolines (which is 'perfectly good' data)
    all_good_02_058.fits.robolines (which is 'perfectly good' data)
    bad_hbeta_wavel_99_999.fits.robolines (which has a somewhat-off found wavelength for H-beta)
    fit_center_flag_99_999.fits.robolines (which has a flag for center correction)
    fit_max_iter_99_999.fits.robolines (which has a flag for having reached max iterations)
    fit_fail_flag_99_999.fits.robolines (which has a 'failed' flag for one line)

    ... and if all goes correctly, the 'good' data should involve the data in
    all the tables except for that in the fit_fail_flag_99_999.fits.robolines
    '''

    # instantiate
    test_instance = scrape_ew_and_errew.Scraper(subdir = test_subdir, verbose = True)

    # return the relevant tables which get pickled and used in the pipeline
    # (not really important here, since we're reading the csvs anyway)
    all_ew_data, good_ew_data_only = test_instance()

    # read in the 'good' csvs
    test_good_1 = pd.read_csv(test_subdir + "all_good_01_102.fits.robolines", header=13, delim_whitespace=True, index_col=False, usecols=np.arange(17))
    test_good_2 = pd.read_csv(test_subdir + "all_good_02_058.fits.robolines", header=13, delim_whitespace=True, index_col=False, usecols=np.arange(17))
    test_good_3 = pd.read_csv(test_subdir + "bad_hbeta_wavel_99_999.fits.robolines", header=13, delim_whitespace=True, index_col=False, usecols=np.arange(17))
    test_good_4 = pd.read_csv(test_subdir + "fit_center_flag_99_999.fits.robolines", header=13, delim_whitespace=True, index_col=False, usecols=np.arange(17))
    
    # read in the bad ones
    test_bad_1 = pd.read_csv(test_subdir + "fit_fail_flag_99_999.fits.robolines", header=13, delim_whitespace=True, index_col=False, usecols=np.arange(17))
    test_bad_2 = pd.read_csv(test_subdir + "fit_max_iter_99_999.fits.robolines", header=13, delim_whitespace=True, index_col=False, usecols=np.arange(17))
    
    # any unique bits from the 'good' csvs should appear
    # among the 'good' data
    for i in range(0,len(test_good_1["dmean"][:])):
        assert np.any(good_ew_data_only["dmean"].isin([test_good_1["dmean"][i]]))
        assert np.any(good_ew_data_only["mean"].isin([test_good_2["mean"][i]]))
        assert np.any(good_ew_data_only["EQW"].isin([test_good_3["EQW"][i]]))
        assert np.any(good_ew_data_only["eta"].isin([test_good_4["eta"][i]]))

    # any unique bits from the 'bad' csvs should not appear
    for j in range(0,len(test_good_1["dmean"][:])):
        assert np.any(good_ew_data_only["flux"].isin([test_bad_1["flux"][j]])) == False
        assert np.any(good_ew_data_only["EQW"].isin([test_bad_2["EQW"][j]])) == False
