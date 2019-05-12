from rrlyrae_metallicity.modules2 import *
from rrlyrae_metallicity.modules2 import error_propagation_and_mapping
import matplotlib
matplotlib.use('Agg')
import numpy as np
import pandas as pd

def test_FeHplotter(test_subdir = config["data_dirs"]["TEST_DIR_LIT_HIGH_RES_FEH"]):
    '''
    Check that the plotting functions are right
    '''

    # instantiate
    test_instance = error_propagation_and_mapping.FeHplotter()

    # test cdf_fcn()
    
    # take input array of random numbers, make sure sort and cumulative addition works
    test_array = np.zeros(100) # make array WITHOUT normally-distributed numbers
    test_array[-5:] = 33
    test_array[-1:] = 100
    test_x, test_y = test_instance.cdf_fcn(test_array)

    # are elements sorted?
    test_diff = np.subtract(test_x, np.roll(test_x, 1))

    # are the CDF values really cumulative?
    # (by contrast, scipy's norm.cdf compares to a normally-distributed function)
    assert np.max(test_y[np.where(test_x == 0)]) == 0.94
    assert np.max(test_y[np.where(test_x == 33)]) == 0.98
    assert test_y[np.where(test_x == 100)][0] == 0.99



    
    ## ## test_instance.cdf_gauss()
    # check this

    # test pickle_plot_info()
    
    fake_star_string = "fake_star"
    scale_distrib = 1.0
    feh_fake_star = np.random.normal(loc = 0.0, scale = scale_distrib, size = 10000)
    feh_1_low, feh_1_50, feh_1_high, feh_2_low, feh_2_mid, feh_2_high = test_instance.pickle_plot_info(name_star = fake_star_string,
                                                                                                       feh_mapped_array = feh_fake_star,
                                                                                                       write_pickle_subdir = config["data_dirs"]["TEST_DIR_PICKLE"])
    # make the plots; these will have to be checked by eye
    test_instance.write_cdf_hist_plot(name_star = fake_star_string,
                                      read_pickle_subdir = config["data_dirs"]["TEST_DIR_PICKLE"],
                                      write_plot_subdir = config["data_dirs"]["TEST_DIR_PLOTS"])
    
    # check to make sure these give values consistent with Gaussian distributions
    # (I'm not sure what better way to check this in an automated way...)
    assert np.subtract(scale_distrib,-feh_1_low) < 0.1
    assert np.subtract(scale_distrib,feh_1_high) < 0.1
    assert feh_1_50 < 0.1
    assert np.subtract(scale_distrib,-feh_2_low) < 0.1
    assert np.subtract(scale_distrib,feh_2_high) < 0.1
    assert feh_2_mid < 0.1

    '''
    print("----------------------------------------")
    print("FeH based on median and sigma brackets:")
    print(feh_1_low)
    print(feh_1_50)
    print(feh_1_high)
    print("FeH based on shortest range:")
    print(feh_2_low)
    print(feh_2_mid)
    print(feh_2_high)
    '''
    
    test_instance.do_bootstrap()


#def test_FeHmapper():
#    '''
#    Check that
#    '''#
#
#    # instantiate
#    test_instance = error_propagation_and_mapping.FeHmapper()
#    test_instance.map_feh_one_star()

    
