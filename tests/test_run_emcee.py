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
from modules import run_emcee
from conf import *
import numpy as np
import glob

# configuration data for reduction
config_red = configparser.ConfigParser() # for parsing values in .init file
# config for reduction to find a, b, c, d
config_red.read(os.path.join(os.path.dirname(__file__), '../conf', 'config_red.ini'))


def test_corner_plot():

    # get a sample of the MCMC posterior data after being read in, and check the column
    # numbers are consisten with the model
    mcmc_sample_abcd = run_emcee.corner_plot(model = "abcd",
                            mcmc_text_output_file_name = config_red["data_dirs"]["TEST_DIR_SRC"] + "test_mcmc_output_abcd.csv",
                            corner_plot_putput_file_name = config_red["data_dirs"]["TEST_DIR_BIN"] + "test_abcd_plot.png")

    mcmc_sample_abcdfghk = run_emcee.corner_plot(model = "abcdfghk",
                            mcmc_text_output_file_name = config_red["data_dirs"]["TEST_DIR_SRC"] + "test_mcmc_output_abcdfghk.csv",
                            corner_plot_putput_file_name = config_red["data_dirs"]["TEST_DIR_BIN"] + "test_abcdfghk_plot.png")

    # assert column numbers are N_coeff + 1 (1 extra from index column)
    assert len(mcmc_sample_abcd.columns) == 5
    assert len(mcmc_sample_abcdfghk.columns) == 9


'''
def test_lnprob():

    run_emcee.lnprob()

    assert 1<2


def test_lnprior():

    run_emcee.lnprior()

    assert 1<2


def test_find_indices():

    run_emcee.find_indices()

    assert 1<2


def test_RunEmcee():

    run_emcee.RunEmcee()

    assert 1<2
'''


def test_function_K():

    coeffs_4_test = np.array([1.5,2.6,3.7,4.8])
    coeffs_8_test = np.array([1.5,2.6,3.7,4.8,0.5,0.32,-0.12,-0.03])
    Bal_test = 0.55
    Feh_test = -0.4

    k_4_test = run_emcee.function_K(coeffs_pass=coeffs_4_test,
                                H_pass=Bal_test,
                                F_pass=Feh_test)

    ## ## CONTINUE HERE ## ##
    k_8_test = run_emcee.function_K(coeffs_pass=coeffs_8_test,
                                H_pass=Bal_test,
                                F_pass=Feh_test)


    assert round(k_4_test, 3) == 0.394


def test_sigma_Km_sqd():

    # case of 4 coefficients
    sigma_Km_4_sqd = run_emcee.sigma_Km_sqd(
                                        coeffs_pass=np.array([2.4,5.3,4.5,4.9]),
                                        Bal_pass=6.3,
                                        err_Bal_pass=0.33,
                                        Feh_pass=-0.14,
                                        err_Feh_pass=0.11
                                        )

    # case of 8 coefficients
    sigma_Km_8_sqd = run_emcee.sigma_Km_sqd(
                                        coeffs_pass=np.array([2.4,5.3,4.5,4.9,0.5,0.32,-0.12,-0.03]),
                                        Bal_pass=6.3,
                                        err_Bal_pass=0.33,
                                        Feh_pass=-0.14,
                                        err_Feh_pass=0.11
                                        )

    assert round(sigma_Km_4_sqd, 3) == 17.456
    assert round(sigma_Km_8_sqd, 3) == 24.786

'''
def test_chi_sqd_fcn():

    # find chi-squared for 4 and 8 coeff models

    coeffs_4_test = np.array([2.4,5.3,4.5,4.9])
    coeffs_8_test = np.array([2.4,5.3,4.5,4.9,0.5,0.32,-0.12,-0.03])
    balmer_ew_test = 6.3
    err_balmer_ew_test = 0.33
    feh_test = -0.14
    err_feh_test = 0.11
    caiik_ew_test = 1.3
    err_caiik_ew_test = 0.27

    # first calculate for single values of EW
    chi_sq_4_test_i = run_emcee.chi_sqd_fcn(xi_pass=balmer_ew_test,
                            yi_pass=feh_test,
                            zi_pass=caiik_ew_test,
                            sig_xi_pass=err_balmer_ew_test,
                            sig_yi_pass=err_feh_test,
                            sig_zi_pass=err_caiik_ew_test,
                            coeffs_pass=coeffs_4_test)

    chi_sq_8_test_i = run_emcee.chi_sqd_fcn(xi_pass=balmer_ew_test,
                            yi_pass=feh_test,
                            zi_pass=caiik_ew_test,
                            sig_xi_pass=err_balmer_ew_test,
                            sig_yi_pass=err_feh_test,
                            sig_zi_pass=err_caiik_ew_test,
                            coeffs_pass=coeffs_8_test)

    print("chi sq")
    print(chi_sq_4_test_i)
    print(chi_sq_8_test_i)

    assert round(chi_sq_4_test_i, 3) == 34.350
    assert round(chi_sq_8_test_i, 3) == 100.777
'''
