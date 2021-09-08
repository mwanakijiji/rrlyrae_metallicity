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

    # assert column numbers are N_coeff + 1 (from index column)
    assert len(mcmc_sample_abcd.columns) == 5
    assert len(mcmc_sample_abcdfghk.columns) == 9


def test_rrmetal():

    #run_emcee.rrmetal()

    assert 1<2


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


def test_function_K():

    run_emcee.function_K()

    assert 1<2


def test_chi_sqd_fcn():

    run_emcee.chi_sqd_fcn()

    assert 1<2
