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

    run_emcee.corner_plot()

    assert 1<2


def test_rrmetal():

    run_emcee.rrmetal()

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
