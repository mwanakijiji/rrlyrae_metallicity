#!/usr/bin/env python
# coding: utf-8

# This makes plots showing the effective temperature retrievals based on synthetic spectra
# produced by R.W.

# Created from parent restacking_scraped_data.ipynb 2021 March 17 by E.S.

import pandas as pd
import os, sys
from astropy.io.fits import getdata
from configparser import ConfigParser, ExtendedInterpolation
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

current_dir = os.path.dirname(__file__)
target_dir = os.path.abspath(os.path.join(current_dir, "../"))
sys.path.insert(0, target_dir)

from . import *
from modules import teff_retrieval
from conf import *

# configuration data for reduction
config_red = ConfigParser(interpolation=ExtendedInterpolation()) # for parsing values in .init file
# config for reduction to find a, b, c, d
config_red.read(os.path.join(os.path.dirname(__file__), '../conf', 'config_red.ini'))


def test_temp_vs_balmer(df_poststack_file_name = config_red["data_dirs"]["TEST_DIR_SRC"]+config_red["file_names"]["TEST_RESTACKED_EW_DATA_W_METADATA_STANDALONE"]):


    m_test, err_m_test, b_test, err_b_test = teff_retrieval.temp_vs_balmer(df_poststack_file_name=df_poststack_file_name)

    assert round(m_test, 2) == 3530.67
    assert round(b_test, 2) == 12.34
