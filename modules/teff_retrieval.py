#!/usr/bin/env python
# coding: utf-8

# This makes plots showing the effective temperature retrievals based on synthetic spectra
# produced by R.W.

# Created from parent restacking_scraped_data.ipynb 2021 March 17 by E.S.

import pandas as pd
from astropy.io.fits import getdata
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from . import *

# name of csv file with EWs as produced by pipeline
ew_good_data_poststack_file_name = "/Users/bandari/Documents/git.repos/rrlyrae_metallicity/ew_products/restacked_ew_w_metadata.csv"

# read in
df_poststack = pd.read_csv(ew_good_data_poststack_file_name)

def line_fit(x_data_pass, y_data_pass):
    '''
    Find line of best fit

    INPUTS:
    x_data_pass: abcissa
    y_data_pass: ordinate

    OUTPUTS:
    m:      slope
    err_m:  error in slope
    b:      y-intercept
    err_b:  error in y-intercept
    '''

    # remove the stuff outside of 6000-7250 K
    #x_data_rrl = x_data_pass.where(np.logical_and(x_data_pass>=5900,x_data_pass<=7350))
    #y_data_rrl = x_data_pass.where(np.logical_and(x_data_pass>=5900,x_data_pass<=7350))
    x_data_rrl = x_data_pass[np.where(np.logical_and(y_data_pass>=5900,y_data_pass<=7350))]
    y_data_rrl = y_data_pass[np.where(np.logical_and(y_data_pass>=5900,y_data_pass<=7350))]

    coeff, cov = np.polyfit(x_data_rrl, y_data_rrl, 1, full=False, cov=True)
    m = coeff[0]
    b = coeff[1]
    err_m = np.sqrt(np.diag(cov))[0]
    err_b = np.sqrt(np.diag(cov))[1]

    print("---------")
    print("Note stuff outside of 6000-7350 K is not being considered")
    print("m:")
    print(m)
    print("err_m:")
    print(err_m)
    print("b:")
    print(b)
    print("err_b:")
    print(err_b)

    return m, err_m, b, err_b


def temp_vs_balmer(df_poststack_file_name = config_red["data_dirs"]["DIR_EW_PRODS"]+config_red["file_names"]["RESTACKED_EW_DATA_W_METADATA"]):

    # read in data
    df_poststack = pd.read_csv(df_poststack_file_name)

    # find linear trends of {net Balmer, Hdelta, and Hgamma} EW with Teff, entire Teff range
    teff = df_poststack["teff"].values.astype(float)
    # fit a straight line: net Balmer
    ews_Balmer = df_poststack["EW_Balmer"].values.astype(float)
    m, err_m, b, err_b = line_fit(x_data_pass=ews_Balmer,y_data_pass=teff)


    # plot: how do Balmer lines scale with Teff?

    plt.clf()
    plt.title("Scaling of lines with Hdelta")
    plt.plot(ews_Balmer,np.add(np.multiply(m,ews_Balmer),b), linestyle='--')
    plt.scatter(ews_Balmer,teff, s=3)
    #plt.ylim([0,70])
    plt.ylabel("Teff (K)")
    plt.xlabel("EW (Angstr)")
    plt.title("Teff from the Balmer EW")
    plt.show()
    #plt.savefig("junk_balmer_rescalings.pdf")
