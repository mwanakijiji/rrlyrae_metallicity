import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os, os.path
from os import listdir
from os.path import isfile, join
import pandas as pd
import sys
from pylab import * 
import glob
from IPython.display import clear_output
from astropy.io import fits
import pickle
from rrlyrae_metallicity.modules2 import *


def graft_feh(pickle_source_dir = config["data_dirs"]["DIR_PICKLE"]):
    ## ## TACK PHASES ONTO LIST OF EWS FROM SPECTRA
    ## ## NEED TO GET RID OF THE 'FAKE' AT SOME POINT
    
    # read in star names first
    ## ## N.b. this is just the RRabs with RRab offsets for now
    real_data_1 = pickle.load( open(pickle_source_dir  + "info_rrab_rrab_offsets.pkl", "rb" ) )

    # arrange the data in a way we can use
    # N.b. This is NOT fake data; I'm just appropriating the old variable name
    ## ## Note the ersatz Layden errors for now; need to revisit this with values from his paper
    fake_data_1 = { "star_name": real_data_1[0]["name_star"],
                "feh_lit": real_data_1[0]["FeH_highres"],
                "feh_layden": real_data_1[0]["FeH_basis"],
                "err_feh_lit": np.zeros(len(real_data_1[0]["FeH_basis"])),
                "err_feh_layden": 0.07*np.ones(len(real_data_1[0]["FeH_basis"]))}
    fake_dataset_1 = pd.DataFrame(data=fake_data_1)

    # loop over each star to read in the calculated metallicities
    final_star_feh = pd.DataFrame(columns=["star_name","final_feh_center","final_feh_lower","final_feh_upper"])
    for t in range(0,len(fake_data_1["star_name"])):
        this_star = fake_data_1["star_name"][t]
        name_star_underscore = str(this_star).replace(" ", "_") # replace space with underscore
        pickle_read_name = pickle_source_dir + "plot_info_" + name_star_underscore + ".pkl" # read the mapped Fe/H values
        with open(pickle_read_name, 'rb') as f:
                name_star,feh_mapped_array,x_vals,y_vals,xvals_interp,cdf_gauss_info,\
                  idx,idx_1sig_low,idx_1sig_high,shortest_xrange_lower,\
                  shortest_xrange_upper,shortest_xrange_halfway = pickle.load(f)
        this_feh_center = shortest_xrange_halfway
        this_feh_lower = shortest_xrange_lower
        this_feh_upper = shortest_xrange_upper
            
        final_star_feh = final_star_feh.append({"star_name_underscore": name_star_underscore,
                                               "final_feh_center": this_feh_center,
                                               "final_feh_lower": this_feh_lower,
                                               "final_feh_upper": this_feh_upper},
                                               ignore_index=True)

    # read in the EW and phase info
    hk_ews = pd.read_csv(config["data_dirs"]["DIR_SRC"]
                         + config["file_names"]["MORE_REALISTIC"])

    # paste the feh values onto the HK table
    # loop over each row of the HK table and assign an FeH based on string in empirical spectrum name
    hk_ews["final_feh_center"] = np.nan
    hk_ews["final_feh_lower"] = np.nan
    hk_ews["final_feh_upper"] = np.nan

    # loop over each star name (of which our program stars are a subset) and paste the FeH values to
    # the HK table rows corresponding to the empirical spectra for that star
    for star_num in range(0,len(final_star_feh["star_name_underscore"])):
        this_star = final_star_feh["star_name_underscore"][star_num]
        print("Retrieving calculated Fe/H value for " + this_star)
        feh_center_this_star = final_star_feh["final_feh_center"][star_num]
        feh_lower_this_star = final_star_feh["final_feh_lower"][star_num]
        feh_upper_this_star = final_star_feh["final_feh_upper"][star_num]

        # loop over each of our program stars; i.e., empirical spectra
        for em_spec_num in range(0,len(hk_ews["empir_spec_name"])):
            # if the star assigned to an FeH value appears in the empirical spectrum name
            if (this_star in hk_ews["empir_spec_name"][em_spec_num]):
                hk_ews["final_feh_center"][em_spec_num] = feh_center_this_star
                hk_ews["final_feh_lower"][em_spec_num] = feh_lower_this_star
                hk_ews["final_feh_upper"][em_spec_num] = feh_upper_this_star

    print(hk_ews)

    # fyi
    #hk_ews.to_csv('junk.csv')
    
    # pickle the table of H,K,phases,Fe/H
    ## ## NEED TO ADD STAR TYPE, TOO
    pickle_write_name = pickle_source_dir + "hk_final_feh_info.pkl"
    with open(pickle_write_name, "wb") as f:
        pickle.dump(hk_ews, f)
    
    return


def winnow(pickle_source_dir = config["data_dirs"]["DIR_PICKLE"]):
    '''
    This removes the program star spectra which are in the bad phase region
    '''

    # read in phase boundaries
    min_good, max_good = phase_regions()
    
    # restore pickle file with all the H,K data
    hk_data = pickle.load( open( pickle_source_dir + config["FILE_NAMES"]["KH_FINAL_PKL"], "rb" ) )
    #hk_data_df = pd.DataFrame(hk_data)
    print(hk_data.keys())
    
    # drop bad phases
    ## ## NOTE THAT THE DROPNA HERE SEEMS TO BE DROPPING ALL ROWS WITH ANY NANS IN IT (SOME OF THE RRC FEHS ARE NANS)
    ## ## ALSO CHECK THAT WERE NOT LOSING V535 OR V445 THROUGH SILLY NAME DIFFERENCES
    hk_data_winnowed = hk_data.where(np.logical_and(hk_data["phase"] > min_good,
                                                    hk_data["phase"] < max_good)).dropna().reset_index()

    #hk_data_winnowed_file_name = "hk_data_winnowed.csv"
    hk_data_winnowed.to_csv(config["data_dirs"]["DIR_BIN"] + config["file_names"]["KH_WINNOWED_FILE_NAME"])
    
    ## ## NEED TO WINNOW BY STAR TYPE, TOO

    return


def yadayada():
    
    return config["reduc_params"]["SMOOTH"]
