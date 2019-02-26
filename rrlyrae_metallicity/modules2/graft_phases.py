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

def graft_phases():
    ## ## TACK PHASES ONTO LIST OF EWS FROM SPECTRA
    ## ## NEED TO GET RID OF THE 'FAKE' AT SOME POINT

        # read in star names first
        ## ## N.b. this is just the RRabs with RRab offsets for now
        real_data_1 = pickle.load( open( "./rrlyrae_metallicity/modules2/pickled_info/info_rrab_rrab_offsets.pkl", "rb" ) )

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
        final_star_feh = pd.DataFrame(columns=["star_name","final_feh_center","final_feh_lower","final_feh_higher"])
        for t in range(0,len(fake_data_1["star_name"])):
            this_star = fake_data_1["star_name"][t]
            name_star_underscore = str(name_star).replace(" ", "_") # replace space with underscore
            pickle_read_name = "./rrlyrae_metallicity/modules2/pickled_info/plot_info_" + name_star_underscore + ".pkl" # read the mapped Fe/H values
            with open(pickle_read_name, 'rb') as f:
                name_star,feh_mapped_array,x_vals,y_vals,xvals_interp,cdf_gauss_info,\
                  idx,idx_1sig_low,idx_1sig_high,shortest_xrange_lower,\
                  shortest_xrange_upper,shortest_xrange_halfway = pickle.load(f)
            this_feh_center = shortest_xrange_halfway
            this_feh_lower = shortest_xrange_lower
            this_feh_higher = shortest_xrange_higher
            
            final_star_feh = final_star_feh.append({"star_name_underscore": name_star_underscore,
                                               "final_feh_center": this_feh_center,
                                               "final_feh_lower": this_feh_lower,
                                               "final_feh_higher": this_feh_higher}, ignore_index=True})

        # read in the EW and phase info
        hkFileName = "rrlyrae_metallicity/src/more_realistic_EWs_w_phase_test.csv"
        hk_ews = pd.read_csv(hkFileName)

        # paste the feh values onto the HK table
        # loop over each row of the HK table and assign an FeH based on string in empirical spectrum name
        hk_ews["final_feh_center"] = np.nan
        hk_ews["final_feh_lower"] = np.nan
        hk_ews["final_feh_higher"] = np.nan
        for em_spec_num in range(0,len(hk_ews["empir_spec_name"])):
            feh_values_we_want = final_star_feh.where(np.isin(hk_ews["empir_spec_name"][em_spec_num],
                                                              final_star_feh["star_name_underscore"])).values()
            hk_ews["final_feh_center"][em_spec_num] = feh_values_we_want["final_feh_center"]
            hk_ews["final_feh_lower"][em_spec_num] = feh_values_we_want["final_feh_lower"]
            hk_ews["final_feh_higher"][em_spec_num] = feh_values_we_want["final_feh_higher"]

        # pickle the table of H,K,phases,Fe/H
        
            
    
    pass


def winnowed():
    ## ## NOW REMOVE THE LINES IN BAD PHASE REGIONS


    
    pass
