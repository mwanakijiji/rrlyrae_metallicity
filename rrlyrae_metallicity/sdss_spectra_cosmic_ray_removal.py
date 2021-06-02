#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# This tests algorithms to remove cosmic rays from multiepoch spectra (in particular from SDSS stripe 82
# spectra, which are too many for manual removal)

# Created 2021 May 10 by E.S.


# In[ ]:


'''
Order of operations:

1.) Flag cosmic rays in spectra that have been normalized once
2.) Mask those pixels, and move spectra with flagged pixels inside absorption lines to another directory 
3.) Inspect latter by hand (discard ones with rays inside lines, put others back in with the others)
4.) Normalize the spectra a second time 
5.) Flag rays as before, but perhaps with slightly tighter constraints
6.) Mask those pixels, and move spectra to another directory (don't bother checking ones by hand)
7.) Normalize one last time
'''


# In[1]:


import pandas as pd
import numpy as np
import glob
import sys
import os
import matplotlib.pyplot as plt
from astropy.stats import sigma_clip

get_ipython().run_line_magic('matplotlib', 'inline')


# In[2]:


# find names of spectra for which continuum has been calculated

# top-level directory for SDSS spectra cosmic ray removal
stem_s82_norm = "/Users/bandari/Documents/git.repos/rrlyrae_metallicity/rrlyrae_metallicity/" +                 "sdss_spectra_cosmic_ray_removal"

# find individual file names
file_list = glob.glob(stem_s82_norm + "/01_normalized_once/*")
# find all parent names (i.e., one name for each target, whether or not multiepoch observations were made)
parent_list = list(set([i.split("g00")[0] for i in file_list]))


# In[ ]:


def main():

    # compile the C spectral normalization script
    compile_normalization.compile_bkgrnd()

    # Take list of unnormalized empirical spectra and noise-churned the
    # spectra, normalize them, and write out normalizations
    create_spec_realizations.create_spec_realizations_main(num = 1, 
                                                           input_spec_list_dir = config["data_dirs"]["DIR_SRC"],
                                                           input_list = config["data_dirs"]["DIR_SRC"] + config["file_names"]["LIST_SPEC_PHASE"],
                                                           unnorm_empirical_spectra_dir = config["data_dirs"]["DIR_RAW_SPEC_DATA"],
                                                           unnorm_noise_churned_spectra_dir = config["data_dirs"]["DIR_SYNTH_SPEC"],
                                                           bkgrnd_output_dir = config["data_dirs"]["DIR_SYNTH_SPEC_NORM"],
                                                           final_dir = config["data_dirs"]["DIR_SYNTH_SPEC_NORM_FINAL"],
                                                           noise_level="None", 
                                                           spec_file_type="ascii.no_header")


# In[3]:


def flag_from_avg(df_empir_pass,df_avg_pass,df_median_pass,sigma_choice=1):
    '''
    Average two spectra and flag points based on their deviation from the average spectrum
    
    INPUTS:
    df_empir_pass: dataframe of empirical spectrum
    df_avg_pass: dataframe of average spectrum
    df_median_pass: dataframe of median spectrum
    sigma_choice: threshold for clipping
    '''
    
    # initialize DataFrame to return
    masked_spec = df_empir_pass.copy(deep=True)
    #masked_spec["flux_masked_1"] = masked_spec["flux"]
    
    # take difference between empirical spectrum and the AVERAGE of the AVERAGE AND MEDIAN spectrum
    # (note this preserves sign information, and (if only 2 spectra are being compared) protects against 
    # misidentification of a cosmic ray in 1 spectrum when the ray is actually in the other)
    #initialize DataFrame for taking an average of some kind
    standin_df = df_avg_pass.copy(deep=True)
    standin_df["median_flux"] = df_median_pass["median_flux"]
    # remove column of wavelengths
    print(standin_df.keys())
    standin_df = standin_df.drop(labels=["wavel"],axis=1)
    # find the mean of a mean and a median
    standin_df["mean_of_stuff"] = standin_df.mean(axis=1) # average of the columns
    
    #avg_flux = np.expand_dims(df_avg_pass["avg_flux"].values,axis=1)
    #median_flux = np.expand_dims(df_median_pass["median_flux"].values,axis=1)
    #print(np.expand_dims(avg_flux,axis=0).shape)
    #print(median_flux.shape)
    #mean_median_combo = np.mean(avg_flux,median_flux)
    masked_spec["diff"] = np.subtract(df_empir_pass["flux"],standin_df["mean_of_stuff"])
    
    # mask deviant points
    # logic: is difference in the array beyond error bounds?
    error_bound = sigma_choice*np.nanstd(masked_spec["diff"])
    logic_1 = np.greater(masked_spec["diff"],error_bound)
    masked_spec["flux_flag_1"] = logic_1 # flag these points as suspect
    
    return masked_spec, error_bound


# In[13]:


# Steps for removing cosmic rays from spectra with >1 for a given object

# Step 1:
# Average two normalized spectra (following the normalization of the raw spectra)
# Flag (#1) points that are N sigma off from the average [use cenfunc=“median”]; keep track of direction

# Step 2:
# Go back to the spectra BEFORE ANY NORMALIZATION and averaging, and remove (mask) the flagged points
# Normalize them anew
# Sigma-clip again, but use larger sigmas than in Step 1, and flag those points with another flag (#2)
# Remove (mask) those points from the spectra

# Step 3:
# Normalize one last time


# find the file names of spectra corresponding to each parent; if there is only 1, ignore; 
# if >= 2, do median comparison to flag it for cosmic rays

for t in range(0,len(parent_list)):
    
    print(t)
    
    matching = list(filter(lambda x: parent_list[t] in x, file_list))
    
    print("-------------------------")
    
    if (len(matching) == 1):
        
        print("Only one match found:")
        print(matching)
    
    elif (len(matching) >= 2):
        
        print(str(len(matching)) + " matches found:")
        print(matching)
        
        # dictionary to hold dataframes
        d = {}
        
        # intialize array to contain all fluxes
        df_dummy = pd.read_csv(matching[0], names=["wavel","flux","noise"], delim_whitespace=True)
        aggregate_flux_array = np.nan*np.ones((len(df_dummy),len(matching)))
        
        # collect spectra in single dictionary
        for p in range(0,len(matching)):
            
            df_single_p = pd.read_csv(matching[p], names=["wavel","flux","noise"], delim_whitespace=True)
            
            #plt.plot(df_single_p["wavel"],df_single_p["flux"])
            
            # sanity check that wavelength abcissa are the same
            if p==0:
                # for checking wavel abcissa is same
                wavel_initial = df_single_p["wavel"].values
            else:
                if len(np.setdiff1d(df_single_p["wavel"].values,wavel_initial) >= 1):
                    print("Hey, the wavelength abcissas are not the same!")
                    sys.exit()

            # put fluxes into array
            aggregate_flux_array[:,p] = df_single_p["flux"].values
            
        # take mean flux of all the spectra
        mean_flux_array = np.mean(aggregate_flux_array,axis=1)
        
        # cast mean spectrum data as DataFrame
        df_mean = pd.DataFrame(mean_flux_array,columns=["avg_flux"])
        df_mean["wavel"] = df_single_p["wavel"] # uses last spectrum read in
        # include median flux too (important for identifying cosmic rays when only 2 spectra are compared)
        median_flux_array = np.median(aggregate_flux_array,axis=1)
        print(median_flux_array)
        df_median = pd.DataFrame(median_flux_array,columns=["median_flux"])
        df_median["wavel"] = df_single_p["wavel"] # uses last spectrum read in
        #mean_flux_array["median_flux"] = pd.Series(median_flux_array.tolist())
        
        for p in range(0,len(matching)):
            # test each empirical spectrum against the mean, and flag points
            df_single_p = pd.read_csv(matching[p], names=["wavel","flux","noise"], delim_whitespace=True)
            flagged_empirical, limit = flag_from_avg(
                                                    df_empir_pass = df_single_p,
                                                    df_avg_pass = df_mean,
                                                    df_median_pass = df_median,
                                                    sigma_choice=5
                                                    )
            
            # if cosmic ray appears to be in an absorption line, discard the spectrum
            ## ## TBD


            '''
            plt.plot(wavel_initial,mean_flux_array,linestyle="--",color="k")
            plt.show()
            plt.clf()
            '''

            fig = plt.figure(figsize=(24,10))
            plt.plot(flagged_empirical["wavel"],np.subtract(flagged_empirical["flux_flag_1"],1),color="gray",alpha=1)
            #.axvline(x=0, ymin=0, ymax=1
            plt.plot(flagged_empirical["wavel"],flagged_empirical["diff"],label="diff")
            plt.plot(df_mean["wavel"],np.add(df_mean["avg_flux"],0.2),label="mean")
            plt.plot(flagged_empirical["wavel"],flagged_empirical["flux"],label="empirical")
            #plt.plot(df_single_p["wavel"].where(test["flux_flag_1"] == True),
            #             df_single_p["flux"].where(test["flux_flag_1"] == True),
            #         label="flagged",color="k",linewidth=4)
            plt.plot([3900,5000],[limit,limit],linestyle="--")
            plt.title(str(os.path.basename(matching[p])))
            plt.legend(loc="lower right")
            plt.savefig("plot_" + str(os.path.basename(matching[p])) + ".png",
                        facecolor="white", edgecolor='white')
            plt.clf()


# In[ ]:




