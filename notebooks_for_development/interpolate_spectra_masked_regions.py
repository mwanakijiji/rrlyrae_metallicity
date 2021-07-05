#!/usr/bin/env python
# coding: utf-8

# In[1]:


# This reads in UN-normalized SDSS spectra and masks and 
# 1.) Interpolates masked regions which are outside absorption lines (to preserve normalization)
# 2.) Masks regions which are located inside lines completely (to avoid messing with the EWs; note there
#      should not be many such cases--- most spectra with cosmic rays in the lines were discarded)

# Created 2021 June 27 by E.S.


# In[11]:


import pandas as pd
#from astropy.io import fits
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import glob
import os

get_ipython().run_line_magic('matplotlib', 'qt')


# In[2]:


stem = "/Users/bandari/Documents/git.repos/rrlyrae_metallicity/rrlyrae_metallicity/sdss_spectra_cosmic_ray_removal/"


# In[3]:


# glob in list of unnormalized files to use

unnorm_file_list = glob.glob(stem + "00_originals_pre_any_normalization/" + "*dat")


# In[4]:


# glob in list of masks

mask_file_list = glob.glob(stem + "01f_masks/" + "*dat*")


# In[5]:


# for each unnormalized spectrum, find the matching mask

dict_files_masks = {"spec_file_name": unnorm_file_list}

df_files_masks = pd.DataFrame.from_dict(dict_files_masks)
df_files_masks["mask_file_name"] = ""

mask_file_array = np.asarray(mask_file_list)

print("Matching masks with spectral files...")

# loop over spectral files to fill in corresponding masks
for file_num in range(0,len(df_files_masks)):
    
    # if spec name is in mask name
    this_spec = os.path.basename(df_files_masks["spec_file_name"].iloc[file_num])
    mask_exists = np.flatnonzero(np.core.defchararray.find(mask_file_array,this_spec)!=-1)
    
    if (mask_exists.size > 0):
        df_files_masks["mask_file_name"].iloc[file_num] = os.path.basename(mask_file_array[mask_exists][0])


# In[79]:


# loop over all spectra
#start = 765
for spec_num in range(10,16):#000):#len(df_files_masks)):
    
    cond_count = 0 # FYI; initialize to count the two types of artifacts: 2= inside and outside absorption lines

    # read in the spectrum and the mask
    print("Spectrum:")
    print(df_files_masks["spec_file_name"][spec_num])
    this_spectrum = pd.read_csv(df_files_masks["spec_file_name"][spec_num], delim_whitespace = True, 
                                names=["wavel","unnorm_flux","noise"])

    # if mask exists (if not, the spectrum was discarded entirely)
    try:
        this_mask = pd.read_csv(stem + "01f_masks/" + df_files_masks["mask_file_name"][spec_num])
    except:
        continue

    # where masked pixels are INSIDE absorption lines, remove wavelength, flux, and noise values from the 
    # data entirely

    # delineate where the absorption lines are
    '''
    Recall
    3933.66-30 # CaII-K
    3970.075 # H-eps
    4101.71 # H-del
    4340.472 # H-gam
    4861.29 # H-beta
    '''
    caii_K_line = np.logical_and(this_spectrum["wavel"] >= 3933.66-30,this_spectrum["wavel"] <= 3933.66+30)
    h_eps_line = np.logical_and(this_spectrum["wavel"] >= 3970.075-30,this_spectrum["wavel"] <= 3970.075+30)
    h_del_line = np.logical_and(this_spectrum["wavel"] >= 4101.71-30,this_spectrum["wavel"] <= 4101.71+30)
    h_gam_line = np.logical_and(this_spectrum["wavel"] >= 4340.472-30,this_spectrum["wavel"] <= 4340.472+30)
    h_beta_line = np.logical_and(this_spectrum["wavel"] >= 4861.29-30,this_spectrum["wavel"] <= 4861.29+30)

    # sum across the arrays
    sum_array = np.sum([np.array(caii_K_line),
                        np.array(h_eps_line),
                        np.array(h_del_line),
                        np.array(h_gam_line),
                        np.array(h_beta_line)],axis=0)
    # convert to boolean (True = 'there is an absorption line here')
    line_bool_array = np.array(sum_array, dtype=bool)
    # inversion to denote regions OUTSIDE lines
    outside_line_bool_array = ~line_bool_array

    #matches_inside_line = 
    #print(np.where(matches_inside_line == True))
    
    # indices of cosmic ray artifacts inside absorption lines...
    idx_2_drop = this_mask.index[np.logical_and(this_mask["flux_flag_1"],line_bool_array)].tolist()
    # ... and outside absorption lines
    idx_2_interp = this_mask.index[np.logical_and(this_mask["flux_flag_1"],outside_line_bool_array)].tolist()

    # if there are cosmic ray artifacts in an absorption line, remove these rows from the spectrum
    if idx_2_drop:
        this_spectrum_dropped = this_spectrum.drop(axis=0, index=idx_2_drop)
        #print("Artifacts found INside absorption lines.")
        print("Dropped index ")
        print(idx_2_drop)
        cond_count+=1
        
    else:
        # note that nothing has actually been dropped in this case
        this_spectrum_dropped = this_spectrum.copy(deep=True)
        #print("No artifacts found inside absorption lines.")
        
    # if there are cosmic ray artifacts outside the absorption lines, interpolate over these
    if idx_2_interp:
        # first, drop the rows corresponding to the artifact
        this_spectrum_dropped = this_spectrum_dropped.drop(axis=0, index=idx_2_interp)
        # then do the interpolation, of the flux and the noise
        wavel_interp = this_mask["wavel"].loc[idx_2_interp]
        flux_interp = np.interp(x=wavel_interp, 
                                xp=this_spectrum_dropped["wavel"], 
                                fp=this_spectrum_dropped["unnorm_flux"])
        noise_interp = np.interp(x=wavel_interp, 
                        xp=this_spectrum_dropped["wavel"], 
                        fp=this_spectrum_dropped["noise"])
        
        # append the interpolated points
        dict_2_append = {"wavel":wavel_interp,"unnorm_flux":flux_interp,"noise":noise_interp}
        df_2_append = pd.DataFrame.from_dict(dict_2_append)
        this_spectrum_final = this_spectrum_dropped.append(df_2_append, ignore_index=True, verify_integrity=True)
        #print("Artifacts found OUTside absorption lines.")
        
        # sort, in case the normalization routine or Robospect are picky
        this_spectrum_final = this_spectrum_final.sort_values(by=["wavel"])
        print("Interp")
        print(df_2_append)
        cond_count+=1
        
    else:
        this_spectrum_final = this_spectrum_dropped.copy(deep=True)
        #print("No artifacts found outside absorption lines.")
        
    '''
    if cond_count==2:
        # write out cleaned (but still unnormalized spectrum)
        # (PLOT)
        plt.clf()
        plt.plot(this_mask["wavel"],100*this_mask["flux_flag_1"])
        plt.plot(this_spectrum_dropped["wavel"],np.add(20.,this_spectrum_dropped["unnorm_flux"]),label="dropped")
        plt.plot(this_spectrum["wavel"],np.add(10.,this_spectrum["unnorm_flux"]),marker="o",markersize=3,label="input")
        plt.plot(this_spectrum_final["wavel"],this_spectrum_final["unnorm_flux"],marker="o",markersize=3,label="output")
        plt.legend()
        plt.show()
    '''
    
    # write out final spectrum, and save plot
    this_spectrum_final["unnorm_flux"] = this_spectrum_final["unnorm_flux"].map(lambda x: '%.3f' % x) # clean up decimals
    this_spectrum_final["noise"] = this_spectrum_final["noise"].map(lambda x: '%.5f' % x) # clean up decimals
    write_data_name = stem + "02a_unnormalized_post_mask/" +                         os.path.basename(df_files_masks["spec_file_name"][spec_num])
    this_spectrum_final.to_csv(write_data_name, sep=" ", header=False, index=False)
    print("Wrote out processed unnormzed spectrum file to " + write_data_name)

