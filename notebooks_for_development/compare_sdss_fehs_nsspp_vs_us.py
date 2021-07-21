#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# This compares Fe/H values as calculated by nSSPP and by us

# Created 2021 July 19 by E.S.


# In[1]:


import pickle
import pandas as pd
import numpy as np
import glob
import os
import re
import matplotlib.pyplot as plt

get_ipython().run_line_magic('matplotlib', 'qt')


# In[2]:


# directory of pickled Fe/H using our abcd calibration (note just first 1k lines of posterior!)
dir_pickled_feh_abcd = "/Users/bandari/Documents/git.repos/rrlyrae_metallicity/rrlyrae_metallicity/calib_application/" +         "bin/pickled_info/escrow_us_abcd_on_sdss"


# In[3]:


# read in nSSPP values

df_nsspp = pd.read_csv("./data/nSSPP82.out", names=["sdss","spectrum", "teff", "logg", 
                                                     "feh_direct_nsspp", "feh_beers"], delim_whitespace=True)


# In[4]:


# excavate all pickle files in the directory
pickle_list = glob.glob(dir_pickled_feh_abcd + "/*.p")

# initialize data frame to hold values
# cols:
# 'inj_feh':        injected [Fe/H]
# 'err_inj_feh'     plus/minus error in injected [Fe/H]
# 'retr_med_feh'    retrieved [Fe/H]
# 'lower_sig_feh'   1-sigma lower bound of [Fe/H]
# 'upper_sig_feh'   1-sigma upper bound of [Fe/H]
# 'logg'            injected logg
# 'Teff'            injected effective temperature Teff

df = pd.DataFrame(columns=["inj_feh", "err_inj_feh", "retr_med_feh",
                            "lower_err_ret_feh", "upper_err_ret_feh", "logg", "teff", "pickle_file_name"]) #, index=range(len(pickle_list)))


for file_name in pickle_list:

    # load each item in pickle file (maybe redundant, since it is one dictionary)
    with open(file_name, "rb") as f:
        data_all = pickle.load(f)

    # calculate errors (just stdev for now, probably want to improve this later)
    feh_retrieved = np.nanmedian(data_all["feh_sample_array"])
    err_feh_retrieved = np.nanstd(data_all["feh_sample_array"])

    # add values to composite table (some keys are vestigial)
    values_to_add = {"inj_feh": data_all["injected_feh"],
                    "err_inj_feh": data_all["err_injected_feh"],
                    "logg": data_all["logg"],
                    "teff": data_all["Teff"],
                    "retr_med_feh": feh_retrieved,
                    "lower_err_ret_feh": err_feh_retrieved,
                    "upper_err_ret_feh": err_feh_retrieved,
                    "pickle_file_name": os.path.basename(file_name)}

    row_to_add = pd.Series(values_to_add, name="x")
    df = df.append(row_to_add)


# In[34]:


# loop through the nSSPP spectra names and [Fe/H] values, find matches with our retrieved [Fe/H]s,
# and put into table

df_nsspp["feh_us_abcd"] = np.nan
for nsspp_spec_num in range(0,len(df_nsspp)):

    # extract the numbers corresponding to the (obj number, group number)
    this_string = df_nsspp["spectrum"][nsspp_spec_num]
    this_obj_number = re.split('spec-|h', this_string)[1]
    this_group_number = re.split('spec-|h', this_string)[2]
    
    # find this (obj number, group number) combination in the DataFrame of retrieved metallicities
    val_interest = df[df["pickle_file_name"].str.contains(this_obj_number + "g" + this_group_number)]
    num_matches = len(val_interest)
    feh_retrieved_vals = val_interest["retr_med_feh"].values
    
    # ... if there is no match, print so
    if len(feh_retrieved_vals)==1:
        df_nsspp["feh_us_abcd"].iloc[nsspp_spec_num] = feh_retrieved_vals[0]
    elif len(feh_retrieved_vals)>1:
        print("ERROR! Too many matches")


# In[37]:


fig, ax1 = plt.subplots(1, 1)

plt.plot([np.nanmin(df_nsspp["feh_us_abcd"]),np.nanmax(df_nsspp["feh_us_abcd"])],
         [np.nanmin(df_nsspp["feh_us_abcd"]),np.nanmax(df_nsspp["feh_us_abcd"])],linestyle="--",color="gray")
plt.scatter(df_nsspp["feh_us_abcd"],df_nsspp["feh_direct_nsspp"])

# annotate
for label_num, txt in enumerate(df_nsspp["feh_us_abcd"]):
    if np.logical_and(np.isfinite(df_nsspp["feh_us_abcd"].iloc[label_num]),np.isfinite(df_nsspp["feh_us_abcd"].iloc[label_num])):
        print(str(df_nsspp["spectrum"].iloc[label_num]))
        print(df_nsspp["feh_us_abcd"].iloc[label_num],df_nsspp["feh_direct_nsspp"].iloc[label_num])
        print("---")
        '''
        ax1.annotate(str(df_nsspp["spectrum"].iloc[label_num]), 
                     xy=(df_nsspp["feh_us_abcd"].iloc[label_num],df_nsspp["feh_direct_nsspp"].iloc[label_num]), 
                     xytext=(df_nsspp["feh_us_abcd"].iloc[label_num],df_nsspp["feh_direct_nsspp"].iloc[label_num]), 
                     xycoords='data')
        '''


plt.xlabel("Retrieved [Fe/H], abcd calibration, 1-to-1")
plt.ylabel("Direct nSSPP [Fe/H]")
plt.show()


# In[7]:


plt.figure(figsize=(10,5))
plt.hist(df["retr_med_feh"],bins=1000)
plt.title("Retrieved [Fe/H] (abcd only), ~2600 SDSS spectra\nMean="+str(np.nanmean(df["retr_med_feh"]))+", Median="+str(np.nanmedian(df["retr_med_feh"])))
plt.show()


# In[19]:


print("Mean: "+str(np.mean(df["retr_med_feh"])))
print("Median: "+str(np.median(df["retr_med_feh"])))


# In[ ]:


'''
import ipdb; ipdb.set_trace()

print(data_all)

# plot retrieved and injected metallicities
# matplotlib to show error bars

fig, axes = plt.subplots(2, 1, figsize=(15, 24), sharex=True)
fig.suptitle("Retrieval comparison, from MCMC file\n" + str(self.mcmc_posteriors_file))

# Fe/H difference
df["feh_diff"] = np.subtract(df["retr_med_feh"],df["inj_feh"])

# introduce scatter in x
scatter_x = 0.1*np.random.rand(len(df["inj_feh"]))
df["inj_feh_scatter"] = np.add(scatter_x,df["inj_feh"])

cmap = sns.color_palette("YlOrBr", as_cmap=True)
#cmap = sns.rocket_palette(rot=-.2, as_cmap=True)

axes[0].plot([-2.5,0.5],[-2.5,0.5],linestyle="--",color="k",zorder=0)

# underplot error bars
axes[0].errorbar(x=df["inj_feh_scatter"],y=df["retr_med_feh"],xerr=df["err_inj_feh"],yerr=df["lower_err_ret_feh"],linestyle="",color="k",zorder=1)

g_abs = sns.scatterplot(
    ax=axes[0],
    data=df,
    x="inj_feh_scatter", y="retr_med_feh",
    hue="teff", size="logg",
    edgecolor="k",
    palette=cmap, sizes=(50, 150),
    zorder=10
)
axes[0].set_ylabel("Retrieved: [Fe/H]$_{r}$")

axes[1].plot([-2.5,0.5],[0,0],linestyle="--",color="k",zorder=0)
g_diff = sns.scatterplot(
    ax=axes[1],
    data=df,
    x="inj_feh_scatter", y="feh_diff",
    hue="teff", size="logg",
    edgecolor="k",
    palette=cmap, sizes=(50, 150),
    legend=False,
    zorder=10
)
axes[1].set_ylabel("Residual: [Fe/H]$_{r}$-[Fe/H]$_{i}$")
axes[1].set_xlabel("Injected: [Fe/H]$_{i}$")
#axes[0].set_ylim([-3.,10.])
#axes[1].set_ylim([-0.45,0.8])

plt.savefig("/Users/bandari/Desktop/junk.pdf")
'''

