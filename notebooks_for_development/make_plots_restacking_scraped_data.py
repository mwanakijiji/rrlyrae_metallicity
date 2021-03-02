#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# This restacks data after having been scraped from Robospect, such that the
# final table has rows of spectra, and cols of absorption lines (among other things)

# Created 2021 Feb. 10 by E.S.


# In[26]:


import pandas as pd
#from astropy.io import fits
from astropy.io.fits import getdata
import matplotlib.pyplot as plt
import numpy as np


# In[2]:


# name of csv file with EWs as produced by pipeline
ew_data_file_name = "/Users/bandari/Documents/git.repos/rrlyrae_metallicity/rrlyrae_metallicity/ew_products/ew_info_good_only.csv"

# read in
df_prestack = pd.read_csv(ew_data_file_name)


# In[4]:


# stem of names of FITS files with needed data in header
fits_dir = "/Users/bandari/Documents/git.repos/rrlyrae_metallicity/rrlyrae_metallicity/realizations_output/"


# In[19]:


# make list of individual spectra for which we have EW data, and
# initialize DataFrame to hold the re-cast data

list_indiv_spectra = list(df_prestack["realization_spec_file_name"].drop_duplicates())

num_indiv_spectra = len(list_indiv_spectra)

df_poststack = pd.DataFrame(columns=["realization_spec_file_name",
                                     "original_spec_file_name",
                                     "FeH", "err_FeH",
                                     "logg", "alpha","Teff",
                                     "EW_Hbeta", "err_EW_Hbeta",
                                     "EW_Hdelta", "err_EW_Hdelta",
                                     "EW_Hgamma", "err_EW_Hgamma",
                                     "EW_Heps", "err_EW_Heps",
                                     "EW_CaIIK", "err_EW_CaIIK"], index=range(num_indiv_spectra)) # initialize

for t in range(0,num_indiv_spectra):
    # loop over all spectra realizations we have measured EWs from to populate the dataframe

    this_spectrum = list_indiv_spectra[t]

    # read in intermediary FITS file to extract values from header
    image, hdr = getdata(fits_dir + this_spectrum.split(".")[0] + ".fits", header=True, ignore_missing_end=True)

    logg = hdr["LOGG"]
    teff = hdr["TEFF"]
    alpha = hdr["ALPHA"]
    feh = hdr["FEH"]
    err_feh = 0.15

    # select data from table relevant to this spectrum
    data_this_spectrum = df_prestack.where(df_prestack["realization_spec_file_name"] == this_spectrum).dropna().reset_index()

    # extract original file name (the one from which realizations are made)
    orig_name = data_this_spectrum["original_spec_file_name"].drop_duplicates().values[0]

    # extract Balmer lines from the table of data from all the spectra
    Hbeta = data_this_spectrum["EQW"].where(data_this_spectrum["line_name"] == "Hbet").dropna().values[0]
    err_Hbeta = data_this_spectrum["uncertaintyEQW"].where(data_this_spectrum["line_name"] == "Hbet").dropna().values[0]

    Hgamma = data_this_spectrum["EQW"].where(data_this_spectrum["line_name"] == "Hgam").dropna().values[0]
    err_Hgamma = data_this_spectrum["uncertaintyEQW"].where(data_this_spectrum["line_name"] == "Hgam").dropna().values[0]

    Hdelta = data_this_spectrum["EQW"].where(data_this_spectrum["line_name"] == "Hdel").dropna().values[0]
    err_Hdelta = data_this_spectrum["uncertaintyEQW"].where(data_this_spectrum["line_name"] == "Hdel").dropna().values[0]

    Heps = data_this_spectrum["EQW"].where(data_this_spectrum["line_name"] == "Heps").dropna().values[0]
    err_Heps = data_this_spectrum["uncertaintyEQW"].where(data_this_spectrum["line_name"] == "Heps").dropna().values[0]

    CaIIK = data_this_spectrum["EQW"].where(data_this_spectrum["line_name"] == "CaIIK").dropna().values[0]
    err_CaIIK = data_this_spectrum["uncertaintyEQW"].where(data_this_spectrum["line_name"] == "CaIIK").dropna().values[0]

    # fill in that row in the dataframe
    df_poststack.iloc[t]["realization_spec_file_name"] = this_spectrum
    df_poststack.iloc[t]["original_spec_file_name"] = orig_name
    df_poststack.iloc[t]["logg"] = logg
    df_poststack.iloc[t]["Teff"] = teff
    df_poststack.iloc[t]["alpha"] = alpha
    df_poststack.iloc[t]["FeH"] = feh
    df_poststack.iloc[t]["err_FeH"] = err_feh
    df_poststack.iloc[t]["EW_Hbeta"] = Hbeta
    df_poststack.iloc[t]["err_EW_Hbeta"] = err_Hbeta
    df_poststack.iloc[t]["EW_Hdelta"] = Hdelta
    df_poststack.iloc[t]["err_EW_Hdelta"] = err_Hdelta
    df_poststack.iloc[t]["EW_Hgamma"] = Hgamma
    df_poststack.iloc[t]["err_EW_Hgamma"] = err_Hgamma
    df_poststack.iloc[t]["EW_Heps"] = Heps
    df_poststack.iloc[t]["err_EW_Heps"] = err_Heps
    df_poststack.iloc[t]["EW_CaIIK"] = CaIIK
    df_poststack.iloc[t]["err_EW_CaIIK"] = err_CaIIK


# In[38]:


# to generate a net Balmer line, make a rescaling of Hgamma
# based on Hdelta

# fit a straight line to Hgam vs Hdel
x_data = df_poststack["EW_Hdelta"].values.astype(float) # Hdel
y_data = df_poststack["EW_Hgamma"].values.astype(float) # Hgam
Hgam = np.copy(y_data)

coeff, cov = np.polyfit(x_data, y_data, 1, full=False, cov=True)
m = coeff[0]
b = coeff[1]
err_m = np.sqrt(np.diag(cov))[0]
err_b = np.sqrt(np.diag(cov))[1]

print("m:")
print(m)
print("b:")
print(b)

# generate a rescaled Hgam, call it rHgam
rHgam = np.divide(np.subtract(Hgam, b), m)


# In[39]:


## BEGIN TEST TO SEE IF RESCALING IS RIGHT
'''
x_data = df_poststack["EW_Hdelta"].values.astype(float) # Hdel
y_data = rHgam # Hgam
Hgam = np.copy(y_data)

coeff, cov = np.polyfit(x_data, y_data, 1, full=False, cov=True)
m = coeff[0]
b = coeff[1]
err_m = np.sqrt(np.diag(cov))[0]
err_b = np.sqrt(np.diag(cov))[1]

print("m:")
print(m)
print("b:")
print(b)
'''
## END TEST TO SEE IF RESCALING IS RIGHT


# In[40]:


# add column of rescaled Hgamma to DataFrame

df_poststack["EW_resc_Hgamma"] = rHgam


# In[63]:


# plot: how do Balmer lines scale with each other?
plt.clf()
plt.title("Scaling of lines with Hdelta")
plt.scatter(df_poststack["EW_Hdelta"],df_poststack["EW_Hbeta"], s=3, label="Hbeta")
plt.scatter(df_poststack["EW_Hdelta"],np.add(df_poststack["EW_Hgamma"],4), s=3, label="Hgamma+4")
plt.scatter(df_poststack["EW_Hdelta"],np.add(df_poststack["EW_Heps"],8), s=3, label="Heps+8")
#plt.ylim([0,15])
plt.xlabel("EW, Hdelta (Angstr)")
plt.ylabel("EW, non-Hdelta (Angstr)")
plt.legend()
plt.savefig("junk_balmer_rescalings.pdf")


# In[77]:


# plot: KH plot
plt.clf()
plt.title("KH plot")
plt.errorbar(df_poststack["EW_resc_Hgamma"],df_poststack["EW_CaIIK"],
             yerr=df_poststack["err_EW_CaIIK"],
             marker="o", markersize=2, mfc="k", mec="k", ecolor="gray", linestyle="")
plt.ylim([0,30])
plt.xlabel("EW, net Balmer (Angstr)")
plt.ylabel("EW, CaIIK (Angstr)")
plt.savefig("junk_KH_plot.pdf")


# In[42]:


# write out data

#df_poststack.to_csv("junk_ew_data_20200216.csv")
