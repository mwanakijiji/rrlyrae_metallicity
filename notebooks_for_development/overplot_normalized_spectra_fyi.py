#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# This reads in spectra from a directory and overplots them, to give a cursory sense of
# how good the normalization routine was

# Created 2021 July 13 by E.S.


# In[1]:


import glob
import pandas as pd
import matplotlib.pyplot as plt



# In[2]:


stem = "/Users/bandari/Documents/git.repos/rrlyrae_metallicity/rrlyrae_metallicity/" +         "calib_application/realizations_output/norm/final/"


# In[3]:


file_list = glob.glob(stem + "*.dat*")


# In[4]:


for i in range(0,len(file_list)):
    
    print(i)
    
    df = pd.read_csv(file_list[i], names=["wavel","norm_flux"])
    
    plt.plot(df["wavel"], df["norm_flux"], alpha=0.3)
    
    plt.xlabel("Wavelength (angs)")
    plt.ylabel("Flux (normalized)")
    
plt.savefig("junk.png")


# In[ ]:




