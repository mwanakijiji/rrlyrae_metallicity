#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# This reads in ascii files and writes out FITS files with info in the headers

# Created 2021 Feb. 9 by E.S.


# In[1]:


import os
import glob
import pyfits
import pandas as pd
from astropy.io import fits


# In[2]:


# directory containing ascii

stem = "/Users/bandari/Documents/git.repos/rrlyrae_metallicity/rrlyrae_metallicity/src/model_spectra/rrmods_all/"


# In[3]:


file_list = glob.glob(stem + "*smo")


# In[4]:


for i in range(0,len(file_list)):
    
    df = pd.read_csv(file_list[i], names=["wavelength", "flux", "error"], delim_whitespace=True)
    
    outfilename = file_list[i].split(".smo")[-2] + ".fits"
    
    # parse specs from filename, given RW convention
    name_string = os.path.basename(file_list[i].split(".smo")[-2])
    
    teff = int(name_string[0:4])
    logg = 0.1*float(name_string[4:6])
    
    # alpha enhancement
    alpha_val = 0.4
    
    if "p" in name_string:
        feh = 0.1*float(name_string[-2:])
    elif "m" in name_string:
        feh = -0.1*float(name_string[-2:])
    
    # set column names
    c1=pyfits.Column(name="wavelength", format='D', array=df["wavelength"])
    c2=pyfits.Column(name="flux", format='D', array=df["flux"])
    c3=pyfits.Column(name="error", format='D', array=df["error"])     
    cols = pyfits.ColDefs([c1, c2, c3])
    
    # make header
    hdr = pyfits.Header()
    hdr["TEFF"] = teff
    hdr.comments["TEFF"] = "Effective temperature (K)"
    hdr["LOGG"] = logg
    hdr.comments["LOGG"] = "Gravity (log(g))"
    hdr["FEH"] = feh
    hdr.comments["FEH"] = "Metallicity ([Fe/H])"
    hdr["ALPHA"] = alpha_val
    hdr.comments["ALPHA"] = "Alpha enhancement ([alpha/Fe])"
    hdr["ORG_NAME"] = os.path.basename(file_list[i])
    hdr.comments["ORG_NAME"] = "Original file name"
    hdr["COMMENT"] = "-------------"
    hdr["COMMENT"] = "Wavelength units: Angstroms"
    hdr["COMMENT"] = "Flux units: F_lambda = 1 erg sec^-1 cm^-2 A^-1"
    hdr["COMMENT"] = "(Note SDSS spectral units are 10^-17 of these)"
    hdr["COMMENT"] = "Error units: (same as flux)"
    hdr["COMMENT"] = "-------------"
    hdr["COMMENT"] = "SPECTRUM v. 2.76"
    hdr["COMMENT"] = "Gray & Corbally (1994) AJ 107(2):742"
    hdr["COMMENT"] = "https://wwwuser.oats.inaf.it/castelli/grids.html"

    # write FITS
    tbhdu = pyfits.BinTableHDU.from_columns(cols)
    hdu = pyfits.PrimaryHDU(data=df, header=hdr)
    thdulist = pyfits.HDUList([hdu,tbhdu])
    thdulist.writeto(outfilename, clobber=True)
    thdulist.close()


# In[84]:


# print column names with
# f = fits.open(fits_table_filename)
# tbdata = f[1].columns

