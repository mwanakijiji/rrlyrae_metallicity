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

# name of csv file with EWs as produced by pipeline
ew_good_data_poststack_file_name = "/Users/bandari/Documents/git.repos/rrlyrae_metallicity/ew_products/restacked_ew_w_metadata.csv"

# read in
df_poststack = pd.read_csv(ew_good_data_poststack_file_name)

def line_fit(x_data_pass, y_data_pass):

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

    return m, b


# plot: how do Balmer lines scale with Teff?

plt.clf()
plt.title("Scaling of lines with Hdelta")
plt.scatter(df_poststack["teff"],df_poststack["EW_Hbeta"], s=3, label="Hbeta")
plt.scatter(df_poststack["teff"],np.add(df_poststack["EW_Hgamma"],6), s=3, label="Hgamma+6")
plt.scatter(df_poststack["teff"],np.add(df_poststack["EW_Hdelta"],12), s=3, label="Hdel+12")
plt.scatter(df_poststack["teff"],np.add(df_poststack["EW_Balmer"],18), s=3, label="Net Balmer+18")
plt.scatter(df_poststack["teff"],np.add(df_poststack["EW_Heps"],24), s=3, label="Heps+24")
#plt.ylim([0,70])
plt.xlabel("Teff (K)")
plt.ylabel("EW (Angstr)")
plt.title("Balmer EW trend with Teff")
plt.legend(ncol=5)
plt.show()
#plt.savefig("junk_balmer_rescalings.pdf")

'''
y_data_metalrich = df_poststack["Teff"].where(df_poststack["FeH"] > -2.9).dropna().values.astype(float)
x_data_Balmer_metalrich = df_poststack["EW_Balmer"].where(df_poststack["FeH"] > -2.9).dropna()

x_data_Balmer_metalrich
y_data_metalrich


# find linear trends of {net Balmer, Hdelta, and Hgamma} EW with Teff, entire Teff range

y_data = df_poststack["Teff"].values.astype(float)

# fit a straight line: net Balmer
x_data_Balmer = df_poststack["EW_Balmer"].values.astype(float)
m_Balmer, b_Balmer = line_fit(x_data_Balmer,y_data)
# same, except that [Fe/H] = -3 is neglected
x_data_Balmer_metalrich = df_poststack["EW_Balmer"].where(df_poststack["FeH"] > -2.9).dropna().values.astype(float)
y_data_metalrich = df_poststack["Teff"].where(df_poststack["FeH"] > -2.9).dropna().values.astype(float)
m_Balmer_metalrich, b_Balmer_metalrich = line_fit(x_data_Balmer_metalrich,y_data_metalrich)

# fit a straight line: Hdelta
x_data_Hdelta = df_poststack["EW_Hdelta"].values.astype(float)
m_Hdelta, b_Hdelta = line_fit(x_data_Hdelta,y_data)

# fit a straight line: Hgamma
x_data_Hgamma = df_poststack["EW_Hgamma"].values.astype(float)
m_Hgamma, b_Hgamma = line_fit(x_data_Hgamma,y_data)


# calculate retrieved Teff and add new columns to DataFrame to make the plotting easier

df_poststack["Teff_retrieved_Balmer"] = np.add(np.multiply(df_poststack["EW_Balmer"],m_Balmer),b_Balmer)
df_poststack["Teff_retrieved_Hdelta"] = np.add(np.multiply(df_poststack["EW_Hdelta"],m_Hdelta),b_Hdelta)
df_poststack["Teff_retrieved_Hgamma"] = np.add(np.multiply(df_poststack["EW_Hgamma"],m_Hgamma),b_Hgamma)
df_poststack["Teff_retrieved_Balmer_metalrich"] = np.add(np.multiply(df_poststack["EW_Balmer"],m_Balmer_metalrich),b_Balmer_metalrich)


colormap = "Reds"


# array of metallicities

feh_values = np.sort(df_poststack["FeH"].drop_duplicates().values)


norm = matplotlib.colors.Normalize(vmin=np.min(feh_values),vmax=np.max(feh_values))


# retrieved Balmer values
# retrieved Balmer values
plt.clf()

colormap="Reds"
norm = matplotlib.colors.Normalize(vmin=np.min(feh_values),vmax=np.max(feh_values))

f, (a0, a1) = plt.subplots(nrows=2, ncols=1, gridspec_kw={'height_ratios': [2, 1]}, sharex=True)

a0.axvspan(6000, 7250, color='y', alpha=0.5, lw=0,zorder=0) # RRLs in instability strip (Catelan 2015)
a1.axvspan(6000, 7250, color='y', alpha=0.5, lw=0,zorder=0)
a0.plot(df_poststack["Teff"],df_poststack["Teff"],zorder=1,linestyle="--",color="k")
a1.plot([np.min(df_poststack["Teff"]),np.max(df_poststack["Teff"])],[0,0],zorder=1,linestyle="--",color="k")

a0.scatter(df_poststack["Teff"],
            df_poststack["Teff_retrieved_Balmer"],
            c=df_poststack["FeH"],
            cmap=colormap, norm=norm, edgecolor="k",zorder=2)

a1.scatter(df_poststack["Teff"],
            np.subtract(df_poststack["Teff_retrieved_Balmer_metalrich"],df_poststack["Teff"]),
            c=df_poststack["FeH"],
            cmap=colormap, norm=norm, edgecolor="k",zorder=2)


# annotation to check the color mapping
#for t in range(0,len(df_poststack["FeH"])):
#    plt.annotate(str(df_poststack["FeH"][t]), (df_poststack["Teff"][t],df_poststack["Teff_retrieved_Balmer"][t]))

# kludge to add legend while mapping colors correctly
for i in range(0,len(feh_values)):
    # indices reversed to get the order descending in the legend
    a0.scatter([0], [0], cmap=colormap, norm=norm, c=feh_values[-i-1],
                edgecolor="k", label="[Fe/H]="+str(feh_values[-i-1]))
    print(feh_values[i])

a0.set_ylabel("Retrieved T$_{eff}$")
a1.set_xlabel("Injected T$_{eff}$")
a1.set_ylabel("Retrieved T$_{eff}$ - Injected T$_{eff}$\n(based on trend for [Fe/H] $\geq$ -2.5)")

f.canvas.draw() # need before legend to render

a0.set_xlim([5500,8000])
a0.set_ylim([5500,8500])

a0.legend(loc="lower right")

plt.show()

print("USE NOTEBOOK VERSION OF THIS! OTHERWISE THE LEGEND DOESN'T HAVE RIGHT HANDLES!")
#plt.savefig("junk.pdf")
#import ipdb; ipdb.set_trace()
f.savefig("junk.pdf")
'''
