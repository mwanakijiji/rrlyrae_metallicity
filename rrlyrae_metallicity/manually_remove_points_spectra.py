#!/usr/bin/env python
# coding: utf-8

# For enabling user to define where points should be masked out

# Created 2021 June 11 by E.S.

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob

# read list of names of files to examine
stem = "/Users/bandari/Documents/git.repos/rrlyrae_metallicity/rrlyrae_metallicity/"
list_file_name = "handleman_all_cleaned.txt"
list_file_df = pd.read_csv(stem + list_file_name, names=["file_name"])

'''
not working!!
# enable clicking on plot
def onclick(event):

    global ix, iy
    ix, iy = event.xdata, event.ydata
    print("x = "+str(ix)+", y = "+str(iy))

    global coords
    coords.append((ix, iy))

    if len(coords) == 2:
        fig.canvas.mpl_disconnect(cid)

    return coords
'''


# loop over each file

# read in the file
df_spec = pd.read_csv(stem + "sdss_spectra_cosmic_ray_removal/01a_all_normalized_once/" +
                    list_file_df["file_name"].iloc[10] + "_000", names=["wavel","flux"], delim_whitespace=True)

# display interactive plot of spectrum
plt.clf()
#fig = plt.figure() #figsize=(20,10))
#ax = fig.add_subplot(111)
plt.plot([3900,5000],[1,1],"k--")
plt.plot(df_spec["wavel"],df_spec["flux"])
plt.show()
#plt.plot()

# define bad regions
'''
break_statement = 0
while not break_statement:
    # let user loop over regions that are bad
    fig.canvas.draw()
    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    print(cid)
    break_statement = input("Done finding bad regions? [N]: ") or 0
'''

# wavelength ranges to mask
#mask_range = input("Input wavelength ranges to mask [[,],[,],...[,]]:")
mask_ranges_low = [int(i) for i in input("Input wavelength range lower limits to mask: ").split(" ")]
mask_ranges_high = [int(i) for i in input("Input wavelength range upper limits to mask: ").split(" ")]

#initialize mask
flagged_interactive = pd.DataFrame(df_spec["wavel"].copy())
flagged_interactive["flux_flag_1"] = False

# loop over continuous masked regions
for t in range(0,len(mask_ranges_low)):

    region_span = np.logical_and(df_spec["wavel"]>mask_ranges_low[t],df_spec["wavel"]<mask_ranges_high[t])
    plt.plot(df_spec["wavel"].where(region_span),df_spec["flux"].where(region_span),"k--")

    # mask out those elements (True: 'bad'; False: 'good; do not mask')
    flagged_interactive["flux_flag_1"].loc[region_span] = True

plt.show()

# write out mask as found interactively
csv_write_name = stem + "sdss_spectra_cosmic_ray_removal/01f_masks/" + "mask1_" + os.path.basename(matching[p])
# write to file (mode=x to avoid overwrite)
flagged_empirical.to_csv(csv_write_name, columns = ["wavel","flux_flag_1"], mode='x')
print("Wrote mask of spectrum as found interactively: " + csv_write_name)


# In[ ]:
