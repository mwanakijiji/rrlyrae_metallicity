from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

def ca_corrxn(corrxn_map):
# READ IN CA ABSORPTION AND MAKE CORRECTIONS

    image, header = fits.getdata(corrxn_map,0,header=True)

    # CALCULATE FINAL FEH VALUES FOR OUR OWN STARS, AND WRITE OUT 



    
    ## ## CAUTION: TEST TO SEE IF THE CONTENT IN THE KEYS IS IN ORDER (I.E., MAKE A PLOT AND SEE IF ITS THE SAME IF DATASETS ARE OVERLAID INDIVIDUALLY)


