import glob
import pandas as pd
import os
import numpy as np
import sys
import matplotlib.pyplot as plt
from rrlyrae_metallicity.modules2 import *
from . import *

class find_feh():
    '''
    Sample the a, b, c, d posteriors as found via the MCMC, to find Fe/H given the equivalent widths
    of the science spectra
    '''

    def __init__(self,
                 posteriors_subdir = config_apply["data_dirs"]["DIR_ABCD_POSTERIORS"],
                 posteriors_filename = config_apply["file_names"]["ABCD_POSTERIORS_FILE_NAME"],
                 verbose = False):

        self.posteriors_subdir = posteriors_subdir
        self.posteriors_filename = posteriors_filename

        # the file containing the MCMC posterior chains of a,b,c,d
        #self.mcmc_chain_file = self.posteriors_subdir + self.posteriors_filename

        # read in the chain data
        self.mcmc_chain = pd.read_csv(self.posteriors_subdir + self.posteriors_filename,
                                      usecols = [1, 2, 3, 4],
                                      names = ["a", "b", "c", "d"],
                                      delim_whitespace = True)
        #self.test = "test"


    def __call__(self):

        pass


    def feh_layden(self,a,b,c,d,H,K):
        '''
        Finds Fe/H given equivalent widths (in angstroms), from
        K = a + b*H + c*[Fe/H] + d*H*[Fe/H]  (Layden 1994 Eqn. 7)
        '''

        feh = np.divide(K-a-b*H,c+d*H)

        return feh


    def sample_feh(self):
        '''
        Find a Fe/H value for a combination of (a,b,c,d)
        from the MCMC chain, and sample from the Balmer and
        CaIIK EWs, given their errors
        '''

        ## find/input EWs for a single spectrum here; use stand-in EWs for the moment
        ersatz_H = 3. # angstroms
        ersatz_errH = 0.2
        ersatz_K = 7.5
        ersatz_errK = 0.3

        gaussian_spread_errH = ersatz_errH ## ## change this in future
        gaussian_spread_errK = ersatz_errK ## ## change this in future

        # number of samples to take within the Gaussian errors around Balmer, CaIIK EWs
        N_EW_samples = 100

        # loop over samples in the MCMC chain
        N_MCMC_samples = 10 ## change this later
        # initialize array
        feh_sample_array = np.nan*np.ones((N_MCMC_samples, N_EW_samples))
        for t in range(0,N_MCMC_samples):

            # loop over each sample within the Gaussian around the EW errors
            for integral_piece in range(0,N_EW_samples):

                # set the offset (note mu=0; this is a relative offset)
                offset_H = np.random.normal(loc = 0.0, scale = gaussian_spread_errH)
                offset_K = np.random.normal(loc = 0.0, scale = gaussian_spread_errK)

                # find one value of Fe/H given those samples in Balmer and CaIIK EWs
                feh_1sample = self.feh_layden(a = self.mcmc_chain["a"][t],
                                  b = self.mcmc_chain["b"][t],
                                  c = self.mcmc_chain["c"][t],
                                  d = self.mcmc_chain["d"][t],
                                  H = ersatz_H + offset_H,
                                  K = ersatz_K + offset_K)

                feh_sample_array[t][integral_piece] = feh_1sample

            print("MCMC sample")
            print(t)

        print("median fe/h")
        print(np.median(feh_sample_array))
        print("std fe/h")
        print(np.std(feh_sample_array))
