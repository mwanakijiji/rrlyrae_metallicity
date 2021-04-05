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
                good_ew_info_file = config_apply["data_dirs"]["DIR_EW_PRODS"]+config_apply["file_names"]["RESTACKED_EW_DATA_GOOD_ONLY"],
                 posteriors_subdir = config_apply["data_dirs"]["DIR_ABCD_POSTERIORS"],
                 posteriors_filename = config_apply["file_names"]["ABCD_POSTERIORS_FILE_NAME"],
                 verbose = False):

        self.good_ew_info_file = good_ew_info_file
        self.posteriors_subdir = posteriors_subdir
        self.posteriors_filename = posteriors_filename

        # the file containing the MCMC posterior chains of a,b,c,d
        self.ew_data = pd.read_csv(self.good_ew_info_file)

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
        # number of samples to take within the Gaussian errors around Balmer, CaIIK EWs
        N_EW_samples = 1

        # loop over samples in the MCMC chain
        N_MCMC_samples = len(self.mcmc_chain)

        # loop over the rows of the table of good EW data, with each row
        # corresponding to a spectrum
        for row_num in range(0,len(self.ew_data)):

            logging.info("Finding Fe/H for spectrum " + self.ew_data.iloc[row_num]["realization_spec_file_name"])

            Balmer_EW = self.ew_data.iloc[row_num]["EW_Balmer"]
            CaIIK_EW = self.ew_data.iloc[row_num]["EW_CaIIK"]
            err_Balmer_EW = self.ew_data.iloc[row_num]["EW_Balmer"]
            err_CaIIK_EW = self.ew_data.iloc[row_num]["EW_CaIIK"]

            # initialize array
            feh_sample_array = np.nan*np.ones((N_MCMC_samples, N_EW_samples))
            for t in range(0,N_MCMC_samples):

                # loop over each sample within the Gaussian around the EW errors
                for integral_piece in range(0,N_EW_samples):

                    ## ## NOTE USE OF ROBO ERROR AS GAUSSIAN WIDTH LEADS TO RIDICULOUSLY WIDE DISTRIBUTIONS
                    # set the offset (note mu=0; this is a relative offset)
                    offset_H = 0 # np.random.normal(loc = 0.0, scale = err_Balmer_EW)
                    offset_K = 0 # np.random.normal(loc = 0.0, scale = err_CaIIK_EW)

                    # find one value of Fe/H given those samples in Balmer and CaIIK EWs
                    feh_1sample = self.feh_layden(a = self.mcmc_chain["a"][t],
                                      b = self.mcmc_chain["b"][t],
                                      c = self.mcmc_chain["c"][t],
                                      d = self.mcmc_chain["d"][t],
                                      H = Balmer_EW + offset_H,
                                      K = CaIIK_EW + offset_K)

                    feh_sample_array[t][integral_piece] = feh_1sample

                print("MCMC sample")
                print(t)


            import matplotlib.pyplot as plt
            plt.hist(np.ravel(feh_sample_array), bins=100)
            plt.title(self.ew_data.iloc[row_num]["realization_spec_file_name"] + \
                    "\nmedian [Fe/H]=" +str(np.median(feh_sample_array)) + \
                    "\n+-" + str(np.std(feh_sample_array)))
            plt.savefig(self.ew_data.iloc[row_num]["realization_spec_file_name"] + ".pdf")
            plt.clf()
            print("median fe/h")
            print(np.median(feh_sample_array))
            print("std fe/h")
            print(np.std(feh_sample_array))
            #import ipdb; ipdb.set_trace()
