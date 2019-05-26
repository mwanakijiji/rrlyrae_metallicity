import glob
import pandas as pd
import os
import numpy as np
import sys
import matplotlib.pyplot as plt
from rrlyrae_metallicity.modules2 import *

class find_feh():
    '''
    Sample the a, b, c, d posteriors to find Fe/H given the equivalent widths
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
                                      usecols = [1,2,3,4],
                                      names = ["a", "b", "c", "d"],
                                      delim_whitespace = True)
        self.test = "test"

        
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
        from the MCMC chain, and sample 
        2. 
        '''

        # calculate Fe/H given EW samples from the MCMC chain
        print(self.test)
        print(type(self.mcmc_chain))

        # loop over samples in the MCMC chain
        for t in range(0,10):

            # stand-in EWs
            ersatz_H = 3. # angstroms
            ersatz_errH = 0.2
            ersatz_K = 7.5
            ersatz_errK = 0.3
            
            this_feh = self.feh_layden(a = self.mcmc_chain["a"][t],
                                  b = self.mcmc_chain["b"][t],
                                  c = self.mcmc_chain["c"][t],
                                  d = self.mcmc_chain["d"][t],
                                  H = ersatz_H,
                                  K = ersatz_K)

            ######################################
            
            # (pasted below from error_propagation*)
            # take one basis set Fe/H (integral over a Gaussian) and find what the mapped value should be 
            N = 100 # number of samples to take within the Gaussian error around Layden's Fe/H value
            gaussian_spread_errH = ersatz_errH ## ## change this in future
            gaussian_spread_errK = ersatz_errK ## ## change this in future
            
            layden_feh = feh_test # this is the discrete value
            feh_mapped_array = np.nan*np.ones((len(m_array),N)) # N_m_samples x N_Layden_samples

            # loop over each sample within the Gaussian around the EW errors
            for integal_piece in range(0,N):
    
                # set the offset (note mu=0; this is a relative offset)
                offset = np.random.normal(loc = 0.0, scale = gaussian_spread)
    
                # loop over all (m,b) combinations found further above
                for sample_num in range(0,len(m_array)):
    
                    feh_mapped_1sample = m_array[sample_num]*layden_feh*(1. + offset) + b_array[sample_num]
                    feh_mapped_array[sample_num][integal_piece] = feh_mapped_1sample

            #########################################
            
            print(self.mcmc_chain["a"][t])
            print(self.mcmc_chain["b"][t])
            print(self.mcmc_chain["c"][t])
            print(self.mcmc_chain["d"][t])
            print(this_feh)
            print('----------')
