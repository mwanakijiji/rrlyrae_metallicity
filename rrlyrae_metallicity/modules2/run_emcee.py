import numpy as np
import emcee
import matplotlib.pyplot as plt
import corner
import time
from os.path import dirname, realpath, sep, pardir
import sys
import pandas as pd
from rrlyrae_metallicity.modules2 import *

'''
This is an emcee wrapper for fitting the Layden '94 metallicity calibration to equivalent widths of
RR Lyrae spectra
'''

def rrmetal(hPass,fPass,pPass):
    '''
    Layden+ 1994 model

    INPUTS:
    hPass: Balmer EW (angstroms)
    fPass: [Fe/H]
    pPass: coefficients [a,b,c,d]

    OUTPUTS:
    kPass: CaIIK EW (angstroms)
    '''
    kPass = pPass[0] + pPass[1]*hPass + pPass[2]*fPass + pPass[3]*hPass*fPass
    return kPass


def chi_sqd_fcn(xiPass, yiPass, ziPass, sig_xiPass, sig_yiPass, sig_ziPass, aPass, bPass, cPass, dPass):
    '''
    Chi-squared

    INPUTS:
    xiPass: Balmer EW (angstroms)
    yiPass: [Fe/H]
    ziPass: CaIIK EW (angstroms)

    OUTPUTS:
    val: chi^2
    '''

    # IDL syntax
    # numerator_i = (ziPass-aPass-bPass*xiPass-cPass*yiPass-dPass*xiPass*yiPass)**2
    # denominator_i = (sig_xiPass**2)* ((bPass+dPass*yiPass)**2)+ (sig_yiPass**2)* ((cPass+dPass*xiPass)^2)+ sig_ziPass**2

    # ... and the awkward Python syntax
    base = np.subtract(np.subtract(np.subtract(np.subtract(ziPass,aPass),np.multiply(bPass,xiPass)),\
                                 np.multiply(cPass,yiPass)),np.multiply(np.multiply(dPass,xiPass),yiPass))
    numerator_i = base**2
    termI = (sig_xiPass**2)
    termII = (np.add(bPass,np.multiply(dPass,yiPass))**2)
    termIII = (sig_yiPass**2)
    termIV = (np.add(cPass,np.multiply(dPass,xiPass)))**2
    termV = (sig_ziPass**2)
    denominator_i = np.add(np.add(np.multiply(termI,termII),np.multiply(termIII,termIV)),termV)
  
    i_element = np.divide(numerator_i,denominator_i)
    val = np.sum(i_element)

    return val


def lnprob(walkerPos, TeffPass, measured_HPass, measured_FPass, measured_KPass, err_measured_HPass, err_measured_FPass, \
           err_measured_KPass):
    '''
    Nat log of probability density

    OUTPUTS:
    ln(prior*like)
    '''

    # walkerPos is the proposed walker position in N-D (likely 4-D) space (i.e., these are the inputs to the model)
    lp = lnprior(walkerPos) # prior
    if not np.isfinite(lp): # afoul of prior
      return -np.inf
    result = -np.divide(1,2*TeffPass)*chi_sqd_fcn(measured_HPass, measured_FPass, measured_KPass, \
                                                err_measured_HPass, err_measured_FPass, err_measured_KPass, \
                                                walkerPos[0], walkerPos[1], walkerPos[2], walkerPos[3])
    return lp + result # ln(prior*like)


def lnprior(theta):
    '''
    Prior

    INPUTS:
    theta: array of parameter values

    OUTPUTS: 0 or -inf (top-hat priors only)
    '''
    aTest, bTest, cTest, dTest = theta
    if (np.abs(aTest) < 40) and (np.abs(bTest) < 5) and (np.abs(cTest) < 20) and (np.abs(dTest) < 10): # top-hat priors
        return 0.0
    return -np.inf


def find_indices(lst, condition):
    '''
    Stand-in equivalent of IDL 'where' function

    INPUTS:
    lst: list
    condition: an anonymous function

    RETURNS:
    
    '''
    
    return [i for i, elem in enumerate(lst) if condition(elem)]


class RunEmcee():
    '''
    Run the emcee MCMC to obtain coefficients a, b, c, d
    '''
    
    def __init__(self,
                 scraped_ew_source_dir = config["data_dirs"]["DIR_BIN"],
                 mcmc_text_output_dir = config["data_dirs"]["DIR_BIN"],
                 corner_plot_putput_dir = config["data_dirs"]["DIR_BIN"]):

        # name of file with final K, H, FeH, and error values (and not the others from the noise-churned spectra)
        self.scraped_ew_filename = scraped_ew_source_dir + config["file_names"]["KH_WINNOWED_FILE_NAME"]

        # name of file of the MCMC output
        self.mcmc_text_output = mcmc_text_output_dir + config["file_names"]["MCMC_OUTPUT"]

        # name of corner plot of the MCMC output
        self.cornerFileString = corner_plot_putput_dir + config["file_names"]["MCMC_CORNER"]

        # read in boundaries of good phase regions
        self.min_good, self.max_good = phase_regions()
        
    def __call__(self):

        # read in EWs, Fe/Hs, phases, errors, etc.
        print("--------------------------")
        print('Reading in data ...')
        print(self.scraped_ew_filename)
        dfChoice = pd.read_csv(self.scraped_ew_filename,
                               delim_whitespace = False) ## ## rename dfChoice and make dfChoice.Spectrum -> dfChoice["Spectrum etc.

        name = dfChoice['empir_spec_name']
        caii = np.divide(dfChoice['K'],1000.) # EWs in table are in milliangstroms
        ecaii = np.divide(dfChoice['err_K'],1000.)
        ave = np.divide(dfChoice['balmer'],1000.)
        eave = np.divide(dfChoice['err_balmer'],1000.)

        ## ## THE BELOW FEH VALUES NEED TO BE CHECKED/FIXED
        feh = dfChoice['final_feh_center']
        efeh = np.subtract(dfChoice['final_feh_center'],
                           dfChoice['final_feh_lower'])
        
        phase = dfChoice['phase']
        #period = dfChoice.type
        #star_type = dataFloats[:,15]

        print("name")
        print(name)
        print("caii")
        print(caii)
        print("ecaii")
        print(ecaii)
        print("ave")
        print(ave)
        print("eave")
        print(eave)
        print("feh")
        print(feh)
        print("efeh")
        print(efeh)
        print("phase")
        print(phase)
        
        # fix some values
        Teff = 0.0586758 # from previous IDL runs (kind of deprecated; just appears as a constant in the MCMC)
        
        # coefficients from first line of Table 8 in Layden+ 1994 (reddening not included), to serve as MCMC starting point
        a_layden = 13.542
        b_layden = -1.072
        c_layden = 3.971
        d_layden = -0.271  
        sigma_a_layden = 0.416
        sigma_b_layden = 0.076
        sigma_c_layden = 0.285
        sigma_d_layden = 0.052

        # starting position, before adding a perturbation
        paramArray_0_Layden = [float(a_layden),
                               float(b_layden),
                               float(c_layden),
                               float(d_layden)]
        sigmas_0_Layden = [float(sigma_a_layden),
                           float(sigma_b_layden),
                           float(sigma_c_layden),
                           float(sigma_d_layden)]

        # remove high metallicity stars
        ## PUT INTO CONFIG FILE  
        metalUpperLimit = 1.0

        # impose conditions using anonymous functions
        good_phase = find_indices(phase, lambda q: (q < self.max_good) & (q > self.min_good))
        good_metal = find_indices(feh, lambda r: (r < metalUpperLimit))
        good_indices = np.intersect1d(good_phase, good_metal) # return common index values

        g_ave = ave[good_indices]
        g_eave = eave[good_indices]
        g_feh = feh[good_indices]
        g_efeh = efeh[good_indices]
        g_caii = caii[good_indices]
        g_ecaii = ecaii[good_indices]


        ################# MCMC setup #################

        print("--------------------------")
        print('Setting up MCMC ...')

        ndim = len(paramArray_0_Layden) # dimensions of space to explore
        nwalkers = 8 # number of chains

        # convert the one starting point into a nwalkers*ndim array with gaussian-offset starting points
        p0 = [np.add(paramArray_0_Layden,
                     np.multiply(paramArray_0_Layden,1e-4*np.random.randn(ndim))) for i in range(nwalkers)]

        # set up sampler
        sampler = emcee.EnsembleSampler(nwalkers,
                                        ndim,
                                        lnprob,
                                        args=[Teff, g_ave, g_feh, g_caii, g_eave, g_efeh, g_ecaii])

        # burn-in
        burnIn = 1e3 # 1e5 seems to be necessary for the slow-converging parameter 'd'
        posAfterBurn, prob, state = sampler.run_mcmc(p0, burnIn)

        # post-burn-in
        startTime = time.time()
        postBurnInLinks = 3e3

        ################# SAVE PROGRESSIVELY TO TEXT FILE #################
        ## ## refer to these code snippets from Foreman-Mackey's website
        # IMPORTANT: sampler will only have memory of the last iteration if storechain flag is set to False

        print("--------------------------")
        print("Saving MCMC chains to text file ...")

        # post-burn-in calculate and save iteratively
        f = open(self.mcmc_text_output, "w")
        f.close()
        progBarWidth = 30
        start_time = time.time()
        for i, result in enumerate(sampler.sample(posAfterBurn, iterations=postBurnInLinks, storechain=True)):
            position = result[0]
            f = open(self.mcmc_text_output, "a") # append
            for k in range(position.shape[0]): # loop over number of chains
                positionString = str(position[k]).strip("[]") # convert to string
                f.write("{0:4d} {1:s}\n".format(k, " ".join(str(p) for p in position[k])))
            n = int((progBarWidth+1) * float(i) / postBurnInLinks) # update progress bar
            sys.stdout.write("\r[{0}{1}]".format("#" * n, " " * (progBarWidth - n)))
            f.close()
        elapsed_time = time.time() - start_time
        sys.stdout.write(" Done!\n")
        sys.stdout.write("{0:s} {1:10d} {2:s}\n".format("Elapsed time: ", int(elapsed_time), "sec"))
        print("--------------------------")
        sys.stdout.write("MCMC chain data written out to")
        sys.stdout.write(str(self.mcmc_text_output))

        # corner plot (requires 'storechain=True' in enumerate above)
        samples = sampler.chain[:, int(burnIn):, :].reshape((-1, ndim))
        fig = corner.corner(samples, labels=["$a$", "$b$", "$c$", "$d$"],
                            quantiles=[0.16, 0.5, 0.84],
                            title_fmt='.2f',
                            show_titles=True, verbose=True, title_kwargs={"fontsize": 12})
        fig.savefig(self.cornerFileString)
        print("--------------------------")
        print("Corner plot of MCMC posteriors written out to")
        print(str(self.cornerFileString))

        # if its necessary to read in MCMC output again
        #data = np.loadtxt(self.mcmc_text_output, usecols=range(1,5))

        # This code snippet from Foreman-Mackey's emcee documentation, v2.2.1 of
        # https://emcee.readthedocs.io/en/stable/user/line.html#results
        a_mcmc, b_mcmc, c_mcmc, d_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                                             zip(*np.percentile(samples, [16, 50, 84], axis=0)))

        print("--------------------------")
        print("Coefficients a, b, c, d, and errors (see corner plot):")
        print(a_mcmc,'\n',b_mcmc,'\n',c_mcmc,'\n',d_mcmc)

        print("--------------------------")
        print("MCMC data written to ")
        print(self.mcmc_text_output)
