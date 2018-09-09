import numpy as np
import emcee
import matplotlib.pyplot as plt
import corner
import time
from os.path import dirname, realpath, sep, pardir
import sys
import pandas as pd


# This is an emcee wrapper for fitting the Layden '94 metallicity calibration to equivalent widths of
# RR Lyrae spectra

# THEN RUN EMCEE WITH FINAL EWS, ERR_EWS, FEHS, ERR_FEHS

class run_emcee():
    
    ##############################################################################
    # STEP 5: RUN EMCEE ON THE SPACE, GET VALUES FOR a, b, c, d (applicable only to A)
    ##############################################################################
    
    def __init__(self):
        print('')
        
    def __call__(self):
        print('')




        
# fcn: Layden+ 1994 model (just for fitting end result)
def rrmetal(hPass,fPass,pPass):
    kPass = pPass[0] + pPass[1]*hPass + pPass[2]*fPass + pPass[3]*hPass*fPass
    return kPass

# fcn: chi-squared
def chi_sqd_fcn(xiPass, yiPass, ziPass, sig_xiPass, sig_yiPass, sig_ziPass, aPass, bPass, cPass, dPass):
    # xiPass: measured H
    # yiPass: measured [Fe/H]
    # ziPass: measured K

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

# fcn: nat log of probability density
def lnprob(walkerPos, TeffPass, measured_HPass, measured_FPass, measured_KPass, err_measured_HPass, err_measured_FPass, \
           err_measured_KPass):
    # walkerPos is the proposed walker position in N-D (likely 4-D) space (i.e., these are the inputs to the model)
    lp = lnprior(walkerPos) # prior
    if not np.isfinite(lp): # afoul of prior
      return -np.inf
    result = -np.divide(1,2*TeffPass)*chi_sqd_fcn(measured_HPass, measured_FPass, measured_KPass, \
                                                err_measured_HPass, err_measured_FPass, err_measured_KPass, \
                                                walkerPos[0], walkerPos[1], walkerPos[2], walkerPos[3])
    return lp + result # ln(prior*like)

# fcn: prior
def lnprior(theta): # theta is an array of parameter values
    aTest, bTest, cTest, dTest = theta
    if (np.abs(aTest) < 20) and (np.abs(bTest) < 5) and (np.abs(cTest) < 10) and (np.abs(dTest) < 2): # top-hat priors
        return 0.0
    return -np.inf

# fcn: equivalent of IDL 'where' function
def find_indices(lst, condition): # condition will be in form of an anonymous function
    return [i for i, elem in enumerate(lst) if condition(elem)]

##########################################################
## END FUNCTION DEFINITIONS
##########################################################

# read in EWs, errors, etc.
#name = np.loadtxt("EWfin.dat.csv", usecols=range(1,2), delimiter=",", skiprows=1, dtype=string)
# dataFloats = np.loadtxt("EWfin.dat.csv", usecols=range(1,15), delimiter=",", skiprows=1) # deprecated line
print('Reading in data ...')

allSpec_1stPass_allSNR = pd.read_csv('ron_data/EWAall.lab', delim_whitespace=True)
allSpec_1stPass_highSNR = pd.read_csv('ron_data/EWASNR25.lab', delim_whitespace=True)
allSpec_2ndPass_allSNR = pd.read_csv('ron_data/EWfin.dat', delim_whitespace=True)
allSpec_2ndPass_highSNR = pd.read_csv('ron_data/EWfinSNR25.dat', delim_whitespace=True)

# the below line was a past way to draw a subset of lines out of one of the other files
# allSpec_1stPass_highSNR = allSpec_1stPass_allSNR.where(allSpec_1stPass_allSNR['Spectrum'].isin(allSpec_2ndPass_highSNR['Spectrum']))



## USER INPUTS HERE
filenameString = "rrl_test_20chains_3e5links_1stPass_highSNR_ABC.dat"
dfChoice = allSpec_1stPass_highSNR.copy() # choice which data table to use
## END USER INPUTS

name = dfChoice.Spectrum
caii = dfChoice.EW_K
ecaii = dfChoice.eEW_K
hdel = dfChoice.EW_D
ehdel = dfChoice.eEW_D
hgam = dfChoice.EW_G
ehgam = dfChoice.eEW_G
hbet = dfChoice.EW_B
ehbet = dfChoice.eEW_B
ave = dfChoice.Have
eave = dfChoice.eHave
feh = dfChoice.Fe_H
efeh = dfChoice.eFe_H
phase = dfChoice.phase
period = dfChoice.type
#star_type = dataFloats[:,15]

# fix some values
Teff = 0.0586758 # from previous IDL runs (kind of deprecated; just appears as a constant in the MCMC)
# coefficients from first line of Table 8 in Layden+ 1994 (reddening not included)
a_layden = 13.542
b_layden = -1.072
c_layden = 3.971
d_layden = -0.271  
sigma_a_layden = 0.416
sigma_b_layden = 0.076
sigma_c_layden = 0.285
sigma_d_layden = 0.052

paramArray_0_Layden = [float(a_layden),float(b_layden),float(c_layden),float(d_layden)] # starting position, before adding a perturbation
sigmas_0_Layden = [float(sigma_a_layden),float(sigma_b_layden),float(sigma_c_layden),float(sigma_d_layden)]

# remove high metallicity stars and bad phases  
phaseUpperLimit = 0.8
phaseLowerLimit = 0.05
metalUpperLimit = 1.0

# impose conditions using anonymous functions
good_phase = find_indices(phase, lambda q: (q<phaseUpperLimit)&(q>phaseLowerLimit))
good_metal = find_indices(feh, lambda r: (r<metalUpperLimit))
good_indices = np.intersect1d(good_phase, good_metal) # return common values

g_ave = ave[good_indices]
g_eave = eave[good_indices]
g_feh = feh[good_indices]
g_efeh = efeh[good_indices]
g_caii = caii[good_indices]
g_ecaii = ecaii[good_indices]


################# MCMC setup #################

print('Setting up MCMC ...')

ndim = len(paramArray_0_Layden) # dimensions of space to explore
nwalkers = 20 # number of chains

# convert the one starting point into a nwalkers*ndim array with gaussian-offset starting points
p0 = [np.add(paramArray_0_Layden,np.multiply(paramArray_0_Layden,1e-4*np.random.randn(ndim))) for i in range(nwalkers)]

# set up sampler
sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=[Teff, g_ave, g_feh, g_caii, \
                                                              g_eave, g_efeh, g_ecaii])

# burn-in
burnIn = 1e5 # 1e5 seems to be necessary for the slow-converging parameter 'd'
posAfterBurn, prob, state = sampler.run_mcmc(p0, burnIn)

# post-burn-in
startTime = time.time()
postBurnInLinks = 3e5

################# SAVE PROGRESSIVELY TO TEXT FILE #################
## ## refer to this code snipped from Foreman-Mackey's website
# IMPORTANT: sampler will only have memory of the last iteration if storechain flag is set to False

print('Saving MCMC chains to text file ...')

# post-burn-in calculate and save iteratively
f = open(filenameString, "w")
f.close()
progBarWidth = 30
start_time = time.time()
for i, result in enumerate(sampler.sample(posAfterBurn, iterations=postBurnInLinks, storechain=True)):
    position = result[0]
    f = open(filenameString, "a") # append
    for k in range(position.shape[0]): # loop over number of chains
        positionString = str(position[k]).strip("[]") # convert to string
        f.write("{0:4d} {1:s}\n".format(k, " ".join(str(p) for p in position[k])))
    n = int((progBarWidth+1) * float(i) / postBurnInLinks) # update progress bar
    sys.stdout.write("\r[{0}{1}]".format("#" * n, " " * (progBarWidth - n)))
    f.close()
elapsed_time = time.time() - start_time
sys.stdout.write(" Done!\n")
sys.stdout.write("{0:s} {1:10d} {2:s}\n".format("Elapsed time: ", int(elapsed_time), "sec"))

# corner plot (requires 'storechain=True' in enumerate above)
samples = sampler.chain[:, burnIn:, :].reshape((-1, ndim))
fig = corner.corner(samples, labels=["$a$", "$b$", "$c$", "$d$"])
fig.savefig("corner.png")


################# WITHOUT SAVING PROGRESSIVELY TO FILE #################
'''
sampler.run_mcmc(posAfterBurn, postBurnInLinks) # if not progressively saving to text file
print(time.time()-startTime)
# remove burn-in
samples = sampler.chain[:, burnIn:, :].reshape((-1, ndim))

# corner plot
fig = corner.corner(samples, labels=["$a$", "$b$", "$c$", "$d$"])
fig.savefig("triangle.png")

print("Mean acceptance fraction: {0:.3f}".format(np.mean(sampler.acceptance_fraction)))
'''

