'''
This is an emcee wrapper for fitting the Layden '94 metallicity
calibration to equivalent widths of RR Lyrae spectra
'''

import sys
import time
import numpy as np
import pandas as pd
import emcee
import corner
from rrlyrae_metallicity.modules2 import *


def rrmetal(h_pass, f_pass, p_pass):
    '''
    Layden+ 1994 model

    INPUTS:
    h_pass: Balmer EW (angstroms)
    f_pass: [Fe/H]
    p_pass: coefficients [a,b,c,d]

    OUTPUTS:
    k_pass: CaIIK EW (angstroms)
    '''
    k_pass = (p_pass[0] +
              p_pass[1]*h_pass +
              p_pass[2]*f_pass +
              p_pass[3]*h_pass*f_pass)

    return k_pass


def chi_sqd_fcn(xi_pass,
                yi_pass,
                zi_pass,
                sig_xi_pass,
                sig_yi_pass,
                sig_zi_pass,
                a_pass,
                b_pass,
                c_pass,
                d_pass):
    '''
    Chi-squared

    INPUTS:
    xi_pass: Balmer EW (angstroms)
    yi_pass: [Fe/H]
    zi_pass: CaIIK EW (angstroms)

    OUTPUTS:
    val: chi^2
    '''

    # IDL syntax
    # numerator_i = (zi_pass-a_pass-b_pass*xi_pass-c_pass*yi_pass-d_pass*xi_pass*yi_pass)**2
    # denominator_i = (sig_xi_pass**2)* ((b_pass+d_pass*yi_pass)**2) +
    #                  (sig_yi_pass**2)*((c_pass+d_pass*xi_pass)^2)+ sig_zi_pass**2

    # ... and the awkward Python syntax
    base = np.subtract(np.subtract(np.subtract(np.subtract(zi_pass,a_pass),np.multiply(b_pass,xi_pass)),\
                                 np.multiply(c_pass,yi_pass)),np.multiply(np.multiply(d_pass,xi_pass),yi_pass))
    numerator_i = base**2
    term_i = (sig_xi_pass**2)
    term_ii = (np.add(b_pass, np.multiply(d_pass, yi_pass))**2)
    term_iii = (sig_yi_pass**2)
    term_iv = (np.add(c_pass, np.multiply(d_pass, xi_pass)))**2
    term_v = (sig_zi_pass**2)
    denominator_i = np.add(np.add(np.multiply(term_i, term_ii),
                                  np.multiply(term_iii, term_iv)), term_v)

    i_element = np.divide(numerator_i, denominator_i)
    val = np.sum(i_element)

    return val


def lnprob(walker_pos,
           Teff_pass,
           measured_H_pass,
           measured_F_pass,
           measured_K_pass,
           err_measured_H_pass,
           err_measured_F_pass,
           err_measured_K_pass):
    '''
    Nat log of probability density

    OUTPUTS:
    ln(prior*like)
    '''

    # walker_pos is the proposed walker position in N-D (likely 4-D) space
    # (i.e., these are the inputs to the model)
    lp = lnprior(walker_pos) # prior
    if not np.isfinite(lp): # afoul of prior
      return -np.inf
    result = -np.divide(1, 2*Teff_pass)*chi_sqd_fcn(measured_H_pass,
                                                    measured_F_pass,
                                                    measured_K_pass,
                                                    err_measured_H_pass,
                                                    err_measured_F_pass,
                                                    err_measured_K_pass,
                                                    walker_pos[0],
                                                    walker_pos[1],
                                                    walker_pos[2],
                                                    walker_pos[3])
    return lp + result # ln(prior*like)


def lnprior(theta):
    '''
    Prior

    INPUTS:
    theta: array of parameter values

    OUTPUTS: 0 or -inf (top-hat priors only)
    '''
    a_test, b_test, c_test, d_test = theta

    # top-hat priors
    if ((np.abs(a_test) < 40) and
        (np.abs(b_test) < 5) and
        (np.abs(c_test) < 20) and
        (np.abs(d_test) < 10)):
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
                 scraped_ew_source_dir=config["data_dirs"]["DIR_BIN"],
                 mcmc_text_output_dir=config["data_dirs"]["DIR_BIN"],
                 corner_plot_putput_dir=config["data_dirs"]["DIR_BIN"]):

        # name of file with final K, H, FeH, and error values (and not the others from the noise-churned spectra)
        self.scraped_ew_filename = (scraped_ew_source_dir +
                                    config["file_names"]["KH_WINNOWED"])

        # name of file of the MCMC output
        self.mcmc_text_output = mcmc_text_output_dir + config["file_names"]["MCMC_OUTPUT"]

        # name of corner plot of the MCMC output
        self.corner_file_string = corner_plot_putput_dir + config["file_names"]["MCMC_CORNER"]

        # read in boundaries of good phase regions
        self.min_good, self.max_good = phase_regions()

    def __call__(self):

        # read in EWs, Fe/Hs, phases, errors, etc.
        print("--------------------------")
        print('Reading in data ...')
        print(self.scraped_ew_filename)

        ## ## make df_choice.Spectrum -> df_choice["Spectrum etc.
        df_choice = pd.read_csv(self.scraped_ew_filename,delim_whitespace=False)

        #THIS IS THE ORIGINAL, SINCE EWS WERE IN MILLIANG
        # EWs in table are in angstroms and are mislabeled as mA (2020 Jan 12)
        name = df_choice['original_spec_file_name']
        #caii = np.divide(df_choice['K'], 1000.)
        caii = df_choice['K']
        #ecaii = np.divide(df_choice['err_K'], 1000.)
        ecaii = df_choice['err_K']
        #ave = np.divide(df_choice['balmer'], 1000.)
        ave = df_choice['balmer']
        eave = df_choice['err_balmer']
        #eave = np.divide(df_choice['err_balmer'], 1000.)
        ## ## THE BELOW FEH VALUES NEED TO BE CHECKED/FIXED
        feh = df_choice['final_feh_center']
        efeh = np.subtract(df_choice['final_feh_center'],
                           df_choice['final_feh_lower'])
        #import ipdb; ipdb.set_trace()


        phase = df_choice['phase']
        #period = df_choice.type
        #star_type = dataFloats[:, 15]

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

        # val from previous IDL runs (kind of deprecated; just
        # appears as a constant in the MCMC)
        Teff = 0.0586758

        # coefficients from first line of Table 8 in Layden+ 1994
        # (reddening not included), to serve as MCMC starting point
        a_layden = 13.542
        b_layden = -1.072
        c_layden = 3.971
        d_layden = -0.271
        sigma_a_layden = 0.416
        sigma_b_layden = 0.076
        sigma_c_layden = 0.285
        sigma_d_layden = 0.052

        # starting position, before adding a perturbation
        param_array_0_Layden = [float(a_layden),
                               float(b_layden),
                               float(c_layden),
                               float(d_layden)]
        sigmas_0_Layden = [float(sigma_a_layden),
                           float(sigma_b_layden),
                           float(sigma_c_layden),
                           float(sigma_d_layden)]

        # remove high metallicity stars
        ## PUT INTO CONFIG FILE
        metal_upper_limit = 1.0

        # impose conditions using anonymous functions
        good_phase = find_indices(phase,
                                  lambda q: (q < self.max_good) & (q > self.min_good))
        good_metal = find_indices(feh,
                                  lambda r: (r < metal_upper_limit))
        # return common indices
        good_indices = np.intersect1d(good_phase, good_metal)

        g_ave = ave[good_indices]
        g_eave = eave[good_indices]
        g_feh = feh[good_indices]
        g_efeh = efeh[good_indices]
        g_caii = caii[good_indices]
        g_ecaii = ecaii[good_indices]


        ################# MCMC setup #################

        print("--------------------------")
        print('Setting up MCMC ...')

        ndim = len(param_array_0_Layden) # dimensions of space to explore
        nwalkers = 8 # number of chains

        # convert the one starting point into a nwalkers*ndim array with gaussian-offset starting points
        p0 = [np.add(param_array_0_Layden,
                     np.multiply(param_array_0_Layden, 1e-4*np.random.randn(ndim))) for i in range(nwalkers)]

        # set up sampler
        sampler = emcee.EnsembleSampler(nwalkers,
                                        ndim,
                                        lnprob,
                                        args=[Teff, g_ave, g_feh, g_caii, g_eave, g_efeh, g_ecaii])

        # burn-in
        burnIn = 1e3 # 1e5 seems to be necessary for the slow-converging parameter 'd'
        posAfterBurn, prob, state = sampler.run_mcmc(p0, burnIn)

        # post-burn-in
        start_time = time.time()
        post_burn_in_links = 3e3

        ################# SAVE PROGRESSIVELY TO TEXT FILE #################
        ## ## refer to these code snippets from Foreman-Mackey's website
        # IMPORTANT: sampler will only have memory of the last iteration if
        # storechain flag is set to False

        print("--------------------------")
        print("Saving MCMC chains to text file ...")

        # post-burn-in calculate and save iteratively
        f = open(self.mcmc_text_output, "w")
        f.close()
        progBarWidth = 30
        start_time = time.time()
        for i, result in enumerate(sampler.sample(posAfterBurn,
                                                  iterations=post_burn_in_links,
                                                  storechain=True)):
            position = result[0]
            f = open(self.mcmc_text_output, "a") # append
            for k in range(position.shape[0]): # loop over number of chains
                position_string = str(position[k]).strip("[]") # convert to string
                f.write("{0:4d} {1:s}\n".format(k, " ".join(str(p) for p in position[k])))
            n = int((progBarWidth+1) * float(i) / post_burn_in_links) # update progress bar
            sys.stdout.write("\r[{0}{1}]".format("#" * n, " " * (progBarWidth - n)))
            f.close()
        elapsed_time = time.time() - start_time
        sys.stdout.write(" Done!\n")
        sys.stdout.write("{0:s} {1:10d} {2:s}\n".format("Elapsed time: ",
                                                        int(elapsed_time), "sec"))
        print("--------------------------")
        sys.stdout.write("MCMC chain data written out to\n")
        sys.stdout.write(str(self.mcmc_text_output))

        # corner plot (requires 'storechain=True' in enumerate above)
        samples = sampler.chain[:, int(burnIn):, :].reshape((-1, ndim))
        fig = corner.corner(samples, labels=["$a$", "$b$", "$c$", "$d$"],
                            quantiles=[0.16, 0.5, 0.84],
                            title_fmt='.2f',
                            show_titles=True,
                            verbose=True,
                            title_kwargs={"fontsize": 12})
        fig.savefig(self.corner_file_string)
        print("--------------------------")
        print("Corner plot of MCMC posteriors written out to")
        print(str(self.corner_file_string))

        # if its necessary to read in MCMC output again
        #data = np.loadtxt(self.mcmc_text_output, usecols=range(1,5))

        # This code snippet from Foreman-Mackey's emcee documentation, v2.2.1 of
        # https://emcee.readthedocs.io/en/stable/user/line.html#results
        a_mcmc, b_mcmc, c_mcmc, d_mcmc = map(lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                                             zip(*np.percentile(samples, [16, 50, 84], axis=0)))

        print("--------------------------")
        print("Coefficients a, b, c, d, and errors (see corner plot):")
        print(a_mcmc, '\n', b_mcmc, '\n', c_mcmc, '\n', d_mcmc)

        print("--------------------------")
        print("MCMC data written to ")
        print(self.mcmc_text_output)
