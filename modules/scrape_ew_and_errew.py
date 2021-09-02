'''
Scrape Robospect output and do some processing of the results
'''

import os
import sys
import glob
import logging
import pandas as pd
import numpy as np
import matplotlib
from astropy.io.fits import getdata
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from . import *

class Scraper():
    '''
    Scrape all the equivalent width info from the Robospect *.fits.robolines files
    '''

    def __init__(self,
                 subdir=config_red["data_dirs"]["DIR_ROBO_OUTPUT"],
                 file_scraped_info=config_red["data_dirs"]["DIR_EW_PRODS"]+config_red["file_names"]["SCRAPED_EW_ALL_DATA"],
                 verbose=False):

        # directory containing the *.fits.robolines
        # files with the EW info
        self.stem = '.' ## ##
        # subdirectory containing the *.c.dat files
        self.subdir = subdir ## ##

        # get list of filenames without the path
        ## ## note the string being sought here is specific to RW's synthetic spectra; this is a weakness here and needs to be fixed later!
        file_list_long = glob.glob(self.subdir+'/'+'*robolines')
        file_list_unsorted = [os.path.basename(x) for x in file_list_long]
        self.file_list = sorted(file_list_unsorted)

        # EW info will get scraped into this
        self.write_out_filename = file_scraped_info

        # return tables of EW data?
        self.verbose = verbose

    def __call__(self):

        def line_order_check(line_centers):
            '''
            Sanity check: are the lines listed in order?
            N.b. This checks the wavelengths using the given line list
            values (and not the fitted centers)
            '''

            logging.info('Verifying line centers...')
            logging.info(line_centers[0])
            glitch_count = int(0) # boolean for bookeeping
            if ((line_centers[0] < 3933.660-10) or
                (line_centers[0] > 3933.660+10)): # CaIIK
                logging.warning('CaIIK line center does not match!')
                glitch_count = int(1) # boolean for bookeeping
            if ((line_centers[1] < 3970.075-10) or
                  (line_centers[1] > 3970.075+10)): # H-epsilon (close to CaIIH)
                logging.warning('H-epsilon center (close to CaIIH) line does not match!')
                glitch_count = int(1) # boolean for bookeeping
            if ((line_centers[2] < 4101.7100-10) or
                  (line_centers[2] > 4101.7100+10)): # H-delta
                logging.warning('H-delta line center does not match!')
                glitch_count = int(1) # boolean for bookeeping
            if ((line_centers[3] < 4340.472-10) or
                  (line_centers[3] > 4340.472+10)): # H-gamma
                logging.warning('H-gamma line center does not match!')
                glitch_count = int(1) # boolean for bookeeping
            if ((line_centers[4] < 4861.290-10) or
                  (line_centers[4] > 4861.290+10)): # H-beta
                logging.warning('H-beta line center does not match!')
                glitch_count = 1 # boolean for bookeeping
            if (glitch_count == int(0)):
                logging.info('CaIIK, H-eps, H-del, H-gam, h_beta line centers are consistent')
            return

        df_master = pd.DataFrame() # initialize

        #print('FILE LIST')
        #print(self.file_list)
        #print(self.file_list[116])

        # loop over all filenames of realizations of empirical spectra, extract line data
        #import ipdb; ipdb.set_trace()
        for t in range(0, len(self.file_list)):

            # read in Robospect output
            logging.info("--------------------")
            logging.info("Reading in Robospect output from directory")
            logging.info(self.subdir)

            '''
            The following parses lines from Robospect *robolines output files,
            which look like the following, as of the v0.76 tag of Robospect:

            ## Units
            ##AA      [ AA           AA             None]        [ AA             AA                None]             [ AA           AA          None]       mAA        mAA             None      None     None       None
            ## Headers
            ##wave_0  [ gaussianMu   gaussianSigma  gaussianAmp] [ uncertaintyMu  uncertaintySigma  uncertaintyAmp]   [ priorMu      priorSigma  priorAmp]   EQW        uncertaintyEQW  chiSqr    flags    blendGroup comment
            3933.6600 [ 3933.618556  1.636451       -0.338310]   [ 0.043767       0.045441          0.008054]         [ 3934.427147  1.754001    0.384793]   1.387738   0.127230        0.004045  0x10020  0          CaII-K
            3970.0750 [ 3969.912002  6.497202       -0.626854]   [ 0.245555       0.237816          0.023196]         [ 3971.262223  4.535872    0.781687]   10.208984  1.331932        0.117392  0x10020  0          H-eps
            4101.7100 [ 4101.728498  6.829899       -0.596311]   [ 0.335244       0.327236          0.025288]         [ 4102.885050  4.878668    0.734648]   10.208852  1.637334        0.220112  0x10020  0          H-del
            4340.4720 [ 4340.374387  7.365172       -0.557777]   [ 0.395447       0.378434          0.025443]         [ 4340.943149  4.961159    0.689719]   10.297539  1.773505        0.300238  0x10020  0          H-gam
            4861.2900 [ 4861.316520  7.570797       -0.505060]   [ 0.441626       0.426212          0.025690]         [ 4861.746895  4.898021    0.635582]   9.584604   1.822847        0.377350  0x10020  0          H-beta
            '''

            df = pd.read_csv(self.subdir+'/'+self.file_list[t],
                             skiprows=19,
                             delim_whitespace=True,
                             index_col=False,
                             usecols=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],
                             names=  ["wavel_stated_center","[1","wavel_found_center","gaussianSigma","gaussianAmp",
                                        "[2","uncertaintyMu","uncertaintySigma","uncertaintyAmp",
                                        "[3","priorMu","priorSigma","priorAmp","EQW","uncertaintyEQW",
                                        "chiSqr","flags","blendGroup","line_name"])
            # remove dummy columns
            df = df.drop(columns=["[1","[2","[3"])
            # remove Robospect delimiter strings from columns and cast contents as floats
            print("df with new parsing")
            print(df)
            #import ipdb; ipdb.set_trace()
            print("file")
            print(self.subdir+'/'+self.file_list[t])
            logging.info("Parsing " + self.file_list[t])
            try:
                # this will fail if there are infs in the EWs
                df["gaussianAmp"] = df["gaussianAmp"].str.replace("]","")
                df["gaussianAmp"] = df["gaussianAmp"].astype(np.float)
                df["uncertaintyAmp"] = df["uncertaintyAmp"].str.replace("]","")
                df["uncertaintyAmp"] = df["uncertaintyAmp"].astype(np.float)
                df["priorAmp"] = df["priorAmp"].str.replace("]","")
                df["priorAmp"] = df["priorAmp"].astype(np.float)
            except:
                # skip this file
                logging.error("Parsing error! " + self.file_list[t])
                continue

            # check lines are in the right order
            # if they are not, a warning is printed in the log
            line_order_check(df['wavel_found_center'])

            # add two cols on the left: the filename, and the name of the line
            #s_length = len(df['mean']) # number of lines (should be 5)

            # file names
            df['robolines_file_name'] = pd.Series(self.file_list[t],
                                        index=df.index)

            # names of empirical spectra realizations (multiple ones
            # correspond to one empirical spectrum)
            # remove .robolines extension
            df['realization_spec_file_name'] = pd.Series(self.file_list[t].split(".robolines")[0],
                                              index=df.index)

            # names of original spectra
            df['original_spec_file_name'] = pd.Series(self.file_list[t].split(".robolines")[0].split("_")[-2],
                                              index=df.index)

            #df['star_name'] = pd.Series(self.file_list[t].split("__")[0], index=df.index)

            # names of the absorption lines
            df['line_name'] = ['CaIIK', 'Heps', 'Hdel', 'Hgam', 'Hbet']

            # print progress
            logging.info('Out of '+str(len(self.file_list))+' files, '+str(t+1)+' scraped...')

            # if this is the first list, start a master copy from it to concatenate stuff to it
            if (t == 0):
                df_master = df.copy()
            else:
                df_master = pd.concat([df_master, df])
                del df # clear variable

        # write to csv, while resetting the indices
        # note THIS TABLE INCLUDES ALL DATA, GOOD AND BAD
        #df_master_reset = df_master.reset_index(drop=True).copy()
        # this is effectively the same, but gets written out
        df_master.reset_index(drop=True).to_csv(self.write_out_filename)
        logging.info("Table of ALL EW info written to " + str(self.write_out_filename))
        #if self.verbose:
        #    return df_master_reset, df_master_reset_drop_bad_spectra
        return

def quality_check(
    read_in_filename = config_red["data_dirs"]["DIR_EW_PRODS"]+config_red["file_names"]["SCRAPED_EW_ALL_DATA"],
    write_out_filename = config_red["data_dirs"]["DIR_EW_PRODS"]+config_red["file_names"]["SCRAPED_EW_DATA_GOOD_ONLY"]):
    '''
    This reads in all the scraped EW data in raw form, removes spectra that have fits
    which are bad based on multiple criteria, and writes out another data_table

    INPUTS:
    read_in_filename: file name of the table with ALL scraped data from Robospect
    write_out_filename: file name of the table with spectra with any bad line fits removed
    '''

    # read in data
    all_data = pd.read_csv(read_in_filename)

    # make new column for 'good' (G) or 'bad' (B) based on the below criteria
    # (initialize all as 'G')
    all_data["quality"] = "G"

    # impose criteria for pruning of data

    # Criterion 1. Remove all rows with a Robospect flag ending with something other than zero
    # (i.e., Robospect found the fit to be bad)
    # make an array consisting of the last character in each spectrum's flag
    red_flag_array = ([u[-1] for u in all_data["flags"]])
    # consider bad flags to be of any flag with a nonzero last character
    where_red_flag = np.where(np.array(red_flag_array) != '0')
    # identify the synthetic spectrum names which have at least one line with a bad fit
    bad_robo_spectra = all_data["realization_spec_file_name"][np.squeeze(where_red_flag)]
    # remove duplicate names
    bad_robo_spectra_uniq = bad_robo_spectra.drop_duplicates()
    # flag as bad the spectra with those names
    all_data["quality"][all_data["realization_spec_file_name"].isin(bad_robo_spectra_uniq)] = "B"

    # Criterion 2. Remove rows where the line centers are not right, using steps similar to above
    # (specifically, if measured line center is more than 10 A away from perfect center)
    where_bad_line_center = np.where(np.abs(np.subtract(all_data["wavel_found_center"],all_data["wavel_stated_center"]) > 10))
    bad_line_center_spectra = all_data["realization_spec_file_name"][np.squeeze(where_bad_line_center,axis=0)] # squeeze necessary to preserve finite size
    bad_line_center_spectra_uniq = bad_line_center_spectra.drop_duplicates()
    all_data["quality"][all_data["realization_spec_file_name"].isin(bad_line_center_spectra_uniq)] = "B"

    # Criterion 3. Remove rows with EWs which are clearly unrealistically large which slipped through other checks
    # (this is particularly an issue with the CaIIK line, which is close to CaIIH)
    # set cutoff at 18 A, based on inspection of >200 Robospect plots of fits to
    # synthetic spectra; all those with CaIIK EW > 18 are clearly not fit right,
    # and all those with EW < 18 look acceptable -E.S.
    where_bad_CaIIK = np.where(np.logical_and(all_data["line_name"] == "CaIIK", all_data["EQW"] > 18))
    bad_CaIIK_spectra = all_data["realization_spec_file_name"][np.squeeze(where_bad_CaIIK)]
    bad_CaIIK_spectra_uniq = bad_CaIIK_spectra.drop_duplicates()
    all_data["quality"][all_data["realization_spec_file_name"].isin(bad_CaIIK_spectra_uniq)] = "B"

    # Criterion 4. Remove bad phases (for empirical data)
    '''
    min_good, max_good = phase_regions()
    #[...]
    #good_indices = np.intersect1d(good_phase, good_metal)
    #[...]
    '''

    # Last step: Write only the rows with a good ("G") flag to file
    # (note that if AT LEAST one absorption line was found to be bad, ALL the
    # data corresponding to that spectrum is removed)
    pruned_data = all_data[all_data.quality == "G"]#.reset_index()
    pruned_data.to_csv(write_out_filename)

    logging.info("--------------------------")
    logging.info('Scraped Robospect output written to')
    logging.info(write_out_filename)

    return pruned_data


def stack_spectra(
    read_in_filename = config_red["data_dirs"]["DIR_EW_PRODS"]+config_red["file_names"]["SCRAPED_EW_DATA_GOOD_ONLY"],
    write_out_filename = config_red["data_dirs"]["DIR_EW_PRODS"]+config_red["file_names"]["RESTACKED_EW_DATA_GOOD_ONLY"],
    fits_dir = config_red["data_dirs"]["DIR_SYNTH_SPEC"],
    objective = "find_abcd"):
    '''
    Takes output of quality_check() and
    1.  calculates errors in EW, using two methods
    2.  generates a scaled, net Balmer line EW
    3.  transposes and stacks data so that the data has *rows* of spectra and *cols* of absorption lines

    INPUTS:
    read_in_filename: file name of scraped Robospect data, after removing bad spectra
    write_out_filename: name of file to contain re-stacked data
    fits_dir: directory of FITS files which we will use for obtaining info
        from the header (Fe/H, etc.)
    objective: if it is to find the calibration, dig up intermediary FITS files
        to get supplementary information on Fe/H etc. from the headers; if it is
        to apply the calibration, just stack the scraped EW information
        ["find_abcd"/"apply_abcd"]
    '''

    def error_scatter_ew(df_pass):
        '''
        Adds a column of errors, as calculated using the method of taking the
        scatter in measured EWs of different lines (as opposed to taking Robospec's
        errors at face value)
        '''

        # get list of original file names with no repeats
        orig_file_array = np.array((df_pass["original_spec_file_name"].drop_duplicates()))

        # add new columns of nans
        df_pass["err_EW_Hbeta_from_EW_variation"] = np.nan
        df_pass["err_EW_Hgamma_from_EW_variation"] = np.nan
        df_pass["err_EW_Hdelta_from_EW_variation"] = np.nan
        df_pass["err_EW_Heps_from_EW_variation"] = np.nan
        df_pass["err_EW_CaIIK_from_EW_variation"] = np.nan

        for orig_file_name_num in range(0,len(orig_file_array)):

            # mask all rows that do not correspond to the original spectrum

            this_orig_spec = orig_file_array[orig_file_name_num]
            df_masked = df_pass.where(df_pass["original_spec_file_name"] == this_orig_spec)

            # find stdev of EWs, as measured for all realizations of those file names
            '''
            df_masked["err_EW_Hbeta_from_EW_variation"] = np.nanstd(df_masked["EW_Hbeta"])
            df_masked["err_EW_Hgamma_from_EW_variation"] = np.nanstd(df_masked["EW_Hgamma"])
            df_masked["err_EW_Hdelta_from_EW_variation"] = np.nanstd(df_masked["EW_Hdelta"])
            df_masked["err_EW_Heps_from_EW_variation"] = np.nanstd(df_masked["EW_Heps"])
            '''

            # insert into columns of input table
            try:
                idx = df_pass.index[df_pass["original_spec_file_name"] == this_orig_spec] # indices
                df_pass.loc[idx, "err_EW_Hbeta_from_EW_variation"] = np.nanstd(df_masked["EW_Hbeta"])
                df_pass.loc[idx, "err_EW_Hgamma_from_EW_variation"] = np.nanstd(df_masked["EW_Hgamma"])
                df_pass.loc[idx, "err_EW_Hdelta_from_EW_variation"] = np.nanstd(df_masked["EW_Hdelta"])
                df_pass.loc[idx, "err_EW_Heps_from_EW_variation"] = np.nanstd(df_masked["EW_Heps"])
                df_pass.loc[idx, "err_EW_CaIIK_from_EW_variation"] = np.nanstd(df_masked["EW_CaIIK"])
            except:
                print("Anomaly in finding scatter in EW measurements in " + str(this_orig_spec))

        return df_pass

    # read in data
    df_prestack = pd.read_csv(read_in_filename)

    # make list of individual spectra for which we have EW data, and
    # initialize DataFrame to hold the re-cast data

    list_indiv_spectra = list(df_prestack["realization_spec_file_name"].drop_duplicates())

    num_indiv_spectra = len(list_indiv_spectra)

    df_poststack = pd.DataFrame(columns=["realization_spec_file_name",
                                         "original_spec_file_name",
                                         "FeH", "err_FeH",
                                         "logg", "alpha","Teff",
                                         "EW_Hbeta", "err_EW_Hbeta_from_robo",
                                         "EW_Hdelta", "err_EW_Hdelta_from_robo",
                                         "EW_Hgamma", "err_EW_Hgamma_from_robo",
                                         "EW_Heps", "err_EW_Heps_from_robo",
                                         "EW_CaIIK", "err_EW_CaIIK_from_robo"], index=range(num_indiv_spectra)) # initialize

    for t in range(0,num_indiv_spectra):
        # loop over all spectra realizations we have measured EWs from to populate the dataframe

        this_spectrum = list_indiv_spectra[t]
        logging.info("Extracting EW data corresponding to " + this_spectrum)

        if (objective == "find_abcd"):
            # read in intermediary FITS file to extract values from header
            logging.info("Reading in intermediary FITS file " + this_spectrum.split(".")[0] + ".fits for header data")
            image, hdr = getdata(fits_dir + this_spectrum.split(".")[0] + ".fits", header=True, ignore_missing_end=True)

            logg = hdr["LOGG"]
            teff = hdr["TEFF"]
            alpha = hdr["ALPHA"]
            feh = hdr["FEH"]
            err_feh = 0.15 # ersatz for now

        # select data from table relevant to this spectrum
        data_this_spectrum = df_prestack.where(df_prestack["realization_spec_file_name"] == this_spectrum).dropna().reset_index()


        # calculate line errors using
        # method 1: values directly from Robospect

        # extract original file name (the one from which realizations are made)
        orig_name = data_this_spectrum["original_spec_file_name"].drop_duplicates().values[0]

        try:
            # extract Balmer lines from the table of data from all the spectra
            Hbeta = data_this_spectrum["EQW"].where(data_this_spectrum["line_name"] == "Hbet").dropna().values[0]
            err_Hbeta = data_this_spectrum["uncertaintyEQW"].where(data_this_spectrum["line_name"] == "Hbet").dropna().values[0]

            Hgamma = data_this_spectrum["EQW"].where(data_this_spectrum["line_name"] == "Hgam").dropna().values[0]
            err_Hgamma = data_this_spectrum["uncertaintyEQW"].where(data_this_spectrum["line_name"] == "Hgam").dropna().values[0]

            Hdelta = data_this_spectrum["EQW"].where(data_this_spectrum["line_name"] == "Hdel").dropna().values[0]
            err_Hdelta = data_this_spectrum["uncertaintyEQW"].where(data_this_spectrum["line_name"] == "Hdel").dropna().values[0]

            Heps = data_this_spectrum["EQW"].where(data_this_spectrum["line_name"] == "Heps").dropna().values[0]
            err_Heps = data_this_spectrum["uncertaintyEQW"].where(data_this_spectrum["line_name"] == "Heps").dropna().values[0]

            CaIIK = data_this_spectrum["EQW"].where(data_this_spectrum["line_name"] == "CaIIK").dropna().values[0]
            err_CaIIK = data_this_spectrum["uncertaintyEQW"].where(data_this_spectrum["line_name"] == "CaIIK").dropna().values[0]

            # fill in that row in the dataframe
            df_poststack.iloc[t]["realization_spec_file_name"] = this_spectrum
            df_poststack.iloc[t]["original_spec_file_name"] = orig_name
            df_poststack.iloc[t]["EW_Hbeta"] = Hbeta
            df_poststack.iloc[t]["err_EW_Hbeta_from_robo"] = err_Hbeta
            df_poststack.iloc[t]["EW_Hdelta"] = Hdelta
            df_poststack.iloc[t]["err_EW_Hdelta_from_robo"] = err_Hdelta
            df_poststack.iloc[t]["EW_Hgamma"] = Hgamma
            df_poststack.iloc[t]["err_EW_Hgamma_from_robo"] = err_Hgamma
            df_poststack.iloc[t]["EW_Heps"] = Heps
            df_poststack.iloc[t]["err_EW_Heps_from_robo"] = err_Heps
            df_poststack.iloc[t]["EW_CaIIK"] = CaIIK
            df_poststack.iloc[t]["err_EW_CaIIK_from_robo"] = err_CaIIK

        except:
            print("Data anomaly; skipping " + this_spectrum)

    # save intermediary table of data, before adding rescaled Balmer line
    intermediary_file_name = write_out_filename + "_intermediary"
    logging.info("Writing out intermediary file of stacked Robospect EWs and rescaled Balmer lines to " + intermediary_file_name)
    df_poststack.to_csv(intermediary_file_name)

    # calculate line errors using
    # method 2: stdev of line EWs
    df_poststack = error_scatter_ew(df_poststack)

    if (objective == "find_abcd"):
        # if the stellar spectra are synthetic, add in that info
        df_poststack.iloc[t]["logg"] = logg
        df_poststack.iloc[t]["Teff"] = teff
        df_poststack.iloc[t]["alpha"] = alpha
        df_poststack.iloc[t]["FeH"] = feh
        df_poststack.iloc[t]["err_FeH"] = err_feh

    # to generate a net Balmer line, make a rescaling of Hgamma
    # based on Hdelta
    logging.info("Making a net Balmer line")

    # fit a straight line to all the Hgam vs Hdel
    x_data = df_poststack["EW_Hdelta"].values.astype(float) # Hdel
    y_data = df_poststack["EW_Hgamma"].values.astype(float) # Hgam
    # better names for clarity below
    Hgamma_all = np.copy(y_data)
    err_Hgamma_all = df_poststack["err_EW_Hgamma_from_robo"].values

    # safety check that both pairs of coordinates are simultaneously finite
    idx_good = np.logical_and(np.isfinite(x_data),np.isfinite(y_data))
    coeff, cov = np.polyfit(x_data[idx_good], y_data[idx_good], 1, full=False, cov=True)
    m = coeff[0]
    b = coeff[1]
    err_m = np.sqrt(np.diag(cov))[0]
    err_b = np.sqrt(np.diag(cov))[1]

    print("m:")
    print(m)
    print("b:")
    print(b)

    # generate a rescaled Hgamma, call it rHgam
    rHgam_all = np.divide(np.subtract(Hgamma_all, b), m)
    # find corresponding error`
    piece1 = np.add(np.power(err_Hgamma_all,2),np.power(err_b,2)).astype(float)
    piece2 = np.power(np.subtract(Hgamma_all,b),2)
    piece3 = np.divide(np.power(err_m,2),np.power(m,2))
    err_rHgam_all = np.multiply(rHgam_all,np.sqrt(np.subtract(np.divide(piece1,piece2),piece3)))

    # add column of rescaled Hgamma to DataFrame
    df_poststack["EW_Balmer"] = rHgam_all
    df_poststack["err_EW_Balmer"] = err_rHgam_all

    # FYI plot: how do Balmer lines scale with each other?
    '''
    plt.clf()
    plt.title("Scaling of lines with Hdelta")
    plt.scatter(df_poststack["EW_Hdelta"],df_poststack["EW_Hbeta"], s=3, label="Hbeta")
    plt.scatter(df_poststack["EW_Hdelta"],np.add(df_poststack["EW_Hgamma"],4), s=3, label="Hgamma+4")
    plt.scatter(df_poststack["EW_Hdelta"],np.add(df_poststack["EW_Heps"],8), s=3, label="Heps+8")
    #plt.ylim([0,15])
    plt.xlabel("EW, Hdelta (Angstr)")
    plt.ylabel("EW, non-Hdelta (Angstr)")
    plt.legend()
    plt.savefig("junk_balmer_rescalings.pdf")

    # FYI plot: KH plot
    plt.clf()
    plt.title("KH plot")
    plt.errorbar(df_poststack["EW_resc_Hgamma"],df_poststack["EW_CaIIK"],
                 yerr=df_poststack["err_EW_CaIIK_from_robo"],
                 marker="o", markersize=2, mfc="k", mec="k", ecolor="gray", linestyle="")
    plt.ylim([0,30])
    plt.xlabel("EW, net Balmer (Angstr)")
    plt.ylabel("EW, CaIIK (Angstr)")
    plt.savefig("junk_KH_plot.pdf")
    '''

    logging.info("------------------------------")
    #logging.info("Data will be written out to file " + write_out_filename)
    #input("Hit [Enter] to continue")
    logging.info("Writing out re-casting of Robospect EWs and rescaled Balmer line to " + write_out_filename)
    df_poststack.to_csv(write_out_filename)

    return
