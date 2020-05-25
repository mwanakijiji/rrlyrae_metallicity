'''
Scrape Robospect output and do some processing of the results
'''

import os
import sys
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pl
from rrlyrae_metallicity.modules2 import *

class Scraper():
    '''
    Scrape all the equivalent width info from the Robospect *.fits.robolines files
    '''

    def __init__(self,
                 subdir=config["data_dirs"]["DIR_ROBO_OUTPUT"],
                 file_scraped_info=config["file_names"]["MCD_LARGE_BAD_REMOVED"],
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
        self.write_out_filename = subdir + file_scraped_info

        # return tables of EW data?
        self.verbose = verbose

    def __call__(self):

        def line_order_check(line_centers):
            '''
            Sanity check: are the lines listed in order?
            N.b. This checks the wavelengths using the given line list
            values (and not the fitted centers)
            '''
            #import ipdb; ipdb.set_trace()
            print(line_centers[0])
            if ((line_centers[0] < 3933.660-10) or
                (line_centers[0] > 3933.660+10)): # CaIIK
                print('Lines not matching!')
                sys.exit()  # ... and abort
            elif ((line_centers[1] < 3970.075-10) or
                  (line_centers[1] > 3970.075+10)): # H-epsilon (close to CaIIH)
                print('Lines not matching!')
                sys.exit()
            elif ((line_centers[2] < 4101.7100-10) or
                  (line_centers[2] > 4101.7100+10)): # H-delta
                print('Lines not matching!')
                sys.exit()
            elif ((line_centers[3] < 4340.472-10) or
                  (line_centers[3] > 4340.472+10)): # H-gamma
                print('Lines not matching!')
                sys.exit()
            elif ((line_centers[4] < 4861.290-10) or
                  (line_centers[4] > 4861.290+10)): # H-beta
                print('Lines not matching!')
                sys.exit()
            return

        df_master = pd.DataFrame() # initialize

        #print('FILE LIST')
        #print(self.file_list)
        #print(self.file_list[116])

        # loop over all filenames of realizations of empirical spectra, extract line data
        #import ipdb; ipdb.set_trace()
        for t in range(0, len(self.file_list)):

            # read in Robospect output
            '''
            df = pd.read_csv(self.subdir+'/'+self.file_list[t],
                             skiprows=19,
                             delim_whitespace=True,
                             index_col=False,
                             usecols=[    0,     2,   3,     6,  7,   11,  13,  14,   15,     16,         18],
                             names=  ["#x0","mean","3","flux","7","EQW","13","14","15","flags","line_name"])
            '''
            print("Reading in Robospect output from directory")
            print(self.subdir+'/')
            print("--------------------")
            df = pd.read_csv(self.subdir+'/'+self.file_list[t],
                             skiprows=19,
                             delim_whitespace=True,
                             index_col=False,
                             usecols=[0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18],
                             names=  ["#x0","[1","mean","gaussianSigma","gaussianAmp",
                                        "[2","uncertaintyMu","uncertaintySigma","uncertaintyAmp",
                                        "[3","priorMu","priorSigma","priorAmp","EQW","uncertaintyEQW",
                                        "chiSqr","flags","blendGroup","line_name"])
            ##names=["#x0","mean","sigma","flux","7","10","11","EQW","14","chi","flags","line_name"]
            ## old command here
            #import ipdb; ipdb.set_trace()
            ## ## TEST IF STATEMENT
            #import ipdb; ipdb.set_trace()
            #import ipdb; ipdb.set_trace()
            #if self.file_list[t] == "600020p02.smo_000.robolines":
            #    import ipdb; ipdb.set_trace()
            '''
            df = pd.read_csv(self.subdir+'/'+self.file_list[t],
                             header=16,
                             delim_whitespace=True,
                             index_col=False,
                             usecols=np.arange(17))
            '''
            #print(df['#x0'])
            #print(self.file_list[t])
            #import ipdb; ipdb.set_trace()
            # check lines are in the right order
            line_order_check(df['#x0'])

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
            print('Out of '+str(len(self.file_list))+' files, '+str(t+1)+' scraped...')

            # if this is the first list, start a master copy from it to concatenate stuff to it
            if (t == 0):
                df_master = df.copy()
            else:
                df_master = pd.concat([df_master, df])
                del df # clear variable

        # write to csv, while resetting the indices
        # note THIS TABLE INCLUDES ALL DATA, GOOD AND BAD
        df_master_reset = df_master.reset_index(drop=True).copy()
        # this is effectively the same, but gets written out
        df_master.reset_index(drop=True).to_csv(self.subdir+config["file_names"]["MCD_ALL_DATA"])

        ## IF WE ARE INTERESTED IN SPECTRA THAT HAVE ALL WELL-FIT LINES
        # remove all rows with a flag ending with something other than zero (i.e., the fit is bad)
        # make an array consisting of the last character in each spectrum's flag
        red_flag_array = ([u[-1] for u in df_master_reset["flags"]])
        # consider bad flags to be of any flag with a nonzero last character
        where_red_flag = np.where(np.array(red_flag_array) != '0')

        # identify the synthetic spectrum names which have at least one line with a bad fit
        bad_synth_spectra = df_master_reset['realization_spec_file_name'][np.squeeze(where_red_flag)]
        # remove duplicate names
        bad_synth_spectra_uniq = bad_synth_spectra.drop_duplicates()
        # keep only the spectra that have all lines well-fit
        df_master_reset_drop_bad_spectra = df_master_reset.where(
            ~df_master_reset['realization_spec_file_name'].isin(bad_synth_spectra_uniq))

        # write to csv
        # note THIS TABLE HAS SPECTRA WITH ANY BAD ROWS REMOVED
        df_master_reset_drop_bad_spectra.to_csv(self.write_out_filename)
        #import ipdb; ipdb.set_trace()
        print("--------------------------")
        print('Scraped Robospect output written to')
        print(self.write_out_filename)

        if self.verbose:
            return df_master_reset, df_master_reset_drop_bad_spectra


class findHK():
    '''
    Read in Robospect EWs of synthetic spectra, rescale them, average them, and plot KH space
    '''

    def __init__(self,
                 source_subdir=config["data_dirs"]["DIR_ROBO_OUTPUT"],
                 phase_subdir=config["data_dirs"]["DIR_SRC"],
                 hk_write_subdir=config["data_dirs"]["DIR_SRC"],
                 plot_write_subdir=config["data_dirs"]["DIR_FYI_INFO"]):

        self.source_subdir = source_subdir
        self.phase_subdir = phase_subdir
        self.scraped_ew_filename = self.source_subdir + \
          config["file_names"]["MCD_LARGE_BAD_REMOVED"]
        self.hk_file_name = hk_write_subdir + \
          config["file_names"]["MORE_REALISTIC"]

        # read in line data
        print(self.scraped_ew_filename)
        self.line_data = pd.read_csv(self.scraped_ew_filename,
                                     delim_whitespace=False)

        # initialize arrays: essential info
        self.original_spec_file_name_array = []
        self.star_name_array = []
        self.H_data_array = []
        self.K_data_array = []
        self.err_H_data_array = []
        self.err_K_data_array = []

        # initialize arrays: other info
        self.Hbet_data_array = [] # Hbeta
        self.err_Hbet_data_array = []
        self.Hgam_data_array = [] # Hgamma
        self.err_Hgam_data_array = []
        self.rHgam_data_array = [] # rescaled Hgamma
        self.err_rHgam_data_array = []
        self.Hdel_data_array = [] # Hdelta
        self.err_Hdel_data_array = []
        self.Heps_data_array = [] # Hepsilon
        self.err_Heps_data_array = []

        # read in boundaries of good phase regions
        # (these phase boundaries are not used for any
        # calculations at this stage; only plotting)
        self.min_good, self.max_good = phase_regions()

        # indicate subdirectory where FYI plot will be written
        self.plot_write_subdir = plot_write_subdir

    def __call__(self):

        # make a list of all UNIQUE, EMPIRICAL spectrum names
        unique_spec_names_pre_index_reset = self.line_data.drop_duplicates(
            subset='original_spec_file_name')['original_spec_file_name']

        # drop row of NaNs and smooth the indexing
        unique_spec_names = unique_spec_names_pre_index_reset.dropna(
            ).reset_index(drop=True)

        # fit a straight line to Hgam vs Hdel
        x_data = self.line_data['EQW'].where(
            self.line_data['line_name'] == 'Hdel'
            ).dropna() # Hdel
        y_data = self.line_data['EQW'].where(
            self.line_data['line_name'] == 'Hgam'
            ).dropna() # Hgam
        Hgam = np.copy(y_data)
        coeff, cov = np.polyfit(x_data, y_data, 1, full=False, cov=True)
        m = coeff[0]
        b = coeff[1]
        err_m = np.sqrt(np.diag(cov))[0]
        err_b = np.sqrt(np.diag(cov))[1]

        # generate a rescaled Hgam, call it rHgam
        rHgam_all = np.divide(np.subtract(Hgam, b), m)

        # prepare data for a plot
        # loop over every EMPIRICAL spectrum and assemble SYNTHETIC data into arrays
        for p in range(0, len(np.array(unique_spec_names.values))):

            print('unique')
            print(unique_spec_names)
            print("Synthetic data being assembled corresponding to spectrum "+\
                  np.array(unique_spec_names)[p])
            # extract all synthetic data corresponding to this empirical spectrum
            data_for_this_empir_spectrum = self.line_data.where(
                self.line_data['original_spec_file_name'][0:-4] == np.array(unique_spec_names)[p]
                )

            # scrape EWs
            raw_Hbet_data = data_for_this_empir_spectrum['EQW'].where(
                self.line_data['line_name'] == 'Hbet'
                )
            raw_Hgam_data = data_for_this_empir_spectrum['EQW'].where(
                self.line_data['line_name'] == 'Hgam'
                )
            raw_Hdel_data = data_for_this_empir_spectrum['EQW'].where(
                self.line_data['line_name'] == 'Hdel'
                )
            raw_Heps_data = data_for_this_empir_spectrum['EQW'].where(
                self.line_data['line_name'] == 'Heps'
                )
            raw_K_data = data_for_this_empir_spectrum['EQW'].where(
                self.line_data['line_name'] == 'CaIIK'
                )

            # scrape Robospect EW errors
            raw_Hbet_err_EW = data_for_this_empir_spectrum['uncertaintyEQW'].where(
                self.line_data['line_name'] == 'Hbet'
                )
            raw_Hgam_err_EW = data_for_this_empir_spectrum['uncertaintyEQW'].where(
                self.line_data['line_name'] == 'Hgam'
                )
            raw_Hdel_err_EW = data_for_this_empir_spectrum['uncertaintyEQW'].where(
                self.line_data['line_name'] == 'Hdel'
                )
            raw_Heps_err_EW = data_for_this_empir_spectrum['uncertaintyEQW'].where(
                self.line_data['line_name'] == 'Heps'
                )
            raw_K_err_EW = data_for_this_empir_spectrum['uncertaintyEQW'].where(
                self.line_data['line_name'] == 'CaIIK'
                )

            # rescale EWs
            Hbet_data_wnans = np.array(np.copy(raw_Hbet_data))
            Hgam_data_wnans = np.array(np.copy(raw_Hgam_data))
            Hdel_data_wnans = np.array(np.copy(raw_Hdel_data))
            Heps_data_wnans = np.array(np.copy(raw_Heps_data))
            K_data_wnans = np.array(np.copy(raw_K_data))
            # rescale Hgam EWs
            rHgam_data_wnans = np.array(np.divide(
                np.subtract(raw_Hgam_data, b), m))

            # rescale EW errors
            Hbet_err_EW_wnans = np.array(np.copy(raw_Hbet_err_EW))
            Hgam_err_EW_wnans = np.array(np.copy(raw_Hgam_err_EW))
            Hdel_err_EW_wnans = np.array(np.copy(raw_Hdel_err_EW))
            Heps_err_EW_wnans = np.array(np.copy(raw_Heps_err_EW))
            K_err_EW_wnans = np.array(np.copy(raw_K_err_EW))
            # rescale Hgam EW errors
            # (this is only approximate, as it assumes zero correlation between m, b)
            # rHgam = (Hgam-b)/m
            # err_rHgam = rHgam*[ ( (err_Hgam+err_b) / (Hgam - b) ) + (err_m/m) ]
            err_piece1 = np.divide(Hgam_err_EW_wnans+err_b,Hgam_data_wnans-b) + \
                        np.divide(err_m,m)
            rHgam_err_EW_wnans = np.multiply(rHgam_data_wnans,err_piece1)

            # remove nans from EWs
            Hbet_data = Hbet_data_wnans[np.isfinite(Hbet_data_wnans)]
            Hgam_data = Hgam_data_wnans[np.isfinite(Hgam_data_wnans)]
            Hdel_data = Hdel_data_wnans[np.isfinite(Hdel_data_wnans)]
            Heps_data = Heps_data_wnans[np.isfinite(Heps_data_wnans)]
            rHgam_data = rHgam_data_wnans[np.isfinite(rHgam_data_wnans)]
            K_data = K_data_wnans[np.isfinite(K_data_wnans)]

            # remove nans from err EWs
            Hbet_err_EW = Hbet_err_EW_wnans[np.isfinite(Hbet_err_EW_wnans)]
            Hgam_err_EW = Hgam_err_EW_wnans[np.isfinite(Hgam_err_EW_wnans)]
            Hdel_err_EW = Hdel_err_EW_wnans[np.isfinite(Hdel_err_EW_wnans)]
            Heps_err_EW = Heps_err_EW_wnans[np.isfinite(Heps_err_EW_wnans)]
            rHgam_err_EW = rHgam_err_EW_wnans[np.isfinite(rHgam_err_EW_wnans)]
            K_err_EW = K_err_EW_wnans[np.isfinite(K_err_EW_wnans)]

            # get the H-K synthetic data together to form individual
            # points in H,K space
            # ( note balmer EW = 0.5*(Hdel + rHgam) )
            balmer_data_allsynthetic_spec = np.nanmean([Hdel_data, rHgam_data], axis=0)
            Hdel_err_sqrd = np.power(Hdel_err_EW,2)
            rHgam_err_sqrd = np.power(rHgam_err_EW,2)
            if np.logical_or(len(np.ravel(Hdel_err_sqrd)) != 1,len(np.ravel(rHgam_err_sqrd)) != 1):
                # if these are not single numbers (that is, they are arrays of more than one member)
                stopgap = input("Something is wrong-- The arrays Hdel_err_sqrd" + \
                                " or rHgam_err_sqrd have length > 1; error bars will not be" + \
                                " calculated correctly. Or else redesign this part of the code.")
            balmer_err_EW_allsynthetic_spec = np.sqrt(Hdel_err_sqrd[0]+rHgam_err_sqrd[0])
            K_data_allsynthetic_spec = np.copy(K_data)
            import ipdb; ipdb.set_trace()
            # the actual points to plot (or record in a table)
            Hbet_data_pt = np.nanmedian(Hbet_data)
            Hgam_data_pt = np.nanmedian(Hgam_data)
            rHgam_data_pt = np.nanmedian(rHgam_data)
            Hdel_data_pt = np.nanmedian(Hdel_data)
            Heps_data_pt = np.nanmedian(Heps_data)
            balmer_data_pt = np.nanmedian(balmer_data_allsynthetic_spec)
            K_data_pt = np.nanmedian(K_data_allsynthetic_spec)

            # the error bars as found by Robospect
            err_Hbet_data = np.nanmedian(Hbet_err_EW_wnans)
            err_Hgam_data = np.nanmedian(Hgam_err_EW_wnans)
            err_Hdel_data = np.nanmedian(Hdel_err_EW_wnans)
            err_Heps_data = np.nanmedian(Heps_err_EW_wnans)
            err_rHgam_data = np.nanmedian(rHgam_err_EW_wnans)
            err_balmer_data = np.nanmedian(balmer_err_EW_allsynthetic_spec)
            err_K_data = np.nanmedian(K_err_EW_wnans)

            # the error bars
            '''
            THE BELOW IS FOR FINDING THE ERRORS BY USING NOISE-CHURNED SPECTRA,
            AND NOT THE ROBOSPECT ERRORS (WE FOUND THIS GIVES ERRORS TOO SMALL)
            print("Hbet std being calculated from " + str(len(Hbet_data)) + " noise-churned spectra.")
            err_Hbet_data = np.nanstd(Hbet_data)
            err_Hgam_data = np.nanstd(Hgam_data)
            err_rHgam_data = np.nanstd(rHgam_data)
            err_Hdel_data = np.nanstd(Hdel_data)
            err_Heps_data = np.nanstd(Heps_data)
            err_balmer_data = np.nanstd(balmer_data_allsynthetic_spec)
            err_K_data = np.nanstd(K_data_allsynthetic_spec)
            '''

            #plt.plot(balmer_data_pt, K_data_pt)
            #plt.errorbar(balmer_data_pt, K_data_pt, yerr=err_K_data, xerr=err_balmer_data)

            # append data to arrays: essential info
            self.original_spec_file_name_array = np.append(self.original_spec_file_name_array,
                                                   np.array(unique_spec_names)[p])
            self.star_name_array = np.append(self.star_name_array,
                                             str(np.array(unique_spec_names)[p]))
            self.H_data_array = np.append(self.H_data_array,
                                          balmer_data_pt)
            self.err_H_data_array = np.append(self.err_H_data_array,
                                              err_balmer_data)
            self.K_data_array = np.append(self.K_data_array,
                                          K_data_pt)
            self.err_K_data_array = np.append(self.err_K_data_array,
                                              err_K_data)

            # append data to arrays: other info
            self.Hbet_data_array = np.append(self.Hbet_data_array,
                                             Hbet_data_pt)
            self.err_Hbet_data_array = np.append(self.err_Hbet_data_array,
                                                 err_Hbet_data)
            self.Hgam_data_array = np.append(self.Hgam_data_array,
                                             Hgam_data_pt)
            self.err_Hgam_data_array = np.append(self.err_Hgam_data_array,
                                                 err_Hgam_data)
            self.rHgam_data_array = np.append(self.rHgam_data_array,
                                              rHgam_data_pt) # rescaled Hgamma
            self.err_rHgam_data_array = np.append(self.err_rHgam_data_array,
                                                  err_rHgam_data)
            self.Hdel_data_array = np.append(self.Hdel_data_array,
                                             Hdel_data_pt)
            self.err_Hdel_data_array = np.append(self.err_Hdel_data_array,
                                                 err_Hdel_data)
            self.Heps_data_array = np.append(self.Heps_data_array,
                                             Heps_data_pt)
            self.err_Heps_data_array = np.append(self.err_Heps_data_array,
                                                 err_Heps_data)

            # clear some variables
            balmer_data_allsynthetic_spec = None
            K_data_allsynthetic_spec = None
            balmer_data_allsynthetic_spec = None
            K_data_allsynthetic_spec = None

        # put everything into a dataframe
        d = {'original_spec_file_name': self.original_spec_file_name_array,
             'star_name': self.star_name_array,
             'Hbet': self.Hbet_data_array,
             'err_Hbet': self.err_Hbet_data_array,
             'Hgam': self.Hgam_data_array,
             'err_Hgam': self.err_Hgam_data_array,
             'Hdel': self.Hdel_data_array,
             'err_Hdel': self.err_Hdel_data_array,
             'Heps': self.Heps_data_array,
             'err_Heps': self.err_Heps_data_array,
             'rHgam': self.rHgam_data_array,
             'err_rHgam': self.err_rHgam_data_array,
             'balmer': self.H_data_array,
             'err_balmer': self.err_H_data_array,
             'K': self.K_data_array,
             'err_K': self.err_K_data_array
            }

        # convert to Pandas DataFrame
        df_collation = pd.DataFrame(data=d)

        # read in phase and star name information
        phase_info = pd.read_csv(self.phase_subdir + \
                                 config["file_names"]["LIST_SPEC_PHASE"], sep=" ")

        # paste phase info into the table of EWs
        phase_array = []
        name_array = []

        # make new column for phase
        df_collation["phase"] = np.nan

        # get the spectrum names from phase_info without the '.dat'
        phase_info_basename = phase_info['Original_spectrum_file_name'].values#.str.split(".",n=1,expand=True)[:][0]

        # loop over each empirical spectrum name and paste the phase and star name into the array
        for q in range(0, len(df_collation['original_spec_file_name'].values)):
            empir_spec_this_one = phase_info['Original_spectrum_file_name'].where(
                phase_info_basename == df_collation['original_spec_file_name'][q]
                ).dropna()
            phase_this_one = phase_info['Phase'].where(
                phase_info_basename == df_collation['original_spec_file_name'][q]
                ).dropna()
            star_name_this_one = phase_info['Star'].where(
                phase_info_basename == df_collation['original_spec_file_name'][q]
                ).dropna()
            #df_collation.iloc[q, "phase"] = phase_this_one.values[0]
            #import ipdb; ipdb.set_trace()
            df_collation.at[q, "phase"] = phase_this_one.values[0]
            df_collation.at[q, "star_name"] = star_name_this_one.values[0]
            # check the name-phase combination is right
            if (df_collation.at[q, "original_spec_file_name"] == empir_spec_this_one.any()):
                continue
            else:
                print("Phases and spectrum names out of order!")
                return

        # drop row of nans (probably redundant)
        df_collation_real = df_collation.dropna().copy(deep=True)

        # write to csv
        df_collation_real.to_csv(self.hk_file_name)
        print('----------------------------------------')
        print('HK data written to ')
        print(str(self.hk_file_name))

        # make plot: each color is a different star, open circles are bad phase region
        data_to_plot = pd.read_csv(self.hk_file_name) # read data back in

        # make list of unique star names
        unique_star_names = df_collation['star_name'].values

        print('unique_star_names')
        print(unique_star_names)

        # plot data points
        cmap = plt.get_cmap(name='jet')
        fig = plt.figure(figsize=(20, 10))
        ax = fig.add_subplot(111)

        # make set of colors/markers I can loop over
        colors = ['red', 'blue', 'orange', 'teal', 'black', 'green', 'purple']*10
        markers = ['o', '^', '>', 's', '<']*10

        '''
        for y in range(0, len(unique_star_names)):
            # loop over every star, overlay the set of points for that star on the plot

            x_data = data_to_plot['balmer'].where(
                data_to_plot['star_name'] == unique_star_names[y])
            y_data = data_to_plot['K'].where(
                data_to_plot['star_name'] == unique_star_names[y])

            err_x_data = data_to_plot['err_balmer'].where(
                data_to_plot['star_name'] == unique_star_names[y])
            err_y_data = data_to_plot['err_K'].where(
                data_to_plot['star_name'] == unique_star_names[y])

            empirical_spectrum_names = data_to_plot['original_spec_file_name'].where(
                data_to_plot['star_name'] == unique_star_names[y])

            #import ipdb; ipdb.set_trace()

            # plot, and keep the same color for each star
            #color_this_star = cmap(float(y)/len(unique_star_names))
            ax.errorbar(x_data,
                        y_data,
                        yerr=err_y_data,
                        xerr=err_x_data,
                        linestyle='',
                        fmt=markers[y],
                        markerfacecolor=colors[y],
                        color=colors[y])


            bad_phase_locs = np.logical_or(data_to_plot['phase'] > self.max_good,
                                           data_to_plot['phase'] < self.min_good)
            x_data_bad_phase = x_data.where(bad_phase_locs)
            y_data_bad_phase = y_data.where(bad_phase_locs)

            # overplot unfilled markers to denote bad phase region
            ax.errorbar(x_data_bad_phase,
                        y_data_bad_phase,
                        linestyle='',
                        fmt=markers[y],
                        markerfacecolor='white',
                        color=colors[y])

            # add star nam
            ax.annotate(unique_star_names[y],
                        xy=(np.array(x_data.dropna())[0],
                            np.array(y_data.dropna())[0]),
                        xytext=(np.array(x_data.dropna())[0],
                                np.array(y_data.dropna())[0]))

            # overplot the name of the empirical spectrum at each data point
            #import ipdb; ipdb.set_trace()
            empirical_spectra_names_this_star = np.array(empirical_spectrum_names.dropna())
            for spec_annotate_num in range(0,len(empirical_spectra_names_this_star)):
                ax.annotate(empirical_spectra_names_this_star[spec_annotate_num],
                        xy=(np.array(x_data.dropna())[spec_annotate_num],
                            np.array(y_data.dropna())[spec_annotate_num]),
                        xytext=(np.array(x_data.dropna())[spec_annotate_num],
                                np.array(y_data.dropna())[spec_annotate_num]),
                        fontsize=6)
        '''
        #import ipdb; ipdb.set_trace()
        # connect lines between each 'star'; that is, with the same gravity and metallicity
        df_20m10 = data_to_plot[data_to_plot['original_spec_file_name'].str.contains('20m10')].sort_values(by=["original_spec_file_name"]).reset_index()
        df_20m15 = data_to_plot[data_to_plot['original_spec_file_name'].str.contains('20m15')].sort_values(by=["original_spec_file_name"]).reset_index()
        df_20m20 = data_to_plot[data_to_plot['original_spec_file_name'].str.contains('20m20')].sort_values(by=["original_spec_file_name"]).reset_index()
        df_20m25 = data_to_plot[data_to_plot['original_spec_file_name'].str.contains('20m25')].sort_values(by=["original_spec_file_name"]).reset_index()
        df_20m30 = data_to_plot[data_to_plot['original_spec_file_name'].str.contains('20m30')].sort_values(by=["original_spec_file_name"]).reset_index()
        df_20p02 = data_to_plot[data_to_plot['original_spec_file_name'].str.contains('20p02')].sort_values(by=["original_spec_file_name"]).reset_index()
        df_25m05 = data_to_plot[data_to_plot['original_spec_file_name'].str.contains('25m05')].sort_values(by=["original_spec_file_name"]).reset_index()
        df_25m10 = data_to_plot[data_to_plot['original_spec_file_name'].str.contains('25m10')].sort_values(by=["original_spec_file_name"]).reset_index()
        df_25m15 = data_to_plot[data_to_plot['original_spec_file_name'].str.contains('25m15')].sort_values(by=["original_spec_file_name"]).reset_index()
        df_25m20 = data_to_plot[data_to_plot['original_spec_file_name'].str.contains('25m20')].sort_values(by=["original_spec_file_name"]).reset_index()
        df_25m25 = data_to_plot[data_to_plot['original_spec_file_name'].str.contains('25m25')].sort_values(by=["original_spec_file_name"]).reset_index()
        df_25m30 = data_to_plot[data_to_plot['original_spec_file_name'].str.contains('25m30')].sort_values(by=["original_spec_file_name"]).reset_index()
        df_30m05 = data_to_plot[data_to_plot['original_spec_file_name'].str.contains('30m05')].sort_values(by=["original_spec_file_name"]).reset_index()
        df_30m10 = data_to_plot[data_to_plot['original_spec_file_name'].str.contains('30m10')].sort_values(by=["original_spec_file_name"]).reset_index()
        df_30m15 = data_to_plot[data_to_plot['original_spec_file_name'].str.contains('30m15')].sort_values(by=["original_spec_file_name"]).reset_index()
        df_30m20 = data_to_plot[data_to_plot['original_spec_file_name'].str.contains('30m20')].sort_values(by=["original_spec_file_name"]).reset_index()
        df_30m25 = data_to_plot[data_to_plot['original_spec_file_name'].str.contains('30m25')].sort_values(by=["original_spec_file_name"]).reset_index()
        df_30m30 = data_to_plot[data_to_plot['original_spec_file_name'].str.contains('30m30')].sort_values(by=["original_spec_file_name"]).reset_index()
        df_30p00 = data_to_plot[data_to_plot['original_spec_file_name'].str.contains('30p00')].sort_values(by=["original_spec_file_name"]).reset_index()
        df_25p00 = data_to_plot[data_to_plot['original_spec_file_name'].str.contains('25p00')].sort_values(by=["original_spec_file_name"]).reset_index()
        df_30p02 = data_to_plot[data_to_plot['original_spec_file_name'].str.contains('30p02')].sort_values(by=["original_spec_file_name"]).reset_index()
        df_25p02 = data_to_plot[data_to_plot['original_spec_file_name'].str.contains('25p02')].sort_values(by=["original_spec_file_name"]).reset_index()
        df_20m05 = data_to_plot[data_to_plot['original_spec_file_name'].str.contains('20m05')].sort_values(by=["original_spec_file_name"]).reset_index()
        df_20p00 = data_to_plot[data_to_plot['original_spec_file_name'].str.contains('20p00')].sort_values(by=["original_spec_file_name"]).reset_index()

        # establish color map
        n = 16
        colors = pl.cm.jet(np.linspace(0,1,n))
        #import ipdb; ipdb.set_trace()
        # definition for making the annotation a bit simpler
        def annotate_fcn(ax_pass, stuff_2_plot):
            for spec_annotate_num in range(0,len(stuff_2_plot)):
                ax_pass.annotate(stuff_2_plot["original_spec_file_name"][spec_annotate_num],
                        xy=(stuff_2_plot["balmer"][spec_annotate_num],stuff_2_plot["K"][spec_annotate_num]),fontsize=6)
        # definition for making the dashed-line plots a bit simpler
        def dashed_line_plot(ax_pass,df_pass):
            ax_pass.errorbar(df_pass["balmer"], df_pass["K"], yerr=df_pass["err_K"], xerr=df_pass["err_balmer"],
                        fmt='o', elinewidth=0.5, ecolor='k', capsize=5, capthick=0.5, linestyle = "--", alpha = 0.2)

        dashed_line_plot(ax,df_20m05)
        annotate_fcn(ax,df_20m05)

        dashed_line_plot(ax,df_20p00)
        annotate_fcn(ax,df_20p00)

        dashed_line_plot(ax,df_20m10)
        annotate_fcn(ax,df_20m10)

        dashed_line_plot(ax,df_20m15)
        annotate_fcn(ax,df_20m15)

        dashed_line_plot(ax,df_20m20)
        annotate_fcn(ax,df_20m20)

        dashed_line_plot(ax,df_20m25)
        annotate_fcn(ax,df_20m25)

        dashed_line_plot(ax,df_20m30)
        annotate_fcn(ax,df_20m30)

        dashed_line_plot(ax,df_20p02)
        annotate_fcn(ax,df_20p02)

        dashed_line_plot(ax,df_25m30)
        annotate_fcn(ax,df_25m30)

        dashed_line_plot(ax,df_25m25)
        annotate_fcn(ax,df_25m25)

        dashed_line_plot(ax,df_25m20)
        annotate_fcn(ax,df_25m20)

        dashed_line_plot(ax,df_25m15)
        annotate_fcn(ax,df_25m15)

        dashed_line_plot(ax,df_25m10)
        annotate_fcn(ax,df_25m10)

        dashed_line_plot(ax,df_25m05)
        annotate_fcn(ax,df_25m05)

        dashed_line_plot(ax,df_25p00)
        annotate_fcn(ax,df_25p00)

        dashed_line_plot(ax,df_25p02)
        annotate_fcn(ax,df_25p02)

        dashed_line_plot(ax,df_30m30)
        annotate_fcn(ax,df_30m30)

        dashed_line_plot(ax,df_30m25)
        annotate_fcn(ax,df_30m25)

        dashed_line_plot(ax,df_30m20)
        annotate_fcn(ax,df_30m20)

        dashed_line_plot(ax,df_30m15)
        annotate_fcn(ax,df_30m15)

        dashed_line_plot(ax,df_30m10)
        annotate_fcn(ax,df_30m10)

        dashed_line_plot(ax,df_30m05)
        annotate_fcn(ax,df_30m05)

        dashed_line_plot(ax,df_30p00)
        annotate_fcn(ax,df_30p00)

        dashed_line_plot(ax,df_30p02)
        annotate_fcn(ax,df_30p02)

        # now remove data for Teff<6000K and Teff>7500K
        df_25m30 = df_25m30.where(np.logical_and(df_25m30["original_spec_file_name"].str[:4].astype(int) >= 6000,
                                      df_25m30["original_spec_file_name"].str[:4].astype(int) <= 7500)).dropna().reset_index()
        df_25m25 = df_25m25.where(np.logical_and(df_25m25["original_spec_file_name"].str[:4].astype(int) >= 6000,
                                      df_25m25["original_spec_file_name"].str[:4].astype(int) <= 7500)).dropna().reset_index()
        df_25m20 = df_25m20.where(np.logical_and(df_25m20["original_spec_file_name"].str[:4].astype(int) >= 6000,
                                      df_25m20["original_spec_file_name"].str[:4].astype(int) <= 7500)).dropna().reset_index()
        df_25m15 = df_25m15.where(np.logical_and(df_25m15["original_spec_file_name"].str[:4].astype(int) >= 6000,
                                      df_25m15["original_spec_file_name"].str[:4].astype(int) <= 7500)).dropna().reset_index()
        df_25m10 = df_25m10.where(np.logical_and(df_25m10["original_spec_file_name"].str[:4].astype(int) >= 6000,
                                      df_25m10["original_spec_file_name"].str[:4].astype(int) <= 7500)).dropna().reset_index()
        df_25m05 = df_25m05.where(np.logical_and(df_25m05["original_spec_file_name"].str[:4].astype(int) >= 6000,
                                      df_25m05["original_spec_file_name"].str[:4].astype(int) <= 7500)).dropna().reset_index()
        df_25p00 = df_25p00.where(np.logical_and(df_25p00["original_spec_file_name"].str[:4].astype(int) >= 6000,
                                      df_25p00["original_spec_file_name"].str[:4].astype(int) <= 7500)).dropna().reset_index()
        df_25p02 = df_25p02.where(np.logical_and(df_25p02["original_spec_file_name"].str[:4].astype(int) >= 6000,
                                      df_25p02["original_spec_file_name"].str[:4].astype(int) <= 7500)).dropna().reset_index()
        df_30m30 = df_30m30.where(np.logical_and(df_30m30["original_spec_file_name"].str[:4].astype(int) >= 6000,
                                      df_30m30["original_spec_file_name"].str[:4].astype(int) <= 7500)).dropna().reset_index()
        df_30m25 = df_30m25.where(np.logical_and(df_30m25["original_spec_file_name"].str[:4].astype(int) >= 6000,
                                      df_30m25["original_spec_file_name"].str[:4].astype(int) <= 7500)).dropna().reset_index()
        df_30m20 = df_30m20.where(np.logical_and(df_30m20["original_spec_file_name"].str[:4].astype(int) >= 6000,
                                      df_30m20["original_spec_file_name"].str[:4].astype(int) <= 7500)).dropna().reset_index()
        df_30m15 = df_30m15.where(np.logical_and(df_30m15["original_spec_file_name"].str[:4].astype(int) >= 6000,
                                      df_30m15["original_spec_file_name"].str[:4].astype(int) <= 7500)).dropna().reset_index()
        df_30m10 = df_30m10.where(np.logical_and(df_30m10["original_spec_file_name"].str[:4].astype(int) >= 6000,
                                      df_30m10["original_spec_file_name"].str[:4].astype(int) <= 7500)).dropna().reset_index()
        df_30m05 = df_30m05.where(np.logical_and(df_30m05["original_spec_file_name"].str[:4].astype(int) >= 6000,
                                      df_30m05["original_spec_file_name"].str[:4].astype(int) <= 7500)).dropna().reset_index()
        df_30p00 = df_30p00.where(np.logical_and(df_30p00["original_spec_file_name"].str[:4].astype(int) >= 6000,
                                      df_30p00["original_spec_file_name"].str[:4].astype(int) <= 7500)).dropna().reset_index()
        df_30p02 = df_30p02.where(np.logical_and(df_30p02["original_spec_file_name"].str[:4].astype(int) >= 6000,
                                      df_30p02["original_spec_file_name"].str[:4].astype(int) <= 7500)).dropna().reset_index()

        # solid line plots, of data which we want to use for the calibration
        ax.errorbar(df_25m30["balmer"], df_25m30["K"], yerr=df_25m30["err_K"], xerr=df_25m30["err_balmer"], linestyle="-", color=colors[0],
                     fmt='o', elinewidth=0.5, ecolor='k', capsize=5, capthick=0.5)
        annotate_fcn(ax,df_25m30)
        ax.errorbar(df_25m25["balmer"], df_25m25["K"], yerr=df_25m25["err_K"], xerr=df_25m25["err_balmer"], linestyle="-", color=colors[1],
                     fmt='o', elinewidth=0.5, ecolor='k', capsize=5, capthick=0.5)
        annotate_fcn(ax,df_25m25)
        ax.errorbar(df_25m20["balmer"], df_25m20["K"], yerr=df_25m20["err_K"], xerr=df_25m20["err_balmer"], linestyle="-", color=colors[2],
                     fmt='o', elinewidth=0.5, ecolor='k', capsize=5, capthick=0.5)
        annotate_fcn(ax,df_25m20)
        ax.errorbar(df_25m15["balmer"], df_25m15["K"], yerr=df_25m15["err_K"], xerr=df_25m15["err_balmer"], linestyle="-", color=colors[3],
                     fmt='o', elinewidth=0.5, ecolor='k', capsize=5, capthick=0.5)
        annotate_fcn(ax,df_25m15)
        ax.errorbar(df_25m10["balmer"], df_25m10["K"], yerr=df_25m10["err_K"], xerr=df_25m10["err_balmer"], linestyle="-", color=colors[4],
                     fmt='o', elinewidth=0.5, ecolor='k', capsize=5, capthick=0.5)
        annotate_fcn(ax,df_25m10)
        ax.errorbar(df_25m05["balmer"], df_25m05["K"], yerr=df_25m05["err_K"], xerr=df_25m05["err_balmer"], linestyle="-", color=colors[5],
                     fmt='o', elinewidth=0.5, ecolor='k', capsize=5, capthick=0.5)
        annotate_fcn(ax,df_25m05)
        ax.errorbar(df_25p00["balmer"], df_25p00["K"], yerr=df_25p00["err_K"], xerr=df_25p00["err_balmer"], linestyle="-", color=colors[6],
                     fmt='o', elinewidth=0.5, ecolor='k', capsize=5, capthick=0.5)
        annotate_fcn(ax,df_25p00)
        ax.errorbar(df_25p02["balmer"], df_25p02["K"], yerr=df_25p02["err_K"], xerr=df_25p02["err_balmer"], linestyle="-", color=colors[7],
                     fmt='o', elinewidth=0.5, ecolor='k', capsize=5, capthick=0.5)
        annotate_fcn(ax,df_25p02)
        ax.errorbar(df_30m30["balmer"], df_30m30["K"], yerr=df_30m30["err_K"], xerr=df_30m30["err_balmer"], linestyle="-", color=colors[8],
                     fmt='o', elinewidth=0.5, ecolor='k', capsize=5, capthick=0.5)
        annotate_fcn(ax,df_30m30)
        ax.errorbar(df_30m25["balmer"], df_30m25["K"], yerr=df_30m25["err_K"], xerr=df_30m25["err_balmer"], linestyle="-", color=colors[9],
                     fmt='o', elinewidth=0.5, ecolor='k', capsize=5, capthick=0.5)
        annotate_fcn(ax,df_30m25)
        ax.errorbar(df_30m20["balmer"], df_30m20["K"], yerr=df_30m20["err_K"], xerr=df_30m20["err_balmer"], linestyle="-", color=colors[10],
                     fmt='o', elinewidth=0.5, ecolor='k', capsize=5, capthick=0.5)
        annotate_fcn(ax,df_30m20)
        ax.errorbar(df_30m15["balmer"], df_30m15["K"], yerr=df_30m15["err_K"], xerr=df_30m15["err_balmer"], linestyle="-", color=colors[11],
                     fmt='o', elinewidth=0.5, ecolor='k', capsize=5, capthick=0.5)
        annotate_fcn(ax,df_30m15)
        ax.errorbar(df_30m10["balmer"], df_30m10["K"], yerr=df_30m10["err_K"], xerr=df_30m10["err_balmer"], linestyle="-", color=colors[12],
                     fmt='o', elinewidth=0.5, ecolor='k', capsize=5, capthick=0.5)
        annotate_fcn(ax,df_30m10)
        ax.errorbar(df_30m05["balmer"], df_30m05["K"], yerr=df_30m05["err_K"], xerr=df_30m05["err_balmer"], linestyle="-", color=colors[13],
                     fmt='o', elinewidth=0.5, ecolor='k', capsize=5, capthick=0.5)
        annotate_fcn(ax,df_30m05)
        ax.errorbar(df_30p00["balmer"], df_30p00["K"], yerr=df_30p00["err_K"], xerr=df_30p00["err_balmer"], linestyle="-", color=colors[14],
                     fmt='o', elinewidth=0.5, ecolor='k', capsize=5, capthick=0.5)
        annotate_fcn(ax,df_30p00)
        ax.errorbar(df_30p02["balmer"], df_30p02["K"], yerr=df_30p02["err_K"], xerr=df_30p02["err_balmer"], linestyle="-", color=colors[15],
                     fmt='o', elinewidth=0.5, ecolor='k', capsize=5, capthick=0.5)
        annotate_fcn(ax,df_30p02)

        #import ipdb; ipdb.set_trace()





        plt.title('KH plot\n(unfilled markers = bad phase region)')
        plt.ylabel('CaIIK EW ($m\AA$?)')
        plt.xlabel('Balmer EW ($m\AA$?)')
        plt.tight_layout()
        plt.savefig(self.plot_write_subdir + config["file_names"]["KH_PLOT_NAME"])

        plt.ylim([0,20])
        plt.savefig(self.plot_write_subdir + "stretched_" + config["file_names"]["KH_PLOT_NAME"])

        plt.close()

        print("HK plots saved as ")
        print(self.plot_write_subdir + config["file_names"]["KH_PLOT_NAME"])
        print(self.plot_write_subdir + "stretched_" + config["file_names"]["KH_PLOT_NAME"])

        # return stuff to enable testing
        return unique_star_names, data_to_plot
