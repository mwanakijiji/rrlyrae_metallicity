'''
Scrape Robospect output and do some processing of the results
'''

import os
import sys
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from rrlyrae_metallicity.modules2 import *

class Scraper():
    '''
    Scrape all the equivalent width info from the Robospect *.fits.robolines files
    '''

    def __init__(self,
                 subdir=config["data_dirs"]["DIR_ROBO_OUTPUT"],
                 file_scraped_info=config["file_names"]["MCD_LARGE_BAD_REMOVED"],
                 verbose=False):

        # directory containing the directory containing the *.fits.robolines
        # files containing the EW info
        self.stem = '.' ## ##
        # subdirectory containing the *.c.dat files
        self.subdir = subdir ## ##

        # get list of filenames without the path
        file_list_long = glob.glob(self.subdir+'/'+'*.fits.robolines')
        file_list_unsorted = [os.path.basename(x) for x in file_list_long]
        self.file_list = sorted(file_list_unsorted)

        # EW info will get scraped into this
        self.write_out_filename = self.subdir + self.file_scraped_info

        # return tables of EW data?
        self.verbose = verbose

    def __call__(self):

        def line_order_check(line_centers):
            '''
            Sanity check: are the lines listed in order?
            N.b. This checks the wavelengths using the given line list
            values (and not the fitted centers)
            '''

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

        # loop over all filenames of realizations of empirical spectra, extract line data
        for t in range(0, len(self.file_list)):

            # read in Robospect output
            df = pd.read_csv(self.subdir+'/'+self.file_list[t],
                             header=13,
                             delim_whitespace=True,
                             index_col=False,
                             usecols=np.arange(17))

            # check lines are in the right order
            line_order_check(df['#x0'])

            # add two cols on the left: the filename, and the name of the line
            #s_length = len(df['mean']) # number of lines (should be 5)

            # file names
            df['file_name'] = pd.Series(self.file_list[t],
                                        index=df.index)

            # names of empirical spectra realizations (multiple ones
            # correspond to one empirical spectrum)
            df['synth_spec_name'] = pd.Series(self.file_list[t].split(".")[0],
                                              index=df.index)

            # names of empirical spectra
            df['empir_spec_name'] = pd.Series(self.file_list[t].split(".")[0][0:-4],
                                              index=df.index)

            #df['star_name'] = pd.Series(self.file_list[t].split("__")[0], index=df.index)

            # names of the abif (sorption lines
            df['line_name'] = ['CaIIK', 'Heps', 'Hdel', 'Hgam', 'Hbet']

            # print progress
            print('Out of '+str(len(self.file_list))+' files, '+str(t)+' scraped...')

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
        bad_synth_spectra = df_master_reset['synth_spec_name'][np.squeeze(where_red_flag)]
        # remove duplicate names
        bad_synth_spectra_uniq = bad_synth_spectra.drop_duplicates()
        # keep only the spectra that have all lines well-fit
        df_master_reset_drop_bad_spectra = df_master_reset.where(
            ~df_master_reset['synth_spec_name'].isin(bad_synth_spectra_uniq)
            )

        # write to csv
        # note THIS TABLE HAS SPECTRA WITH ANY BAD ROWS REMOVED
        df_master_reset_drop_bad_spectra.to_csv(self.write_out_filename)

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
        self.empir_spec_name_array = []
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
            subset='empir_spec_name')['empir_spec_name']

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
        m, b = np.polyfit(x_data, y_data, 1) # might want errors later, too

        # generate a rescaled Hgam, call it rHgam
        rHgam_all = np.divide(np.subtract(Hgam, b), m)

        # prepare data for a plot
        # loop over every EMPIRICAL spectrum and assemble SYNTHETIC data into arrays
        for p in range(0, len(np.array(unique_spec_names.values))):

            print("Synthetic data being assembled corresponding to spectrum "+\
                  np.array(unique_spec_names)[p])
            # extract all synthetic data corresponding to this empirical spectrum
            data_for_this_empir_spectrum = self.line_data.where(
                self.line_data['empir_spec_name'][0:-4] == np.array(unique_spec_names)[p]
                )

            # scrape data
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

            # rescale
            Hbet_data_wnans = np.array(np.copy(raw_Hbet_data))
            Hgam_data_wnans = np.array(np.copy(raw_Hgam_data))
            Hdel_data_wnans = np.array(np.copy(raw_Hdel_data))
            Heps_data_wnans = np.array(np.copy(raw_Heps_data))
            K_data_wnans = np.array(np.copy(raw_K_data))
            # rescale Hgam EWs
            rHgam_data_wnans = np.array(np.divide(
                np.subtract(raw_Hgam_data, b), m))

            # remove nans
            Hbet_data = Hbet_data_wnans[np.isfinite(Hbet_data_wnans)]
            Hgam_data = Hgam_data_wnans[np.isfinite(Hgam_data_wnans)]
            Hdel_data = Hdel_data_wnans[np.isfinite(Hdel_data_wnans)]
            Heps_data = Heps_data_wnans[np.isfinite(Heps_data_wnans)]
            rHgam_data = rHgam_data_wnans[np.isfinite(rHgam_data_wnans)]
            K_data = K_data_wnans[np.isfinite(K_data_wnans)]

            # get the H-K synthetic data together
            # (note balmer EW = 0.5*(Hdel + rHgam) )
            balmer_data_allsynthetic_spec = np.nanmean([Hdel_data, rHgam_data], axis=0)
            K_data_allsynthetic_spec = np.copy(K_data)

            # the actual points to plot (or record in a table)
            Hbet_data_pt = np.nanmedian(Hbet_data)
            Hgam_data_pt = np.nanmedian(Hgam_data)
            rHgam_data_pt = np.nanmedian(rHgam_data)
            Hdel_data_pt = np.nanmedian(Hdel_data)
            Heps_data_pt = np.nanmedian(Heps_data)
            balmer_data_pt = np.nanmedian(balmer_data_allsynthetic_spec)
            K_data_pt = np.nanmedian(K_data_allsynthetic_spec)

            # the error bars
            err_Hbet_data = np.nanstd(Hbet_data)
            err_Hgam_data = np.nanstd(Hgam_data)
            err_rHgam_data = np.nanstd(rHgam_data)
            err_Hdel_data = np.nanstd(Hdel_data)
            err_Heps_data = np.nanstd(Heps_data)
            err_balmer_data = np.nanstd(balmer_data_allsynthetic_spec)
            err_K_data = np.nanstd(K_data_allsynthetic_spec)

            #plt.plot(balmer_data_pt, K_data_pt)
            #plt.errorbar(balmer_data_pt, K_data_pt, yerr=err_K_data, xerr=err_balmer_data)

            # append data to arrays: essential info
            self.empir_spec_name_array = np.append(self.empir_spec_name_array,
                                                   np.array(unique_spec_names)[p])
            self.star_name_array = np.append(self.star_name_array,
                                             str(np.array(unique_spec_names)[p])[0:-3])
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
        d = {'empir_spec_name': self.empir_spec_name_array,
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

        # read in phase information
        phase_info = pd.read_csv(self.phase_subdir + \
                                 config["file_names"]["LIST_SPEC_PHASE"], sep=" ")

        # paste phase info into the table of EWs
        phase_array = []
        name_array = []

        # make new column for phase
        df_collation["phase"] = np.nan

        # get the spectrum names from phase_info without the '.dat'
        phase_info_basename = phase_info['Spectrum'].str.split(".",
                                                               n=1,
                                                               expand=True)[:][0]

        # loop over each empirical spectrum name and paste the phase into the array
        for q in range(0, len(df_collation['empir_spec_name'].values)):
            empir_spec_this_one = phase_info['Spectrum'].where(
                phase_info_basename == df_collation['empir_spec_name'][q]
                ).dropna()
            phase_this_one = phase_info['Phase'].where(
                phase_info_basename == df_collation['empir_spec_name'][q]
                ).dropna()
            #df_collation.iloc[q, "phase"] = phase_this_one.values[0]
            df_collation.at[q, "phase"] = phase_this_one.values[0]
            # check the name-phase combination is right
            if (df_collation.at[q, "empir_spec_name"] == empir_spec_this_one.any().split(".")[0]):
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
        unique_star_names = data_to_plot.drop_duplicates(
            subset=['star_name'])['star_name'].values

        # plot data points
        cmap = plt.get_cmap(name='jet')
        fig = plt.figure(figsize=(20, 10))
        ax = fig.add_subplot(111)

        # make set of colors/markers I can loop over
        colors = ['red', 'blue', 'orange', 'teal', 'black', 'green', 'purple']*10
        markers = ['o', '^', '>', 's', '<']*10
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

            # add star name
            ax.annotate(unique_star_names[y],
                        xy=(np.array(x_data.dropna())[0],
                            np.array(y_data.dropna())[0]),
                        xytext=(np.array(x_data.dropna())[0],
                                np.array(y_data.dropna())[0]))



        plt.title('KH plot\n(unfilled markers = bad phase region)')
        plt.ylabel('CaIIK EW (milliangstrom)')
        plt.xlabel('Balmer EW (milliangstrom)')
        plt.tight_layout()
        plt.savefig(self.plot_write_subdir + config["file_names"]["KH_PLOT_NAME"])
        plt.close()

        # return stuff to enable testing
        return unique_star_names, data_to_plot
