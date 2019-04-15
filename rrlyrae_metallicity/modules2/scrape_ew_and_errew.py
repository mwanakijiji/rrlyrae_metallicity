import glob
import pandas as pd
import os
import numpy as np
import sys
import matplotlib.pyplot as plt
from rrlyrae_metallicity.modules2 import *

class Scraper():
    '''
    Scrape all the equivalent width info from the Robospect *.fits.robolines files

    INPUTS:
    subdir: subdirectory containing the Robospect output
    '''
        
    def __init__(self, subdir = config["data_dirs"]["DIR_ROBO_OUTPUT"], verbose=False):
        
        # directory containing the directory containing the *.fits.robolines files containing the EW info
        self.stem = '.' ## ##
        # subdirectory containing the *.c.dat files
        self.subdir = subdir ## ##
        
        # get list of filenames without the path
        fileListLong = glob.glob(self.subdir+'/'+'*.fits.robolines')
        fileListUnsorted = [os.path.basename(x) for x in fileListLong]
        self.fileList = sorted(fileListUnsorted)

        # EW info will get scraped into this
        self.writeOutFilename = self.subdir+config["file_names"]["MCD_LARGE_BAD_REMOVED"]

        # return tables of EW data?
        self.verbose = verbose
        
    def __call__(self):

        def line_order_check(line_centers):
            '''
            Sanity check: are the lines listed in order?
            N.b. This checks the wavelengths using the given line list values (and not the fitted centers)
            '''
            
            if ((line_centers[0] < 3933.660-10) or (line_centers[0] > 3933.660+10)): # CaIIK
                print('Lines not matching!')
                sys.exit()  # ... and abort
            elif ((line_centers[1] < 3970.075-10) or (line_centers[1] > 3970.075+10)): # H-epsilon (close to CaIIH)
                print('Lines not matching!')
                sys.exit()
            elif ((line_centers[2] < 4101.7100-10) or (line_centers[2] > 4101.7100+10)): # H-delta
                print('Lines not matching!')
                sys.exit()
            elif ((line_centers[3] < 4340.472-10) or (line_centers[3] > 4340.472+10)): # H-gamma
                print('Lines not matching!')
                sys.exit()
            elif ((line_centers[4] < 4861.290-10) or (line_centers[4] > 4861.290+10)): # H-beta
                print('Lines not matching!')
                sys.exit()
            return

        dfMaster = pd.DataFrame() # initialize

        # loop over all filenames of realizations of empirical spectra, extract line data
        for t in range(0,len(self.fileList)):

            # read in Robospect output
            df = pd.read_csv(self.subdir+'/'+self.fileList[t], header=13, delim_whitespace=True, index_col=False, usecols=np.arange(17))
    
            # check lines are in the right order
            line_order_check(df['#x0'])
    
            # add two cols on the left: the filename, and the name of the line
            sLength = len(df['mean']) # number of lines (should be 5)

            # file names
            df['file_name'] = pd.Series(self.fileList[t], index=df.index)

            # names of empirical spectra realizations (multiple ones correspond to one empirical spectrum)
            df['synth_spec_name'] = pd.Series(self.fileList[t].split(".")[0], index=df.index)

            # names of empirical spectra
            df['empir_spec_name'] = pd.Series(self.fileList[t].split(".")[0][0:-4], index=df.index)
            
            #df['star_name'] = pd.Series(self.fileList[t].split("__")[0], index=df.index)

            # names of the absorption lines
            df['line_name'] = ['CaIIK', 'Heps', 'Hdel', 'Hgam', 'Hbet']
    
            # print progress
            print('Out of '+str(len(self.fileList))+' files, '+str(t)+' scraped...')
    
            # if this is the first list, start a master copy from it to concatenate stuff to it
            if (t==0):
                dfMaster = df.copy()
            else:
                dfMaster = pd.concat([dfMaster,df])
                del df # clear variable

        # write to csv, while resetting the indices
        # note THIS TABLE INCLUDES ALL DATA, GOOD AND BAD
        dfMaster_reset = dfMaster.reset_index(drop=True).copy() # this gets shown further down in this notebook
        dfMaster.reset_index(drop=True).to_csv(self.subdir+config["file_names"]["MCD_ALL_DATA"]) # this is effectively the same, but gets written out

        ## IF WE ARE INTERESTED IN SPECTRA THAT HAVE ALL WELL-FIT LINES
        # remove all rows with a flag ending with something other than zero (i.e., the fit is bad)
        # make an array consisting of the last character in each spectrum's flag
        redFlagArray = ([u[-1] for u in dfMaster_reset.flags])
        # consider bad flags to be of any flag with a nonzero last character
        whereRedFlag = np.where(np.array(redFlagArray) != '0')
        
        # identify the synthetic spectrum names which have at least one line with a bad fit
        badSynthSpectra = dfMaster_reset['synth_spec_name'][np.squeeze(whereRedFlag)]
        # remove duplicate names
        badSynthSpectra_uniq = badSynthSpectra.drop_duplicates()
        # keep only the spectra that have all lines well-fit
        dfMaster_reset_dropBadSpectra = dfMaster_reset.where(~dfMaster_reset['synth_spec_name'].isin(badSynthSpectra_uniq))
        
        # write to csv
        # note THIS TABLE HAS SPECTRA WITH ANY BAD ROWS REMOVED
        dfMaster_reset_dropBadSpectra.to_csv(self.writeOutFilename) # this is effectively the same, but gets written out

        print("--------------------------")
        print('Scraped Robospect output written to')
        print(self.writeOutFilename)

        if self.verbose:
            return dfMaster_reset, dfMaster_reset_dropBadSpectra


class findHK():
    '''
    Read in Robospect EWs of synthetic spectra, rescale them, average them, and plot KH space
    '''
    
    def __init__(self,
                 source_subdir = config["data_dirs"]["DIR_ROBO_OUTPUT"],
                 phase_subdir = config["data_dirs"]["TEST_DIR_SRC"],
                 plot_write_subdir = config["data_dirs"]["DIR_FYI_INFO"],
                 verbose=False):

        self.source_subdir = source_subdir
        self.phase_subdir = phase_subdir
        self.scrapedEWfilename = self.source_subdir + config["file_names"]["MCD_LARGE_BAD_REMOVED"]
        self.hkFileName = config["data_dirs"]["DIR_SRC"] + config["file_names"]["MORE_REALISTIC"]
        
        # read in line data
        print(self.scrapedEWfilename)
        self.line_data = pd.read_csv(self.scrapedEWfilename, delim_whitespace=False)
        
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
        self.min_good, self.max_good = phase_regions()

        # indicate subdirectory where FYI plot will be written
        self.plot_write_subdir = plot_write_subdir
        
    def __call__(self):
        
        # make a list of all UNIQUE, EMPIRICAL spectrum names
        uniqueSpecNames_preIndexReset = self.line_data.drop_duplicates(subset='empir_spec_name')['empir_spec_name']

        # drop row of NaNs and smooth the indexing
        uniqueSpecNames = uniqueSpecNames_preIndexReset.dropna().reset_index(drop=True) 
        
        # fit a straight line to Hgam vs Hdel
        x_data = self.line_data['EQW'].where(self.line_data['line_name'] == 'Hdel').dropna() # Hdel
        y_data = self.line_data['EQW'].where(self.line_data['line_name'] == 'Hgam').dropna() # Hgam
        Hgam = np.copy(y_data)
        m,b = np.polyfit(x_data, y_data, 1) # might want errors later, too 
        
        # generate a rescaled Hgam, call it rHgam
        rHgam_all = np.divide(np.subtract(Hgam,b),m)
        
        # prepare data for a plot
        # loop over every EMPIRICAL spectrum and assemble SYNTHETIC data into arrays
        for p in range(0,len(np.array(uniqueSpecNames.values))):
            print(len(np.array(uniqueSpecNames.values)))

            print("Synthetic data being assembled corresponding to spectrum "+np.array(uniqueSpecNames)[p])
            # extract all synthetic data corresponding to this empirical spectrum
            data_for_this_empir_spectrum = self.line_data.where(self.line_data['empir_spec_name'][0:-4] == np.array(uniqueSpecNames)[p])

            print('1')
            # scrape data
            raw_Hbet_data = data_for_this_empir_spectrum['EQW'].where(self.line_data['line_name'] == 'Hbet')
            raw_Hgam_data = data_for_this_empir_spectrum['EQW'].where(self.line_data['line_name'] == 'Hgam')
            raw_Hdel_data = data_for_this_empir_spectrum['EQW'].where(self.line_data['line_name'] == 'Hdel')
            raw_Heps_data = data_for_this_empir_spectrum['EQW'].where(self.line_data['line_name'] == 'Heps')
            raw_K_data = data_for_this_empir_spectrum['EQW'].where(self.line_data['line_name'] == 'CaIIK')

            print('2')
            # rescale and remove nans
            Hbet_data_wnans = np.array(np.copy(raw_Hbet_data))
            Hgam_data_wnans = np.array(np.copy(raw_Hgam_data))
            Hdel_data_wnans = np.array(np.copy(raw_Hdel_data))
            Heps_data_wnans = np.array(np.copy(raw_Heps_data))    
            K_data_wnans = np.array(np.copy(raw_K_data))
            rHgam_data_wnans = np.array(np.divide(np.subtract(raw_Hgam_data,b),m)) # rescale Hgam EWs

            print('3')
            Hbet_data = Hbet_data_wnans[np.isfinite(Hbet_data_wnans)] # remove nans
            Hgam_data = Hgam_data_wnans[np.isfinite(Hgam_data_wnans)]
            Hdel_data = Hdel_data_wnans[np.isfinite(Hdel_data_wnans)]
            Heps_data = Heps_data_wnans[np.isfinite(Heps_data_wnans)]
            rHgam_data = rHgam_data_wnans[np.isfinite(rHgam_data_wnans)]
            K_data = K_data_wnans[np.isfinite(K_data_wnans)]

            print('4')
            # get the H-K synthetic data together
            balmer_data_allsynthetic_spec = np.nanmean([Hdel_data,rHgam_data], axis=0) # Balmer EW = 0.5*(Hdel + rHgam)
            K_data_allsynthetic_spec = np.copy(K_data)

            print('5')
            # the actual points to plot (or record in a table)
            Hbet_data_pt = np.nanmedian(Hbet_data)
            Hgam_data_pt = np.nanmedian(Hgam_data)
            rHgam_data_pt = np.nanmedian(rHgam_data)
            Hdel_data_pt = np.nanmedian(Hdel_data)
            Heps_data_pt = np.nanmedian(Heps_data)
            balmer_data_pt = np.nanmedian(balmer_data_allsynthetic_spec)
            K_data_pt = np.nanmedian(K_data_allsynthetic_spec)

            print('6')
            # the error bars
            err_Hbet_data = np.nanstd(Hbet_data)
            err_Hgam_data = np.nanstd(Hgam_data)
            err_rHgam_data = np.nanstd(rHgam_data)
            err_Hdel_data = np.nanstd(Hdel_data)
            err_Heps_data = np.nanstd(Heps_data)
            err_balmer_data = np.nanstd(balmer_data_allsynthetic_spec)
            err_K_data = np.nanstd(K_data_allsynthetic_spec)
    
            #plt.plot(balmer_data_pt,K_data_pt)
            #plt.errorbar(balmer_data_pt, K_data_pt, yerr=err_K_data, xerr=err_balmer_data)

            print('7')
            # append data to arrays: essential info
            self.empir_spec_name_array = np.append(self.empir_spec_name_array,np.array(uniqueSpecNames)[p])
            self.star_name_array = np.append(self.star_name_array,str(np.array(uniqueSpecNames)[p])[0:-3])
            self.H_data_array = np.append(self.H_data_array,balmer_data_pt)
            self.err_H_data_array = np.append(self.err_H_data_array,err_balmer_data)
            self.K_data_array = np.append(self.K_data_array,K_data_pt)
            self.err_K_data_array = np.append(self.err_K_data_array,err_K_data)
    
            # append data to arrays: other info
            self.Hbet_data_array = np.append(self.Hbet_data_array,Hbet_data_pt)
            self.err_Hbet_data_array = np.append(self.err_Hbet_data_array,err_Hbet_data)
            self.Hgam_data_array = np.append(self.Hgam_data_array,Hgam_data_pt)
            self.err_Hgam_data_array = np.append(self.err_Hgam_data_array,err_Hgam_data)
            self.rHgam_data_array = np.append(self.rHgam_data_array,rHgam_data_pt) # rescaled Hgamma
            self.err_rHgam_data_array = np.append(self.err_rHgam_data_array,err_rHgam_data)
            self.Hdel_data_array = np.append(self.Hdel_data_array,Hdel_data_pt)
            self.err_Hdel_data_array = np.append(self.err_Hdel_data_array,err_Hdel_data)
            self.Heps_data_array = np.append(self.Heps_data_array,Heps_data_pt)
            self.err_Heps_data_array = np.append(self.err_Heps_data_array,err_Heps_data)

            print('7')
            # clear some variables
            balmer_data_allsynthetic_spec=None 
            K_data_allsynthetic_spec=None
            balmer_data_allsynthetic_spec=None 
            K_data_allsynthetic_spec=None
            
        # put everything into a dataframe

        print('8')
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
        df_collation = pd.DataFrame(data=d)

        print('9')
        # read in a text file containing phase information
        ## ## (NO- THIS SHOULD BE READ IN AT THE BEGINNING OF THE PIPELINE)
        phase_info = pd.read_csv(self.phase_subdir + config["file_names"]["DETACHED_PHASES"])
        
        # paste phase info into the table of EWs
        phase_array = []
        feh_array = []
        err_feh_array = []
        name_array = []

        print('10')
        # loop over each empirical spectrum name and populate the arrays
        for q in range(0,len(df_collation['empir_spec_name'].values)):    
            name_this_one = phase_info['Spectrum'].where(phase_info['Spectrum'] == df_collation['empir_spec_name'][q]).dropna()
            phase_this_one = phase_info['phase'].where(phase_info['Spectrum'] == df_collation['empir_spec_name'][q]).dropna()
            feh_this_one = phase_info['FeH'].where(phase_info['Spectrum'] == df_collation['empir_spec_name'][q]).dropna()
            err_feh_this_one = phase_info['eFeH'].where(phase_info['Spectrum'] == df_collation['empir_spec_name'][q]).dropna()
            name_array = np.append(name_array,name_this_one)
            phase_array = np.append(phase_array,phase_this_one)
            feh_array = np.append(feh_array,feh_this_one)
            err_feh_array = np.append(err_feh_array,err_feh_this_one)
            #print('-----------')
            #print(name_this_one)
            #print(phase_this_one)
            #print(feh_this_one)
            #print(err_feh_this_one)

        print('11')
        print(df_collation)
        df_collation_real = df_collation.dropna().copy(deep=True) # drop row of nans (probably redundant)
        print(df_collation_real)
        ## ## THIS IS AN ERSATZ FOR NOW; REAL PHASE VALUES NEED TO BE READ IN AT BEGINNING
        #phase_array = 
        df_collation_real['phase'] = phase_array
        df_collation_real['FeH'] = feh_array
        df_collation_real['eFeH'] = err_feh_array

        print('11a')
        # write to csv
        df_collation_real.to_csv(self.hkFileName)
        print('----------------------------------------')
        print('HK data written to ' + self.hkFileName)

        print('11b')
        # make plot: each color is a different star, open circles are bad phase region
        data_to_plot = pd.read_csv(self.hkFileName) # read data back in

        print('11c')
        # make list of unique star names 
        unique_star_names = data_to_plot.drop_duplicates(subset=['star_name'])['star_name'].values

        print('11d')
        # plot data points
        cmap = plt.get_cmap(name='jet')
        fig = plt.figure(figsize=(20,10))
        ax = fig.add_subplot(111)

        print('12')
        # loop over every star, overlay the set of points for that star on the plot
        colors = ['red', 'blue', 'orange', 'teal', 'black', 'green', 'purple']*10 # make set of colors/markers I can loop over
        markers = ['o', '^', '>', 's', '<']*10
        for y in range(0,len(unique_star_names)):
    
            x_data = data_to_plot['balmer'].where(data_to_plot['star_name'] == unique_star_names[y])
            y_data = data_to_plot['K'].where(data_to_plot['star_name'] == unique_star_names[y])
    
            err_x_data = data_to_plot['err_balmer'].where(data_to_plot['star_name'] == unique_star_names[y])
            err_y_data = data_to_plot['err_K'].where(data_to_plot['star_name'] == unique_star_names[y])
    
            # plot, and keep the same color for each star
            color_this_star = cmap(float(y)/len(unique_star_names))
            ax.errorbar(x_data,y_data,yerr=err_y_data,xerr=err_x_data,linestyle='',fmt=markers[y],markerfacecolor=colors[y],color=colors[y])

            bad_phase_locs = np.logical_or(data_to_plot['phase'] > self.max_good, data_to_plot['phase'] < self.min_good)
            x_data_badPhase = x_data.where(bad_phase_locs)
            y_data_badPhase = y_data.where(bad_phase_locs)
    
            # overplot unfilled markers to denote bad phase region
            ax.errorbar(x_data_badPhase,y_data_badPhase,linestyle='',fmt=markers[y],markerfacecolor='white',color=colors[y])
    
            # add star name
            ax.annotate(unique_star_names[y], xy=(np.array(x_data.dropna())[0], 
                                          np.array(y_data.dropna())[0]), 
                xytext=(np.array(x_data.dropna())[0], np.array(y_data.dropna())[0]))

            
            
        plt.title('KH plot\n(unfilled markers = bad phase region)')
        plt.ylabel('CaIIK EW (milliangstrom)')
        plt.xlabel('Balmer EW (milliangstrom)')
        plt.tight_layout()
        plt.savefig(self.plot_write_subdir + config["file_names"]["KH_PLOT_NAME"])
        plt.close()

        # return stuff to enable testing
        return unique_star_names, data_to_plot
