#!/usr/bin/python
'''
This module takes a list of spectra and generates normalized realizations using Gaussian Error.

This module takes a list of spectra. First it generates 100 realizations of the spectra using the error bar to move the points around
simulating Gaussian noise. Then it runs Ken Carrell's bkgrnd routine to determine the normalization and finally creates the normalized spectra.

@package create_spec_realizations
@author deleenm
@version \e \$Revision$
@date \e \$Date$

Usage: create_spec_realizations.py
'''

# -----------------------------
# Standard library dependencies
# -----------------------------
import argparse
import os
from subprocess import Popen,PIPE
import sys
## ## import test.so # experimenting with C extensions
# -------------------
# Third-party imports
# -------------------
from astropy.io import fits
from astropy.table import Table
import numpy as np
from modules2 import *

# --------------------
# Function Definitions
# --------------------
def create_norm_spec(name_list,
                     normdir,
                     finaldir):
    '''
    Create final normalized spectra, using the output from the bkgrnd routine (which puts out wavelength, flux, and continuum flux, but
    not the actual normalized flux)
    
    Arguments:
        name_list: List of Realization file names (no path info)
        normdir: bkgrnd ascii files
        finaldir: The final directory for files.
    Returns:
       A list of final file names
    '''
    
    new_name_list = list()
    
    for spec in name_list: # loop through spectrum realizations
        
        spec_name = os.path.join(normdir,spec) # spectrum realization file name (as output by bkgrnd), with relative path info
        spec_tab = read_bkgrnd_spec(spec_name) # astropy table containing a spectrum's 1.) wavelength, 2.) flux, 3.) background flux
        new_name = os.path.join(finaldir,spec) # output file name of final, normalized spectrum, with relative path info
        new_name_list.append(new_name) # add to list

        try:
            outfile = open(new_name,'w') # open file to write normalized spectrum to
        except IOError:
            print("File {} could not be opened!".format(new_name))
        for j in range(len(spec_tab['wavelength'])):
            outfile.write("{} {:.4f}\n".format(spec_tab['wavelength'][j],spec_tab['flux'][j]/spec_tab['bckgrnd_flux'][j])) # do the division to normalize and write out
        
        outfile.close()
    
    return(new_name_list)  
    
def generate_realizations(spec_name,outdir,num):
    '''
    Calculates a Number of Realizations of a given spectrum using Gaussian Errorbars
    
    Arguments:
        spec_name: The spectrum filename
        outdir: The working directory
        num: Number of realizations to generate
    Returns:
       A list of filenames for the realization spectra.
    '''
    spec_tab = read_spec(spec_name) # astropy table containing an empirical spectrum's 1.) wavelength, 2.) flux, 3.) error
    
    basename = os.path.basename(spec_name) # shave off path stem

    # generate realizations
    new_name_list = list()
    for i in range(num):
        new_name = "{}_{:03d}".format(basename,i) # basename of spectrum realization
        new_name_list.append(new_name) # don't need path info in spec_name list
        new_name = os.path.join(outdir,new_name) # name of spectrum realization, with path
        new_flux = np.random.standard_normal(len(spec_tab))*spec_tab['error'] + spec_tab['flux'] # add Gaussian error to the empirical flux
        try:
            outfile = open(new_name,'w')
        except IOError:
            print("File {} could not be opened!".format(new_name))
        for j in range(len(new_flux)):
            print("Writing out realization file " + os.path.basename(new_name))
            outfile.write("{} {:.2f}\n".format(spec_tab['wavelength'][j],new_flux[j]))
        outfile.close()
    return(new_name_list)
    
def read_bkgrnd_spec(spec_name):
    '''
    Reads in ascii spectra created by bckgrnd and returns numpy arrays of wavelength, flux, bckgrnd_flux
    
    Arguments:
        spec_name: The spectrum filename. If Ascii file should have 3 columns: wavelength, flux, bckgrnd_flux
    Returns:
       A numpy Table with three columns: waveleread_bknght, flus, bckgrnd_flux
       wavelength: Numpy array of wavelengths
       flux: Numpy array of fluxes
       bckgrnd_flux: Numpy array of flux error
    '''
    
    spec_tab = Table.read(spec_name,format='ascii.no_header',names=['wavelength','flux','bckgrnd_flux'])
        
    return(spec_tab)

def read_list(input_list):
    '''
    Reads in list of spectrum names and returns a table of filenamse
    
    Arguments:
        input_list: The spectrum filename
    Returns:
       Numpy array of filenames
    '''   
       
    filenames_arr = np.genfromtxt(input_list,'str',usecols = (0)) # col 0 contains the file names
    return(filenames_arr)

def read_spec(spec_name):
    '''
    Reads in ascii empirical spectra and returns numpy arrays of wavelength, flux, and error.
    
    Arguments:
        spec_name: The spectrum filename. If Ascii file should have 3 columns: wavelength, flux, error no headers
    Returns:
       A numpy Table with three columns: wavelenght, flus, error
       wavelength: Numpy array of wavelengths
       flux: Numpy array of fluxes
       error: Numpy array of flux error
    '''

    spec_tab = Table.read(spec_name,format='ascii.no_header',names=['wavelength','flux','error'])
        
    return(spec_tab)

def write_bckgrnd_input(name_list,indir,normdir):
    '''
    Create input file for the bckgrnd program
    
    Arguments:
        name_list: List of Realization file names (no path info)
        indir: The working directory with files
        normdir: The output directory for normalized files
    Returns:
       A string with the background input filename
    '''
    
    #Check to see if inputfile is already there
    bckgrnd_input = os.path.join(indir,"bckgrnd_input.txt")
    if os.path.isfile(bckgrnd_input) is True:
        os.remove(bckgrnd_input)
    try:
        outfile = open(bckgrnd_input,'w')
    except IOError:
            print("File {} could not be opened!".format(bckgrnd_input))
    
    
    #Write the header (in_dir out_dir)
    outfile.write("{} {}\n".format(indir,normdir))
    for j in range(len(name_list)):
        outfile.write("{}\n".format(name_list[j]))
    outfile.close()
    return(bckgrnd_input)

# -------------
# Main Function
# -------------
def create_spec_realizations_main(num = 100,
                                  input_spec_list_dir = config["data_dirs"]["DIR_SRC"],
                                  unnorm_noise_churned_spectra_dir = config["data_dirs"]["DIR_SYNTH_SPEC"],
                                  bkgrnd_output_dir = config["data_dirs"]["DIR_SYNTH_SPEC_NORM"],
                                  final_dir = config["data_dirs"]["DIR_SYNTH_SPEC_NORM_FINAL"],
                                  verb=False):

    print("--------------------------")
    print("Making "+str(num)+" realizations of each empirical spectrum")
    
    # Read list of empirical spectra
    input_list = input_spec_list_dir + config["file_names"]["LIST_SPEC_PHASE"]
    list_arr = read_list(input_list)
    
    # Check to make sure outdir (to receive realizations of spectra) exists
    outdir = unnorm_noise_churned_spectra_dir
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
           
    # Create realizations for each spectrum
    name_list = list() # initialize
    for i in range(len(list_arr)): # make spectrum realizations and list of their filenames
        name_list.extend(generate_realizations(list_arr[i],outdir,num))
        
    # Create input list of spectrum realization filenames
    bkg_input_file = write_bckgrnd_input(name_list,outdir,bkgrnd_output_dir)
    
    # Normalize each spectrum realization (smoothing parameter is set in __init__)
    bkgrnd = Popen([get_setuptools_script_dir() + "/bkgrnd","--smooth "+str(smooth_val),"--sismoo 1", "--no-plot", "{}".format(bkg_input_file)],stdout=PIPE,stderr=PIPE)
    (out,err) = bkgrnd.communicate() # returns tuple (stdout,stderr)
    
    if verb == True: ## ## decode messages; are they used later? why take this step?
        print(out.decode("utf-8"))
        print(err.decode("utf-8"))
        
    # Normalize spectrum realizations
    final_list = create_norm_spec(name_list, bkgrnd_output_dir, finaldir) # write files of normalized fluxes, and return list of those filenames
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generates normalized spectra realizations using Gaussian Error')
    parser.add_argument('input_list',help='List of spectra to process.')
    parser.add_argument('-o',default='tmpdir', metavar='Output_Dir', help='Output directory (Default tmpdir).')
    parser.add_argument('-n',type=int,default=100,metavar='Num',help='Number of Realizations (Default 100).')
    parser.add_argument('-v',action='store_true',help='Turn on verbosity')   
    #Put this in a dictionary    
    args = vars(parser.parse_args())
    ret = create_spec_realizations_main(args['input_list'],args['o'],args['n'],args['v'])

##
#@mainpage
 #@copydetails  create_spec_realizations
    
