#!/usr/bin/python
'''
This module takes a list of spectra and normalizes them, without churning the spectrum noise.

@package create_spec_realizations
@author deleenm (of parent create_spec_realizations)
@version \e \$Revision$
@date \e \$Date$

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
from rrlyrae_metallicity.modules2 import *

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
       
    filenames_arr = np.genfromtxt(input_list, 'str', skip_header=1, usecols = (0)) # col 0 contains the file names
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
def normalize_simple(input_spec_list_dir = config_apply["data_dirs"]["DIR_SRC"],
                     unnorm_science_spectra_dir = config_apply["data_dirs"]["DIR_SCI_SPECTRA"],
                     bkgrnd_output_dir = config_apply["data_dirs"]["DIR_SCI_SPEC_NORM"],
                     final_dir = config_apply["data_dirs"]["DIR_SCI_SPEC_NORM_FINAL"],
                     verb = False):
    '''
    INPUTS:
    input_spec_list_dir: directory containing a file containing a list of names of science
        spectra which we want to normalize and find the Fe/H from
    unnorm_science_spectra_dir: directory which actually contains the science spectra
    bkgrnd_output_dir: directory to contain output of bkgrnd (spectra and fit continuua)
    final_dir: directory to contain normalized spectrum realizations

    OUTPUTS:
    (text files written)
    '''

    print("--------------------------")
    print("Normalizing spectra")
    
    # Read list of science spectra to normalize and apply calibration to
    input_list = input_spec_list_dir + config_apply["file_names"]["LIST_SPEC_APPLY"]
    list_arr = read_list(input_list)

    # make list of stemless filenames for the normalized spectra
    name_list = [s + "_norm" for s in list_arr]

    # create list of full filenames of normalized spectra which will be written
    bkg_input_file = write_bckgrnd_input(name_list = name_list,
                                         indir = unnorm_science_spectra_dir,
                                         normdir = bkgrnd_output_dir)
    
    # Normalize each spectrum (smoothing parameter is set in __init__)
    bkgrnd = Popen([str(config_apply["data_dirs"]["DIR_BIN"]) + "bkgrnd", "--smooth "+str(["reduc_params"]["SMOOTH"]),
                    "--sismoo 1", "--no-plot", "{}".format(bkg_input_file)], stdout=PIPE, stderr=PIPE)
    (out,err) = bkgrnd.communicate() # returns tuple (stdout,stderr)
    
    if verb == True:
        print(out.decode("utf-8"))
        print(err.decode("utf-8"))
        
    # write files of normalized fluxes, and return list of those filenames
    final_list = create_norm_spec(name_list, bkgrnd_output_dir, finaldir)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Generates normalized spectra realizations using Gaussian Error')
    parser.add_argument('input_list', help='List of spectra to process.')
    parser.add_argument('-o', default='tmpdir', metavar='Output_Dir', help='Output directory (Default tmpdir).')
    parser.add_argument('-n', type=int, default=100, metavar='Num', help='Number of Realizations (Default 100).')
    parser.add_argument('-v', action='store_true', help='Turn on verbosity')   
    #Put this in a dictionary    
    args = vars(parser.parse_args())
    ret = normalize_simple(args['input_list'],args['o'],args['n'],args['v'])
