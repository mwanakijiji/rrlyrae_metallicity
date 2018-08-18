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

# -----------------
# Class Definitions
# -----------------

# --------------------
# Function Definitions
# --------------------
def create_norm_spec(name_list,normdir,finaldir):
    '''
    Create final normalized spectra
    
    Arguments:
        name_list: List of Realization file names (no path info)
        normdir: bkgrnd ascii files
        finaldir: The final directory for files.
    Returns:
       A list of final file names
    '''
    
    new_name_list = list()
    
    for spec in name_list:
        spec_name = os.path.join(normdir,spec)
        spec_tab = read_bkgrnd_spec(spec_name)
        new_name = os.path.join(finaldir,spec)
        new_name_list.append(new_name)
        try:
            outfile = open(new_name,'w')
        except IOError:
            print("File {} could not be opened!".format(new_name))
        for j in range(len(spec_tab['wavelength'])):
            outfile.write("{} {:.4f}\n".format(spec_tab['wavelength'][j],spec_tab['flux'][j]/spec_tab['bckgrnd_flux'][j]))
    
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
    spec_tab = read_spec(spec_name)

    basename = os.path.basename(spec_name)

    new_name_list = list()
    #Generate Realizations
    for i in range(num):
        new_name = "{}_{:03d}".format(basename,i)
        new_name_list.append(new_name) #I don't need path info in spec_name list
        new_name = os.path.join(outdir,new_name)
        new_flux = np.random.standard_normal(len(spec_tab))*spec_tab['error'] + spec_tab['flux']
        try:
            outfile = open(new_name,'w')
        except IOError:
            print("File {} could not be opened!".format(new_name))
        for j in range(len(new_flux)):
            outfile.write("{} {:.2f}\n".format(spec_tab['wavelength'][j],new_flux[j]))
        outfile.close()
    return(new_name_list)
    
def read_bkgrnd_spec(spec_name):
    '''
    Reads in ascii spectra created by bckgrnd and returns numpy arrays of wavelength, flux, bckgrnd__flux
    
    Arguments:
        spec_name: The spectrum filename. If Ascii file should have 3 columns: wavelength, flux, bckgrnd_flux
    Returns:
       A numpy Table with three columns: wavelenght, flus, bckgrnd_flux
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
    Reads in ascii spectra and returns numpy arrays of wavelength, flux, and error.
    
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
def create_spec_realizations_main(input_list,outdir,num=100,verb=False):
    #Read list of spectra
    list_arr = read_list(input_list)
    #Check to make sure outdir exists
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
           
    #Create realizations for each spectrum
    name_list = list()
    for i in range(len(list_arr)):
        name_list.extend(generate_realizations(list_arr[i],outdir,num))

    #Check to see if outdir/norm exits
    normdir = os.path.join(outdir,'norm')
    if not os.path.isdir(normdir):
        os.mkdir(normdir)

    #Normalize each realization spectrum
    #Create input file
    bkg_input_file = write_bckgrnd_input(name_list,outdir,normdir)

    #Run bckgrnd
    ## ## NOTE SMOOTHING IS HARDCODED HERE; MAKE THIS VALUE AN INHERITED DEFAULT OF 22
    bkgrnd = Popen(["./bkgrnd","--smooth 22","--sismoo 1", "--no-plot", "{}".format(bkg_input_file)],stdout=PIPE,stderr=PIPE)

    ## ## start attempt to use *.so file as executable
    #gcc -o yourexecutable objects/yourobject1.o objects/yourobject2.o -lLibrary
    #bkgrnd = Popen(["./bin/RRab-metallicity/bkgrnd13.cpython-35m-darwin.so","--smooth 22","--sismoo 1", "--no-plot", "{}".format(bkg_input_file)],stdout=PIPE,stderr=PIPE)
    ## ## end attempt
    
    (out,err) = bkgrnd.communicate()
    
    if verb == True:
        print(out.decode("utf-8"))
        print(err.decode("utf-8"))
    
    #Normalize spectrum
    #Check to see if outdir/final exists
    finaldir = os.path.join(outdir,'final')
    if not os.path.isdir(finaldir):
        os.mkdir(finaldir)
        
    final_list = create_norm_spec(name_list, normdir, finaldir)

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
    
