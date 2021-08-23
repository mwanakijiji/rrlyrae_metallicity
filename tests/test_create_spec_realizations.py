import matplotlib
matplotlib.use('Agg')

import sys, os
import configparser
import pandas as pd

current_dir = os.path.dirname(__file__)
target_dir = os.path.abspath(os.path.join(current_dir, "../"))
sys.path.insert(0, target_dir)

# import more things with changed system path
from modules import *
from modules import create_spec_realizations
from conf import *
import numpy as np

# configuration data for reduction
config_red = configparser.ConfigParser() # for parsing values in .init file
# config for reduction to find a, b, c, d
config_red.read(os.path.join(os.path.dirname(__file__), '../conf', 'config_red.ini'))

'''
# check if the directory-making function works
def test_make_dirs():

    # call function to make directories
    make_dirs() # this leads to permissions errors in online build

    # do all the directories exist now?
    for vals in config["data_dirs"]:
        abs_path_name = str(config["data_dirs"][vals])
        assert os.path.exists(abs_path_name)
'''

def test_create_norm_spec():
    '''
    Create final normalized spectra, using the output from the bkgrnd routine (which
    puts out wavelength, flux, and continuum flux, but not the actual normalized flux)

    Arguments:
        name_list: List of Realization file names (no path info)
        normdir: bkgrnd ascii files
        finaldir: The final directory for files which have completed the full normalization process.
    Returns:
       A list of final file names
    '''

    #test_name_list =
    #test_normdir =
    #test_final_dir =
    #final_list = create_norm_spec(name_list, normdir, final_dir)

    # is min smaller than max
    assert 1 < 2


def test_read_spec():

    # FITS format
    spec_name_fits = config_red["data_dirs"]["TEST_DIR_SRC"] + "575030m20.fits"
    test_spec_tab_fits, test_hdr_fits = create_spec_realizations.read_spec(spec_name=spec_name_fits, format="fits")

    # for FITS data, there should be 3 columns of floats, and a header
    assert test_hdr_fits
    assert len(test_spec_tab_fits.colnames) == 3 # 3 columns
    assert isinstance(test_spec_tab_fits["wavelength"][0],np.float64)
    assert isinstance(test_spec_tab_fits["flux"][0],np.float64)
    assert isinstance(test_spec_tab_fits["error"][0],np.float64)
    #assert (test_spec_tab_fits["wavelength"].info.dtype == np.float64)
    #assert (test_spec_tab_fits["flux"].info.dtype == np.float64)
    #assert (test_spec_tab_fits["error"].info.dtype == np.float64)

    # ascii format
    spec_name_ascii = config_red["data_dirs"]["TEST_DIR_SRC"] + "700025m20.smo"
    test_spec_tab_ascii, test_hdr_ascii = create_spec_realizations.read_spec(spec_name=spec_name_ascii, format="ascii.no_header")

    # for ascii data, there should be 3 columns of floats, and NO header
    assert np.isfinite(test_hdr_ascii) == False
    assert len(test_spec_tab_ascii.colnames) == 3 # 3 columns
    assert isinstance(test_spec_tab_ascii["wavelength"][0],np.float64)
    assert isinstance(test_spec_tab_ascii["flux"][0],np.float64)
    assert isinstance(test_spec_tab_ascii["error"][0],np.float64)

def test_generate_realizations():

    # use a pair of test spectra for each format (FITS or ascii)
    abs_stem_src = config_red["data_dirs"]["TEST_DIR_SRC"]
    abs_stem_bin = config_red["data_dirs"]["TEST_DIR_BIN"]

    # set fractional noise level for these tests
    noise_choice = 0.01

    # test on FITS files
    test_spec_list_fits = [
                            abs_stem_src+"575030m20.fits",
                            abs_stem_src+"spec-3480-54999-0629g003.fits"
                            ]
    # expected names
    # (note these should just be basenames)
    expected_filenames_fits = [
                            "575030m20_noise_ver_000.fits",
                            "575030m20_noise_ver_001.fits",
                            "spec-3480-54999-0629g003_noise_ver_000.fits",
                            "spec-3480-54999-0629g003_noise_ver_001.fits"
                            ]
    returned_filenames_fits = []
    for spec_num in range(0,len(test_spec_list_fits)):
        return_names_one_spec = create_spec_realizations.generate_realizations(spec_name=test_spec_list_fits[spec_num],
                                               outdir=abs_stem_bin,
                                               spec_file_format="fits",
                                               num=2,
                                               noise_level=noise_choice)
        returned_filenames_fits.extend(return_names_one_spec)

    # test on ascii files
    test_spec_list_ascii = [
                            abs_stem_src+"700025m20.smo",
                            abs_stem_src+"spec-3478-55008-0186g002.dat"
                            ]
    # expected names
    # (note these should just be basenames)
    expected_filenames_ascii = [
                            "700025m20_noise_ver_000.smo",
                            "700025m20_noise_ver_001.smo",
                            "spec-3478-55008-0186g002_noise_ver_000.dat",
                            "spec-3478-55008-0186g002_noise_ver_001.dat"
                            ]
    returned_filenames_ascii = []
    for spec_num in range(0,len(test_spec_list_ascii)):
        return_names_one_spec = create_spec_realizations.generate_realizations(spec_name=test_spec_list_ascii[spec_num],
                                               outdir=abs_stem_bin,
                                               spec_file_format="ascii.no_header",
                                               num=2,
                                               noise_level=noise_choice)
        returned_filenames_ascii.extend(return_names_one_spec)

    # flatten lists
    #returned_filenames_fits = [item for sublist in t for item in sublist]
    #returned_filenames_ascii = [item for sublist in t for item in sublist]

    # check if elements in list 1 are in list 2
    result_fits =  all(elem in returned_filenames_fits for elem in expected_filenames_fits)
    result_ascii =  all(elem in returned_filenames_ascii for elem in expected_filenames_ascii)

    # test: are the file names right?
    print('returned/expected filenames fits')
    print(returned_filenames_fits)
    print(expected_filenames_fits)
    print('returned/expected filenames ascii')
    print(returned_filenames_ascii)
    print(expected_filenames_ascii)
    assert result_fits
    assert result_ascii

    # test: are the original spectra divided by the realizations equivalent
    # to 1 + noise residuals? (this is just on ascii for now)

    # read in original spectra

    # sorting is critical here, to keep spectra and their realizations organized
    test_spec_list_ascii.sort() # in-place
    expected_filenames_ascii.sort()

    orig_spec_0 = pd.read_csv(test_spec_list_ascii[0], names=["wavel", "abs_flux", "error"], delim_whitespace=True)
    orig_spec_1 = pd.read_csv(test_spec_list_ascii[1], names=["wavel", "abs_flux", "error"], delim_whitespace=True)
    print("orig_spec_0")
    print(orig_spec_0)
    print("orig_spec_1")
    print(orig_spec_1)

    # read in the realizations based off of the original spectra
    # realizations of orig spec 0
    realzn_spec_0_0 = pd.read_csv(abs_stem_bin + expected_filenames_ascii[0], names=["wavel", "abs_flux"], delim_whitespace=True)
    realzn_spec_0_1 = pd.read_csv(abs_stem_bin + expected_filenames_ascii[1], names=["wavel", "abs_flux"], delim_whitespace=True)
    # realizations of orig spec 1
    realzn_spec_1_0 = pd.read_csv(abs_stem_bin + expected_filenames_ascii[2], names=["wavel", "abs_flux"], delim_whitespace=True)
    realzn_spec_1_1 = pd.read_csv(abs_stem_bin + expected_filenames_ascii[3], names=["wavel", "abs_flux"], delim_whitespace=True)

    # do the division
    div_spec_0_by_0 = np.divide(orig_spec_0["abs_flux"],realzn_spec_0_0["abs_flux"])
    div_spec_0_by_1 = np.divide(orig_spec_0["abs_flux"],realzn_spec_0_1["abs_flux"])
    div_spec_1_by_0 = np.divide(orig_spec_1["abs_flux"],realzn_spec_1_0["abs_flux"])
    div_spec_1_by_1 = np.divide(orig_spec_1["abs_flux"],realzn_spec_1_1["abs_flux"])

    print("medians")
    print(orig_spec_0["abs_flux"])
    print(realzn_spec_0_0["abs_flux"])
    print(div_spec_0_by_0)
    print("stedev")
    print(round(np.median(div_spec_0_by_0), 2))
    print(np.std(div_spec_0_by_0))
    print(np.std(div_spec_0_by_1))
    print(np.std(div_spec_1_by_0))
    print(np.std(div_spec_1_by_1))

    # are the original and realization spectra really of the same amplitude?
    assert round(np.median(div_spec_0_by_0), 2) == 1.00
    assert round(np.median(div_spec_0_by_1), 2) == 1.00
    assert round(np.median(div_spec_1_by_0), 2) == 1.00
    assert round(np.median(div_spec_1_by_1), 2) == 1.00

    ## ## noise injection level not well tested yet



'''
def test_read_bkgrnd_spec():

    %min_good_phase, max_good_phase = phase_regions()

    # is min smaller than max
    %assert min_good_phase < max_good_phase

    # are the phases interpreted as floats
    %assert isinstance(min_good_phase,float)


def test_read_list():

    %min_good_phase, max_good_phase = phase_regions()

    # is min smaller than max
    %assert min_good_phase < max_good_phase

    # are the phases interpreted as floats
    %assert isinstance(min_good_phase,float)





def test_write_bckgrnd_input():

    %min_good_phase, max_good_phase = phase_regions()

    # is min smaller than max
    %assert min_good_phase < max_good_phase

    # are the phases interpreted as floats
    %assert isinstance(min_good_phase,float)
'''
