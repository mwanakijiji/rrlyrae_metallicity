import matplotlib
matplotlib.use('Agg')

import sys, os
import configparser

current_dir = os.path.dirname(__file__)
target_dir = os.path.abspath(os.path.join(current_dir, "../"))
sys.path.insert(0, target_dir)

# import more things with changed system path
from modules import *
from modules import create_spec_realizations
from conf import *

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
    spec_name_fits = config_red["TEST_DIR_SRC"] + "575030m20.fits"
    test_spec_tab_fits, test_hdr_fits = read_spec(spec_name=spec_name_fits, format="fits")

    # for FITS data, there should be 3 columns of floats, and a header
    assert np.isfinite(test_hdr_fits) == True
    assert len(test_spec_tab_fits.colnames) == 3 # 3 solumns
    assert (test_spec_tab_fits["wavelength"].info.dtype == np.float64)
    assert (test_spec_tab_fits["flux"].info.dtype == np.float64)
    assert (test_spec_tab_fits["error"].info.dtype == np.float64)

    # ascii format
    spec_name_ascii = config_red["TEST_DIR_SRC"] + "700025m20.smo"
    test_spec_tab_ascii, test_hdr_ascii = read_spec(spec_name=spec_name_ascii, format="ascii.no_header")

    # for ascii data, there should be 3 columns of floats, and NO header
    assert np.isfinite(test_hdr_ascii) == False
    assert len(test_spec_tab_ascii.colnames) == 3 # 3 solumns
    assert (test_spec_tab_ascii["wavelength"].info.dtype == np.float64)
    assert (test_spec_tab_ascii["flux"].info.dtype == np.float64)
    assert (test_spec_tab_ascii["error"].info.dtype == np.float64)

def test_generate_realizations():

    # use a pair of test spectra for each format (FITS or ascii)
    abs_stem_src = config_red["data_dirs"]["TEST_DIR_SRC"]
    abs_stem_bin = config_red["data_dirs"]["TEST_DIR_BIN"]

    # test on FITS files
    test_spec_list_fits = [
                            abs_stem_src+"575030m20.fits",
                            abs_stem_src+"spec-3480-54999-0629g003.fits"
                            ]
    # expected names
    expected_filenames_fits = [
                            abs_stem_bin+"575030m20.fits_000",
                            abs_stem_bin+"575030m20.fits_001",
                            abs_stem_bin+"spec-3480-54999-0629g003.fits_000",
                            abs_stem_bin+"spec-3480-54999-0629g003.fits_001"
                            ]
    for spec_num in range(0,len(test_spec_list_fits)):
        return_filenames_fits = create_spec_realizations.generate_realizations(spec_name=test_spec_list_fits[spec_num],
                                               outdir=abs_stem_bin,
                                               spec_file_format="fits",
                                               num=2,
                                               noise_level=0.01)
    # test on ascii files
    test_spec_list_ascii = [
                            abs_stem_src+"700025m20.smo",
                            abs_stem_src+"spec-3478-55008-0186g002.dat"
                            ]
    # expected names
    expected_filenames_ascii = [
                            abs_stem_bin+"700025m20.smo_000",
                            abs_stem_bin+"700025m20.smo_001",
                            abs_stem_bin+"spec-3478-55008-0186g002.dat_000",
                            abs_stem_bin+"spec-3478-55008-0186g002.dat_001"
                            ]
    for spec_num in range(0,len(test_spec_list_ascii)):
        return_filenames_ascii = generate_realizations(spec_name=test_spec_list_ascii[spec_num],
                                               outdir=abs_stem_bin,
                                               spec_file_format="ascii.no_header",
                                               num=2,
                                               noise_level=0.01)

    # check if elements in list __ are in list __
    result_fits =  all(elem in return_filenames_fits  for elem in expected_filenames_fits)
    result_ascii =  all(elem in return_filenames_fits  for elem in expected_filenames_ascii)

    # originals divided by the realizations should be equivalent
    # to 1 + noise residuals

    # are the phases interpreted as floats
    assert result_fits
    assert result_ascii

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
