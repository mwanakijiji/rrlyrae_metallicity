import configparser
from rrlyrae_metallicity.modules2 import *


# check if the directory-making function works
def test_make_dirs():

    # call function to make directories
    make_dirs()

    # do all the directories exist now?
    for vals in config["data_dirs"]:
        abs_path_name = str(config["data_dirs"][vals])
        assert os.path.exists(abs_path_name)


# test if the phase region boundaries are being read in correctly
def test_phase_regions():

    min_good_phase, max_good_phase = phase_regions()

    # is min smaller than max
    assert min_good_phase < max_good_phase

    # are the phases interpreted as floats
    assert isinstance(min_good_phase,float)
