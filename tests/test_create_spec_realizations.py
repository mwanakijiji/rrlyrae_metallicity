import matplotlib
matplotlib.use('Agg')

import sys, os

current_dir = os.path.dirname(__file__)
target_dir = os.path.abspath(os.path.join(current_dir, "../"))
print(current_dir)
print(target_dir)
sys.path.insert(0, target_dir)

from rrlyrae_metallicity.rrlyrae_metallicity import modules2
from rrlyrae_metallicity.rrlyrae_metallicity import *
from rrlyrae_metallicity.rrlyrae_metallicity.modules2 import *

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

def test_phase_regions():

    min_good_phase, max_good_phase = phase_regions()

    # is min smaller than max
    assert min_good_phase < max_good_phase

    # are the phases interpreted as floats
    assert isinstance(min_good_phase,float)

#
def test_create_norm_spec():

    final_list = create_norm_spec(name_list, bkgrnd_output_dir, final_dir)

    # is min smaller than max
    assert 1 < 2


def test_generate_realizations():

    %min_good_phase, max_good_phase = phase_regions()

    # is min smaller than max
    %assert min_good_phase < max_good_phase

    # are the phases interpreted as floats
    %assert isinstance(min_good_phase,float)


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


def test_read_spec():

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
