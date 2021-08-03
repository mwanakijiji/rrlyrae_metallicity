import matplotlib
matplotlib.use('Agg')

import sys, os
from unittest.mock import patch, call

current_dir = os.path.dirname(__file__)
target_dir = os.path.abspath(os.path.join(current_dir, "../"))
print(current_dir)
print(target_dir)
sys.path.insert(0, target_dir)

from rrlyrae_metallicity.rrlyrae_metallicity import modules2
from rrlyrae_metallicity.rrlyrae_metallicity import *
from rrlyrae_metallicity.rrlyrae_metallicity.modules2 import *

@patch('builtins.print')
def test_compile_bkgrnd():
    # does bkgrnd compile?

    final_list = compile_bkgrnd(compiled_bkgrnd_file_path_abs_pass=compiled_bkgrnd_file_path_abs,
                                cc_bkgrnd_file_path_abs_pass=cc_bkgrnd_file_path_abs)

    # check that what would be printed is what we expect upon
    # successful compilation
    assert mocked_print.mock_calls == [call("--------------------------")]
