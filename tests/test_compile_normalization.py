import matplotlib
matplotlib.use('Agg')
#from io import StringIO
import sys, os
#from unittest.mock import patch, call

current_dir = os.path.dirname(__file__)
target_dir = os.path.abspath(os.path.join(current_dir, "../"))
print(current_dir)
print(target_dir)
sys.path.insert(0, target_dir)

#from rrlyrae_metallicity.rrlyrae_metallicity import modules2.compile_normalization
#from modules2 import compile_normalization
from rrlyrae_metallicity.rrlyrae_metallicity import *
from rrlyrae_metallicity.rrlyrae_metallicity.modules2 import compile_normalization
from rrlyrae_metallicity.rrlyrae_metallicity.modules2 import *

#mock = MagicMock()

#@patch('builtins.print')
def test_compile_bkgrnd():
    # does bkgrnd compile?

    import os
    cwd = os.getcwd()
    print(cwd)
    root_dir_test = "/home/runner/work/rrlyrae_metallicity/rrlyrae_metallicity/"
    #root_dir_test = $GITHUB_WORKSPACE
    compile_status = modules2.compile_normalization.compile_bkgrnd(
        compiled_bkgrnd_file_path_abs_pass = root_dir_test+"src/"+os.path.basename(compiled_bkgrnd_file_path_abs),
        cc_bkgrnd_file_path_abs_pass = "src/bkgrnd.cc"
        )

    #with patch('sys.stdout.readlines()[-1]', new = StringIO()) as fake_out:
#        modules2.compile_normalization.compile_bkgrnd()
    #    self.assertEqual(fake_out.getvalue(), "--------------------------")

    # check that compilation is successful via boolean value
    assert compile_status
