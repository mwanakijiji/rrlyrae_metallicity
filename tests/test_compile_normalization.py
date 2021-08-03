import matplotlib
matplotlib.use('Agg')
from io import StringIO
import sys, os
from unittest.mock import patch, call

current_dir = os.path.dirname(__file__)
target_dir = os.path.abspath(os.path.join(current_dir, "../"))
print(current_dir)
print(target_dir)
sys.path.insert(0, target_dir)

#from rrlyrae_metallicity.rrlyrae_metallicity import modules2.compile_normalization
#from modules2 import compile_normalization
from rrlyrae_metallicity.rrlyrae_metallicity import *
from rrlyrae_metallicity.rrlyrae_metallicity.modules2 import compile_normalization

#mock = MagicMock()

@patch('builtins.print')
def test_compile_bkgrnd(self):
    # does bkgrnd compile?

    #final_list = modules2.compile_normalization.compile_bkgrnd()

    with patch('sys.stdout.readlines()[-1]', new = StringIO()) as fake_out:
        modules2.compile_normalization.compile_bkgrnd()
        self.assertEqual(fake_out.getvalue(), "--------------------------")

    # check that what would be printed is what we expect upon
    # successful compilation
    #assert mock_print.mock_calls == [call("--------------------------")]
