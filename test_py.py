import subprocess
import shlex
from subprocess import call
from subprocess import Popen

normzn_compile = shlex.split("g++ -o bkgrnd bkgrnd.cc")
test = subprocess.Popen(normzn_compile)

normzn_run = shlex.split("./bkgrnd --smooth 22 input_file")
test2 = subprocess.Popen(normzn_run, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

print('Finished')
