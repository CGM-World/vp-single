# generates latin hypercube initial parameter files to be used to generate systems with varied pixel resolution elements
# last updated:  2022-02-07
# Bryson Stemock

from pyDOE import lhs
import numpy as np
import os
import sys
import random
from progressbar import ProgressBar, SimpleProgress, Bar, ETA

# Set seed values
seed_value = 42
os.environ['PYTHONHASHSEED']=str(seed_value)
random.seed(seed_value)
np.random.seed(seed_value)

pbar = ProgressBar(widgets=[SimpleProgress(), Bar(), ETA()])

try:
    n = int(sys.argv[1])    # number of spectra to simulate
except ValueError:
    sys.exit("ERROR: Please use an integer number of spectra to simulate.")

params = lhs(4, n)          # Latin Hypercube Sampling of parameter space

# delete residual files from a previous run
old_run_files = ["noise_log.txt", "init_pars.txt"]
for old_file in old_run_files:
    if os.path.exists(old_file):
        os.system("rm " + old_file)

for i in pbar(range(len(params))):
    snr = round(params[i][0] * 81. + 9., 3)     # signal-to-noise ratio
    logN = round(params[i][1] * 2.9 + 11.2, 7)  # logN
    b = round(params[i][2] * 9. + 1., 7)        # b parameter
    v = round(params[i][3] * 10. - 5., 7)       # velocity
    randint_2796 = str(-int(random.randrange(1, 100000, 1)))
    randint_2803 = str(-int(random.randrange(1, 100000, 1)))

    pars = [str(v), str(logN), str(b)]          # give the machine v, logN, b

    # record initial parameters
    ipfile = open("init_pars.txt", "a+")        # append to label file
    ipfile.write("\t".join(pars) + "\n")
    ipfile.close()

    noisefile = open("noise_log.txt", "a+")         # append snr, 2796 seed, and 2803 seed to noise log file
    noisefile.write("\t".join([str(snr), randint_2796, randint_2803]) + "\n")
    noisefile.close()
