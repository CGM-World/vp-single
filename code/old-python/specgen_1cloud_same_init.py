# generates 1 cloud spectra from premade input files
# Don't forget to adjust the R factor in sysanal.inp before running this code!
# last updated: 2022-12-08
# Bryson Stemock and Avery Lee

import os
import sys
import h5py
import pandas as pd
import numpy as np
import random
from astropy.io import ascii
from progressbar import ProgressBar, SimpleProgress, Bar, ETA

# Set seed values
seed_value = 42
os.environ['PYTHONHASHSEED']=str(seed_value)
random.seed(seed_value)
np.random.seed(seed_value)

rebin = "/fs1/project/cgm_world/code/fortran/rebinspec/rebinspec"
ref2796 = "/fs1/project/cgm_world/code/fortran/input_files/MgII2796.ref"
ref2803 = "/fs1/project/cgm_world/code/fortran/input_files/MgII2803.ref"
specsynth = "/fs1/project/cgm_world/code/fortran/specsynth/specsynth"
sysanal = "/fs1/project/cgm_world/code/fortran/sysanal/sysanal"
pbar = ProgressBar(widgets=[SimpleProgress(), Bar(), ETA()])

# Ensure needed files exist
if not os.path.exists("ions.table") or not os.path.exists("sysanal.inp") or not os.path.exists("zabs.dat"):
    sys.exit("ERROR: You need ions.table, zabs.dat, and sysanal.inp in your working directory to run sysanal.")
if not os.path.exists("init_pars.txt") or not os.path.exists("noise_log.txt"):
    sys.exit("ERROR: You need init_pars.txt and noise_log.txt in your working directory to generate systems.")

# Remind the user to adjust the R factor in sysanal.inp
check_flag = 0
while check_flag == 0:
    check = input("Did you adjust the R factor in sysanal.inp? [y/n]\n")
    if check == "y" or check == "Y":
        check_flag = 1
    elif check == "n" or check == "N":
        sys.exit("Exiting...")
    else:
        print("Invalid input.")

c_km = 2.99792458e5             # speed of light (in km/s) used in sysanal

# # delete residual files from a previous run
# old_run_files = ["labels.txt", "MgII2796data.txt", "MgII2803data.txt"]
# for old_file in old_run_files:
#     if os.path.exists(old_file):
#         os.system("rm " + old_file)

init_pars = np.loadtxt("init_pars.txt")
noise_log = np.loadtxt("noise_log.txt")

# Define lists we will append to in order to create our data frames later
labels_builder = []
noise_log_builder = []
MgII2796_builder = []
MgII2803_builder = []

for i in pbar(range(init_pars.shape[0])):

    # Set the parameters
    v = init_pars[i][0]                         # velocity
    logN = init_pars[i][1]                      # logN
    b = init_pars[i][2]                         # b parameter
    snr = noise_log[i][0]                       # signal-to-noise ratio
    randint_2796 = str(int(noise_log[i][1]))    # MgII2796 noise seed
    randint_2803 = str(int(noise_log[i][2]))    # MgII2803 noise seed

    # Generate specsynth.par
    parfile = open("specsynth.par", "w+")
    parfile.write("\n".join(["C***************** input parameters for specsynth *****************************",
                             str(500), "3.0", "45000.", str(1.0), str(1), str(3.0), str(3), str(snr)]))
    parfile.close()

    # Create MgII2796 input file
    inp = open("specsynth.inp", "w+")
    inp.write("1.000\n" + "\t".join(["MgII2796", str(v), str(logN), str(b)]))
    inp.close()
    os.system(specsynth + " specsynth.inp specsynth.par " + randint_2796)
    os.system("mv specsynth.out MgII2796")
    os.system("rm specsynth.inp specsynth.ticks")

    # Create MgII2803 input file
    inp = open("specsynth.inp", "w+")
    inp.write("1.000\n" + "\t".join(["MgII2803", str(v), str(logN), str(b)]))
    inp.close()
    os.system(specsynth + " specsynth.inp specsynth.par " + randint_2803)
    os.system("mv specsynth.out MgII2803")
    os.system("rm specsynth.inp specsynth.par specsynth.ticks")

    # Center absorption line & shift velocity recorded in labels.txt accordingly
    os.system(sysanal)

    # No detection == exclude it from the training sample
    if not os.path.exists("ew_regions.dat"):
        continue

    # Append our parameters to labels
    labels_builder.append(list(map(float, [str(v), str(logN), str(b)])))

    # Append our snr and random seeds to noise log
    noise_log_builder.append(list(map(float, [str(snr), randint_2796, randint_2803])))

    # rebin and clean spectra for the network
    os.system(rebin + " " + ref2796 + " MgII2796 1 0")
    os.system(rebin + " " + ref2803 + " MgII2803 1 0")
    MgII2796 = ascii.read("MgII2796.rebin.txt")
    MgII2803 = ascii.read("MgII2803.rebin.txt")

    # Append current data row to MgII2796
    MgII2796_builder.append(list(map(float, MgII2796["col3"])))

    # Append current data row to MgII2803
    MgII2803_builder.append(list(map(float, MgII2803["col3"])))

    # Remove unneeded files
    os.system("rm MgII2796 MgII2803 *.rebin.txt *.aod *.ewlim *.ewreg *.ews *.orig *.vmom ew_regions.*")

# Convert our builder lists to pandas dataframes
labels_df = pd.DataFrame(labels_builder)
noise_log_df = pd.DataFrame(noise_log_builder)
MgII2796_df = pd.DataFrame(MgII2796_builder)
MgII2803_df = pd.DataFrame(MgII2803_builder)

# Add dataframes to our output file, data.h5
with h5py.File('data.h5', 'w') as hdf:
    hdf.create_dataset('labels', data=labels_df)
    hdf.create_dataset('noise_log', data=noise_log_df)
    hdf.create_dataset('MgII2796', data=MgII2796_df)
    hdf.create_dataset('MgII2803', data=MgII2803_df)
