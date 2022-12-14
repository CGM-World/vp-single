# generates 1 cloud spectra from premade input files
# Don't forget to adjust the R factor in sysanal.inp before running this code!
# last updated:  2022-03-02
# Bryson Stemock

import os
import sys
import numpy as np
from astropy.io import ascii
from progressbar import ProgressBar, SimpleProgress, Bar, ETA

rebin = "/fs1/project/cgm_world/code/fortran/rebinspec/rebinspec"
ref2796 = "/fs1/project/cgm_world/code/fortran/input_files/MgII2796.ref"
ref2803 = "/fs1/project/cgm_world/code/fortran/input_files/MgII2803.ref"
specsynth = "/fs1/project/cgm_world/code/fortran/specsynth/specsynth"
sysanal = "/fs1/project/cgm_world/code/fortran/sysanal/sysanal"
pbar = ProgressBar(widgets=[SimpleProgress(), Bar(), ETA()])

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

for i in pbar(range(init_pars.shape[0])):
    v = init_pars[i][0]                         # velocity
    logN = init_pars[i][1]                      # logN
    b = init_pars[i][2]                         # b parameter
    snr = noise_log[i][0]                       # signal-to-noise ratio
    randint_2796 = str(int(noise_log[i][1]))    # MgII2796 noise seed
    randint_2803 = str(int(noise_log[i][2]))    # MgII2803 noise seed

    parfile = open("specsynth.par", "w+")  # generate specsynth.par
    parfile.write("\n".join(["C***************** input parameters for specsynth *****************************",
                             str(500), "3.0", "45000.", str(1.0), str(1), str(3.0), str(3), str(snr)]))
    parfile.close()

    inp = open("specsynth.inp", "w+")  # create MgII2796 input file
    inp.write("1.000\n" + "\t".join(["MgII2796", str(v), str(logN), str(b)]))
    inp.close()

    os.system(specsynth + " specsynth.inp specsynth.par " + randint_2796)  # run specsynth

    os.system("mv specsynth.out MgII2796")
    os.system("rm specsynth.inp specsynth.ticks")

    inp = open("specsynth.inp", "w+")  # create MgII2803 input file
    inp.write("1.000\n" + "\t".join(["MgII2803", str(v), str(logN), str(b)]))
    inp.close()

    os.system(specsynth + " specsynth.inp specsynth.par " + randint_2803)  # run specsynth

    os.system("mv specsynth.out MgII2803")
    os.system("rm specsynth.inp specsynth.par specsynth.ticks")

    # center the absorption line and shift the velocity recorded in labels.txt accordingly
    os.system(sysanal)

    if not os.path.exists("ew_regions.dat"):    # no detection == exclude it from the training sample
        continue

    pars = [str(v), str(logN), str(b)]          # give the machine v, logN, b

    # record initial parameters
    ipfile = open("init_pars_.txt", "a+")       # append to label file
    ipfile.write("\t".join(pars) + "\n")
    ipfile.close()

    # record system parameters
    labelfile = open("labels.txt", "a+")        # append to label file
    labelfile.write("\t".join(pars) + "\n")
    labelfile.close()

    noisefile = open("noise_log_.txt", "a+")    # append snr, 2796 seed, and 2803 seed to noise log file
    noisefile.write("\t".join([str(snr), randint_2796, randint_2803]) + "\n")
    noisefile.close()

    # rebin and clean spectra for the network
    os.system(rebin + " " + ref2796 + " MgII2796 1 0")     # rebin spectra
    os.system(rebin + " " + ref2803 + " MgII2803 1 0")

    MgII2796, MgII2803 = ascii.read("MgII2796.rebin.txt"), ascii.read("MgII2803.rebin.txt")

    flux1 = " ".join(MgII2796["col3"].astype("str"))
    datafile1 = open("MgII2796data.txt", "a+")  # Append to MgII2796 data file
    datafile1.write(flux1 + "\n")
    datafile1.close()

    flux2 = " ".join(MgII2803["col3"].astype("str"))
    datafile2 = open("MgII2803data.txt", "a+")  # Append to MgII2803 data file
    datafile2.write(flux2 + "\n")
    datafile2.close()

    os.system("rm MgII2796 MgII2803 *.rebin.txt *.aod *.ewlim *.ewreg *.ews *.orig *.vmom ew_regions.*")

os.system("mv labels.txt labels_ex.txt")
os.system("mv MgII2796data.txt MgII2796data_ex.txt")
os.system("mv MgII2803data.txt MgII2803data_ex.txt")
os.system("mv init_pars_.txt init_pars_ex.txt")
os.system("mv noise_log_.txt noise_log_ex.txt")
