# generates an input number n of absorption line spectra for a given line list
# last updated: 2022-02-09
# Bryson Stemock and Avery Lee

import os
import sys
import h5py
import random
import numpy as np
import pandas as pd
from pyDOE import lhs
from progressbar import ProgressBar, SimpleProgress, Bar, ETA

sys.path.insert(1, "/fs1/project/cgm_world/code")

from specsynth import Absorber, Instrument, get_ew_spec, get_abs_regions, dblt_checker, cleanspec
from functions import get_atomic, get_constants

# Set seed values
seed_value = 42
os.environ['PYTHONHASHSEED'] = str(seed_value)
random.seed(seed_value)
np.random.seed(seed_value)

# Generate Latin Hypercube if an integer is given in the command line
try:
    n = int(sys.argv[1])                                        # number of spectra to simulate
    pars = lhs(4, n)                                            # Latin Hypercube Sampling of parameter space
    v = np.around(pars[:, 0] * 10.0 - 5.0, 7)                   # rest-frame velocity
    logN = np.around(pars[:, 1] * 2.9 + 11.2, 7)                # base 10 log of column density
    b = np.around(pars[:, 2] * 9.0 + 1.0, 7)                    # Doppler b parameter
    snr = np.around(pars[:, 3] * 81.0 + 9.0, 3)                 # signal-to-noise ratio
    seed2796 = np.random.randint(0, 4294967295, size=v.size)    # random seeds for noise generation in MgII2796 spectra
    seed2803 = np.random.randint(0, 4294967295, size=v.size)    # random seeds for noise generation in MgII2803 spectra
except ValueError:
    # Read in hdf5 file and extract parameters if a path is given
    try:
        with h5py.File(sys.argv[1], "r") as f:
            v = f["labels"][:, 0]                               # rest-frame velocity
            logN = f["labels"][:, 1]                            # base 10 log of column density
            b = f["labels"][:, 2]                               # Doppler b parameter
            snr = f["noise_log"][:, 0]                          # signal-to-noise ratio
            seed2796 = f["noise_log"][:, 1].astype(int)         # random seeds for noise generation in MgII2796 spectra
            seed2803 = f["noise_log"][:, 2].astype(int)         # random seeds for noise generation in MgII2803 spectra
    except FileNotFoundError:
        sys.exit("ERROR: File not found. Please provide an integer number of spectra to generate or provide " +
                 "the path to an hdf5 data file to recreate spectra using its stored parameters.")

# General setup
pbar = ProgressBar(widgets=[SimpleProgress(), Bar(), ETA()])                            # progress bar
zabs = 0.000000                                                                         # absorber redshift
vel_range = [-250, 250]                                                                 # velocity window
ATOM = get_atomic("/fs1/project/cgm_world/code/atoms.dat")                              # atomic data
con = get_constants("/fs1/project/cgm_world/code/const.dek")                            # constants
HIRES = Instrument(con, R=45000, presel=3.0, rdnoise=3.0, slit=1.0, n=3.0, resfac=3.0)  # instrument

# Define lists we will append to in order to create our data frames later
labels_builder = []
noise_log_builder = []
MgII2796_builder = []
MgII2803_builder = []

for i in pbar(range(v.shape[0])):

    # Generate a spectrum for each transition
    MgII2796 = Absorber(con, ATOM, trans="MgII2796", vel_range=vel_range, Inst=HIRES, seed=seed2796[i],
                        snr=snr[i], zabs=zabs, v=[v[i]], logN=[logN[i]], b=[b[i]])
    MgII2803 = Absorber(con, ATOM, trans="MgII2803", vel_range=vel_range, Inst=HIRES, seed=seed2803[i],
                        snr=snr[i], zabs=zabs, v=[v[i]], logN=[logN[i]], b=[b[i]])

    # Check for 5 sigma detections of dominant transition
    ew_spec2796, ew_sig2796 = get_ew_spec(con, MgII2796, HIRES)
    abs_vels2796 = get_abs_regions(MgII2796, ew_spec2796, ew_sig2796, sigma_threshold=5.0)
    if len(abs_vels2796) == 0:
        continue

    # Check for aligned 3 sigma detections in non-dominant transitions
    ew_spec2803, ew_sig2803 = get_ew_spec(con, MgII2803, HIRES)
    abs_flags2803 = get_abs_regions(MgII2803, ew_spec2803, ew_sig2803, sigma_threshold=3.0,
                                    dominant=False, region_vels=abs_vels2796)
    abs_vels2796 = dblt_checker(abs_vels2796, abs_flags2803)
    if len(abs_vels2796) == 0:
        continue

    # Set pixels outside of absorbing regions to unity
    MgII2796, MgII2803 = cleanspec(MgII2796, abs_vels2796), cleanspec(MgII2803, abs_vels2796)

    # Append to builder lists
    labels_builder.append(list(map(float, [str(v[i]), str(logN[i]), str(b[i])])))           # parameters to labels
    noise_log_builder.append(list(map(float, [str(snr[i]), seed2796[i], seed2803[i]])))     # snr, seeds to noise_log
    MgII2796_builder.append(list(map(float, MgII2796.f_norm)))                              # MgII2796 flux to MgII2796
    MgII2803_builder.append(list(map(float, MgII2803.f_norm)))                              # MgII2803 flux to MgII2803

# Convert our builder lists to pandas dataframes
labels_df = pd.DataFrame(labels_builder)
noise_log_df = pd.DataFrame(noise_log_builder)
MgII2796_df = pd.DataFrame(MgII2796_builder)
MgII2803_df = pd.DataFrame(MgII2803_builder)

# Add dataframes to our output file, data.h5
with h5py.File("data.h5", "w") as hdf:
    hdf.create_dataset("labels", data=labels_df, compression="gzip", compression_opts=9)
    hdf.create_dataset("noise_log", data=noise_log_df, compression="gzip", compression_opts=9)
    hdf.create_dataset("MgII2796", data=MgII2796_df, compression="gzip", compression_opts=9)
    hdf.create_dataset("MgII2803", data=MgII2803_df, compression="gzip", compression_opts=9)
