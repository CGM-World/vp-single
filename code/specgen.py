# Generates an input number n of absorption line spectra for a given line list
# last updated: 2023-04-07
# Bryson Stemock and Avery Lee

# Imports
import os
import sys
import h5py
import time
import random
import multiprocessing
import numpy as np
import pandas as pd
from pyDOE import lhs
from progressbar import ProgressBar, SimpleProgress, Bar, ETA

# Set path and import from our other files
sys.path.insert(1, "/fs1/project/cgm_world/code")
from specsynth import Absorber, Instrument, get_ew_spec, get_abs_regions, dblt_checker, cleanspec
from functions import get_atomic, get_constants

    ## SPECGEN ##

# SETUP #

# Start timer so that we know how long the multiprocessed version runs
tic = time.perf_counter()

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

# If a file path is given as an argument, read as hdf5 file and extract parameters
except ValueError:
    try:
        with h5py.File(sys.argv[1], "r") as f:
            v = f["labels"][:, 0]                           # rest-frame velocity
            logN = f["labels"][:, 1]                        # base 10 log of column density
            b = f["labels"][:, 2]                           # Doppler b parameter
            snr = f["noise_log"][:, 0]                      # signal-to-noise ratio
            seed2796 = f["noise_log"][:, 1].astype(int)     # random seeds for noise generation in MgII2796 spectra
            seed2803 = f["noise_log"][:, 2].astype(int)     # random seeds for noise generation in MgII2803 spectra
    except FileNotFoundError:
        sys.exit("ERROR: File not found. Please provide an integer number of spectra to generate or provide " +
                 "the path to an hdf5 data file to recreate spectra using its stored parameters.")

# Declare constants to be used
zabs = 0.000000                                                                             # absorber redshift
vel_range = [-250, 250]                                                                     # velocity window
ATOM = get_atomic("/fs1/project/cgm_world/code/atoms.dat")                                  # atomic data
con = get_constants("/fs1/project/cgm_world/code/const.dek")                                # constants
HIRES = Instrument(con, R=45000, presel=3.0, rdnoise=3.0, slit=1.0, n=3.0, resfac=3.0)      # instrument

# Define lists we will append to in order to create our data frames later
full_cube_labels_builder = []
full_cube_noise_log_builder = []
labels_builder = []
noise_log_builder = []
MgII2796_builder = []
MgII2803_builder = []

# GENERATION #

# Function to generate data. Processes will share this workload
def generate(arr, result, index):

    # Define builder lists
    full_cube_labels_b = []
    full_cube_noise_log_b = []
    labels_b = []
    noise_log_b = []
    MgII2796_b = []
    MgII2803_b = []

    for i in arr:

        # Append to full_cube builder lists
        full_cube_labels_b.append(list(map(float, [str(v[i]), str(logN[i]), str(b[i])])))
        full_cube_noise_log_b.append(list(map(float, [str(snr[i]), seed2796[i], seed2803[i]])))

        # Generate a spectrum for each transition
        MgII2796 = Absorber(con, ATOM, trans="MgII2796", vel_range=vel_range, Inst=HIRES, seed=seed2796[i], snr=snr[i], zabs=zabs, v=[v[i]], logN=[logN[i]], b=[b[i]])
        MgII2803 = Absorber(con, ATOM, trans="MgII2803", vel_range=vel_range, Inst=HIRES, seed=seed2803[i], snr=snr[i], zabs=zabs, v=[v[i]], logN=[logN[i]], b=[b[i]])

        # Check for 5 sigma detections of dominant transition
        ew_spec2796, ew_sig2796 = get_ew_spec(con, MgII2796, HIRES)
        abs_vels2796 = get_abs_regions(MgII2796, ew_spec2796, ew_sig2796, sigma_threshold=5.0)

        # Check for aligned 3 sigma detections in non-dominant transitions
        ew_spec2803, ew_sig2803 = get_ew_spec(con, MgII2803, HIRES)
        abs_flags2803 = get_abs_regions(MgII2803, ew_spec2803, ew_sig2803, sigma_threshold=3.0, dominant=False, region_vels=abs_vels2796)
        abs_vels2796 = dblt_checker(abs_vels2796, abs_flags2803)

        # Throw out systems that aren't detectable above the noise threshold
        if len(abs_vels2796) == 0:
            continue

        # Set pixels outside of absorbing regions to unity
        # MgII2796, MgII2803 = cleanspec(MgII2796, abs_vels2796), cleanspec(MgII2803, abs_vels2796)

        # Append to the other builder lists
        labels_b.append(list(map(float,[str(v[i]), str(logN[i]), str(b[i])])))          # parameters to labels
        noise_log_b.append(list(map(float, [str(snr[i]), seed2796[i], seed2803[i]])))   # snr and seeds to noise_log
        MgII2796_b.append(list(map(float, MgII2796.f_norm)))                            # MgII2796 ransition flux to transition
        MgII2803_b.append(list(map(float, MgII2803.f_norm)))                            # MgII2803 ransition flux to transition

    # Save the results
    result[index] = [full_cube_labels_b, full_cube_noise_log_b, labels_b, noise_log_b, MgII2796_b, MgII2803_b]

# Split the work of generation among a specified number of processes
def run_processes():
    process_count = 56                              # Matches number of cores on Discovery
    processes = [None] * process_count              # Array to hold our processes
    manager = multiprocessing.Manager()             # Helps us share data between processes
    results = manager.dict()                        # Used to store each process's results
    length = (int)(v.shape[0] / process_count)      # Helps us divide the work evenly

    # Create processes to generate even-sized chunks of data
    for i in range(process_count):

        # Most processes generate ith chunk of given length
        if (i != process_count - 1):
            arr = list(range((i*length),((i+1)*length)))
            processes[i] = multiprocessing.Process(target=generate, args=(arr,results,i,))

        # Last process processed differently; it generates remainder of data
        else:
            arr = list(range((i*length), v.shape[0]))
            processes[i] = multiprocessing.Process(target=generate, args=(arr,results,i,))

        # Start process
        processes[i].start()

    # Wait until all processes finish
    for process in processes:
        process.join()

    # Append each process's results to make one list for each attribute
    for i in range(process_count):
        if results[i] is not None:
            full_cube_labels_builder.extend(results[i][0])
            full_cube_noise_log_builder.extend(results[i][1])
            labels_builder.extend(results[i][2])
            noise_log_builder.extend(results[i][3])
            MgII2796_builder.extend(results[i][4])
            MgII2803_builder.extend(results[i][5])

run_processes()

# SAVING RESULTS #

# Convert our builder lists to pandas dataframes
full_cube_labels_df = pd.DataFrame(full_cube_labels_builder)
full_cube_noise_log_df = pd.DataFrame(full_cube_noise_log_builder)
labels_df = pd.DataFrame(labels_builder)
noise_log_df = pd.DataFrame(noise_log_builder)
MgII2796_df = pd.DataFrame(MgII2796_builder)
MgII2803_df = pd.DataFrame(MgII2803_builder)

# Add dataframes to our output file, data.h5
with h5py.File("data.h5", "w") as hdf:
    hdf.create_dataset("data/labels", data=labels_df, compression="gzip", compression_opts=9)
    hdf.create_dataset("data/noise_log", data=noise_log_df, compression="gzip", compression_opts=9)
    hdf.create_dataset("data/MgII2796", data=MgII2796_df, compression="gzip", compression_opts=9)
    hdf.create_dataset("data/MgII2803", data=MgII2803_df, compression="gzip", compression_opts=9)
    hdf.create_dataset("full_cube/labels", data=full_cube_labels_df, compression="gzip", compression_opts=9)
    hdf.create_dataset("full_cube/noise_log", data=full_cube_noise_log_df, compression="gzip", compression_opts=9)

# Stop timer and print time
toc = time.perf_counter()
print("Multiprocessing time: ", (toc - tic))
