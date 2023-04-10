# General ML python functions
# Bryson Stemock
# Last updated:  2020-06-24

import os
import sys
import fileinput
import numpy as np
import pandas as pd
from random import randrange
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score
from sklearn.metrics import mean_squared_error

rebinspec = "/matrix/cwc/Programs/sarith/rebinspec/rebinspec"      # actual rebinspec code
ref2796 = "/media/MAT-DAT/cwc/ML-VPfit/Rebinspec/MgII2796.ref"
ref2803 = "/media/MAT-DAT/cwc/ML-VPfit/Rebinspec/MgII2803.ref"
specsynth = "/media/MAT-DAT/cwc/ML-VPfit/specsysnth/specsynth"
sysanal = "/matrix/cwc/Programs/sysanal/sysanal"
minfit4 = "/matrix/cwc/Swinburne/Swin-Programs/Programs/minfit/minfit4/minfit"
minfit4S = "/matrix/cwc/Programs/minfit/minfit4S/minfit"


# Standardize labels (to have zero mean and unity std)
def standardize_labels(labels_train, labels_obs):
    m, s = np.zeros(3), np.zeros(3)
    for i in range(3):
        s[i] = np.std(labels_train[:, i])
        m[i] = np.mean(labels_train[:, i])
        labels_train[:, i] -= m[i]
        labels_train[:, i] /= s[i]
        labels_obs[:, i] -= m[i]    # we need to use the mean/std of the training set because obviously
        labels_obs[:, i] /= s[i]    # this is what the network sees, otherwise an offest is expected.
    return labels_train, labels_obs, m, s


def standardize_spectra(specs_train):
    for i in range(specs_train.shape[0]):
        specs_train[i] = (specs_train[i] - np.mean(specs_train[i])) / np.std(specs_train[i])
        specs_train[i] = (specs_train[i] - min(specs_train[i])) / (max(specs_train[i]) - min(specs_train[i]))
    return specs_train


def find_nearest(array, values):    # function to find index of closest values in array
    diff = np.sum((array - values)**2., axis=1)
    idx = diff.argmin()
    return idx


# Plot two similar spectra from the training and observed samples. For a particular index
# in the training sample, set idx equal to that index. For a random index, leave idx=None.
def plot_nearest(labels_train, labels_obs, specs_train, specs_obs, idx=None):
    if idx is None:
        any_ind = np.random.randint(0, specs_obs.shape[0])  # pick up a random index from the observed spectra
    else:
        any_ind = idx
    print("spectrum index", any_ind)
    near_ind = find_nearest(labels_train, labels_obs[any_ind])
    train_any = labels_train[near_ind]
    obs_any = labels_obs[any_ind]

    num_pix = specs_train.shape[1]
    fig, ax = plt.subplots(2, sharex="all", figsize=(23, 10))
    ax[0].step(range(num_pix), specs_train[near_ind])
    ax[0].set_title("train_param= (" + str("%.2f" % train_any[0]) + ", "+str("%.2f" % train_any[1]) +
                    ", " + str("%.2f" % train_any[2]) + ") ", fontsize=30)
    ax[1].step(range(num_pix),  specs_obs[any_ind])
    ax[1].set_title("obs_param= (" + str("%.2f" % obs_any[0]) + ", "+str("%.2f" % obs_any[1]) +
                    ", " + str("%.2f" % obs_any[2]) + ") ", fontsize=30)
    for i in range(2):
        ax[i].tick_params(labelsize=16)
        ax[i].set_ylim([0, 1.2])
    plt.show()
    return


# Print r2 score and mean squared error and plot prediction results
def display_results(model, labels_test, specs_test, m, s):
    pred_test = model.predict(specs_test)
    fig, ax = plt.subplots(1, 3, figsize=(16, 4))
    for i in range(3):
        tr, pr = labels_test[:, i]*s[i] + m[i], pred_test[:, i]*s[i] + m[i]
        ax[i].plot(tr, pr, 'o')
        ax[i].plot(tr, tr, 'k-')
        ax[i].set_xlim([tr.min(), tr.max()])
        ax[i].set_ylim([tr.min(), tr.max()])
        ax[i].tick_params(labelsize=14)
        print(" param %d %.4f %.4f" % (i, r2_score(tr, pr), mean_squared_error(tr, pr)))
    plt.tight_layout()
    plt.show()
    return


def makespec(v, logN, b, R=45000, p=3.0, snr=0, randint2796=-int(randrange(1, 100000, 1)),
             randint2803=-int(randrange(1, 100000, 1)), zabs="1.000000", rebin=False, clean=False, nuke=False):
    # all inputs are strings
    with open("zabs.dat", "w+") as zfile:
        zfile.write(str(zabs))

    parfile = open("specsynth.par", "w+")  # generate specsynth.par
    parfile.write("\n".join(["C***************** input parameters for specsynth *****************************",
                             str(500), str(p), str(R), str(1.0), str(1), str(3.0), str(3), str(snr)]))
    parfile.close()

    inp = open("specsynth.inp", "w+")  # create MgII2796 input file
    inp.write("1.000\n" + "\t".join(["MgII2796", str(v), str(logN), str(b)]))
    inp.close()

    os.system(specsynth + " specsynth.inp specsynth.par " + str(randint2796))  # run specsynth
    os.system("mv specsynth.out MgII2796")
    os.system("rm specsynth.inp specsynth.ticks")

    inp = open("specsynth.inp", "w+")  # create MgII2803 input file
    inp.write("1.000\n" + "\t".join(["MgII2803", str(v), str(logN), str(b)]))
    inp.close()

    os.system(specsynth + " specsynth.inp specsynth.par " + str(randint2803))  # run specsynth
    os.system("mv specsynth.out MgII2803")
    os.system("rm specsynth.inp specsynth.par specsynth.ticks")

    if clean:
        rebin_binaries = "1 0"
    elif not clean:
        rebin_binaries = "0 0"
    else:
        rebin_binaries = ""

    if rebin:
        os.system("/matrix/cwc/Programs/sysanal/sysanal")                  # run sysanal
        os.system(rebinspec + " " + ref2796 + " MgII2796 " + rebin_binaries)        # rebin spectra w/o cleaning
        os.system(rebinspec + " " + ref2803 + " MgII2803 " + rebin_binaries)
        os.system("mv MgII2796.rebin.txt MgII2796")                             # rename rebinned spectra
        os.system("mv MgII2803.rebin.txt MgII2803")

    if nuke and rebin:
        os.system("rm MgII2796.* MgII2803.* ew_regions.* sysanal.aod sysanal.ews sysanal.vmom zabs.dat")
    return


def fill_lines(output_file_, max_clouds_):
    max_pars_ = max_clouds_ * 6
    for line_ in fileinput.input(output_file_, inplace=True):
        line_ = line_.rstrip('\r\n')
        fit_vals_ = [x.strip() for x in line_.split("\t")]
        while len(fit_vals_) < max_pars_:
            fit_vals_.append("0.0")
        print("\t".join(fit_vals_))
    return


def find_min_max_pixel(index, lines):   # takes an index and f.readlines() as lines from some file f = open("ex.txt")
    line = lines[index]
    flux = [x.strip() for x in line.split(" ") if x != ""]
    i_min, i_max = 0, 0
    for i_min in range(len(flux)):
        if flux[i_min] != "1.0":
            break
    for i_max in reversed(range(len(flux))):
        if flux[i_max] != "1.0":
            break
    return i_min, i_max


def read_vpfit_par_file_1cloud(filename, minfit_flag=False, no_index_flag=False, round=False, delimiter="\t"):
    v, logN, b = [], [], []
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            pars = [x.strip() for x in line.split(delimiter)]
            if not minfit_flag and not no_index_flag:
                v.append(float(pars[1])), logN.append(float(pars[2])), b.append(float(pars[3]))
            elif minfit_flag and not no_index_flag:
                v.append(float(pars[1])), logN.append(float(pars[3])), b.append(float(pars[5]))
            elif no_index_flag and not minfit_flag:
                v.append(float(pars[0])), logN.append(float(pars[1])), b.append(float(pars[2]))
            else:
                sys.exit("ERROR: minfit_flag and no_index_flag must be boolean and both can't be true simultaneously.")
    v, logN, b = np.asarray(v), np.asarray(logN), np.asarray(b)
    if not round:
        return v, logN, b
    else:
        v, logN, b = np.round(v, int(round)), np.round(logN, int(round)), np.round(b, int(round))
        return v, logN, b


def read_cloudy_par_file_phase1(filename):
    logZ, lognH, b_opt, logT, f_HI, logNHI = [], [], [], [], [], []
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            pars = [x.strip() for x in line.split("\t")]
            logZ.append(float(pars[0])), lognH.append(float(pars[1])), b_opt.append(float(pars[2]))
            logT.append(float(pars[3])), f_HI.append(float(pars[4])), logNHI.append(float(pars[5]))
    return (np.asarray(logZ), np.asarray(lognH), np.asarray(b_opt),
            np.asarray(logT), np.asarray(f_HI), np.asarray(logNHI))


def read_cloudy_par_file_phase2(filename):
    logZ, lognH, b_opt, logT, f_HI, logNHI = [], [], [], [], [], []
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            pars = [x.strip() for x in line.split("\t")]
            logZ.append(float(pars[6])), lognH.append(float(pars[7])), b_opt.append(float(pars[8]))
            logT.append(float(pars[9])), f_HI.append(float(pars[10])), logNHI.append(float(pars[11]))
    return (np.asarray(logZ), np.asarray(lognH), np.asarray(b_opt),
            np.asarray(logT), np.asarray(f_HI), np.asarray(logNHI))


def get_constants(path):
    df = pd.read_table(path, index_col="constants", comment="#", delim_whitespace=True)
    hbar = 0.5 * df.loc["h"] / df.loc["pi"]
    sig = df.loc["pi"] * df.loc["pi"] * df.loc["kerg"] * df.loc["kerg"] * df.loc["kerg"] * df.loc["kerg"] /\
        (60.0 * hbar * hbar * hbar * df.loc["c"] * df.loc["c"])
    constants = ["deg2rad", "rad2deg", "hbar", "rbohr", "fine", "sig", "asol", "weinlam", "weinfre", "pc"]
    values = [
        df.loc["pi"] / 180.0,
        180.0 / df.loc["pi"],
        hbar,
        hbar * hbar / (df.loc["me"] * df.loc["e"] * df.loc["e"]),
        df.loc["e"] * df.loc["e"] / (hbar * df.loc["c"]),
        sig,
        4.0 * sig / df.loc["c"],
        df.loc["h"] * df.loc["c"] / (df.loc["kerg"] * 4.965114232),
        2.821439372 * df.loc["kerg"] / df.loc["h"],
        3.261633 * df.loc["lyr"]
    ]
    df_ = pd.DataFrame(data=values, index=constants, columns=df.columns)
    df = df.append(df_)
    return pd.Series(df["values"])


def get_atomic(path):
    df = pd.read_table(path, index_col="trans", delim_whitespace=True)
    return df
