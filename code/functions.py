# General python functions
# Bryson Stemock
# Last updated:  2020-06-24

import os
import smtplib
from email.mime.text import MIMEText
import pandas as pd
from astropy.coordinates import SkyCoord

def Jname_to_decimals(Jname):
    # Written:      2020-06-24
    # Last updated: 2020-06-24
    """Takes input Jname with the following format: "J001306+000431" and returns the ra and dec in decimal degrees.
    Ex.  >>>ra, dec = Jname_to_decimals("J001306+000431")
            ra = 3.275, dec = 0.07527777777777778"""
    ra_string = Jname[1:3] + "h" + Jname[3:5] + "m" + Jname[5:7] + "s"
    dec_string = Jname[7:10] + "d" + Jname[10:12] + "m" + Jname[12:14] + "s"
    x = SkyCoord(ra_string + " " + dec_string)
    ra, dec = x.ra.degree, x.dec.degree
    return ra, dec


def delete_lines(original_file, line_numbers):
    """Takes a filename and a list of all 0 indexed line numbers to be removed and removes those lines from the file."""
    is_skipped = False
    current_index = 0
    dummy_file = original_file + ".bak"
    # Open original file in read only mode and dummy file in write mode
    with open(original_file, "r") as read_obj, open(dummy_file, "w") as write_obj:
        # Line by line copy data from original file to dummy file
        for line in read_obj:
            # If current line number matches the given line number then skip copying
            if current_index not in line_numbers:
                write_obj.write(line)
            else:
                is_skipped = True
            current_index += 1
    # If any line is skipped then rename dummy file as original file
    if is_skipped:
        os.remove(original_file)
        os.rename(dummy_file, original_file)
    else:
        os.remove(dummy_file)
    return


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
    df = pd.concat((df, df_))
    return pd.Series(df["values"])


def get_atomic(path):
    df = pd.read_table(path, index_col="trans", delim_whitespace=True)
    return df


def wave_to_vel(con, waves, wave_cen, zabs):
    return con["ckms"] * (waves / (1.0 + zabs) - wave_cen) / wave_cen


def vel_to_wave(con, vels, wave_cen, zabs):
    return (1.0 + zabs) * ((vels * wave_cen) / con["ckms"] + wave_cen)


# def send_finished_code_email(scriptname, sent="bstemock@gmail.com", to="bstemock@nmsu.edu", subject="CODE FINISHED"):
#     msg = MIMEText(scriptname + " has finished running.")
#     msg["Subject"] = subject
#     msg["From"] = sent
#     msg["To"] = to
#
#     s = smtplib.SMTP("localhost")
#     s.sendmail(sent, [to], msg.as_string())
#     s.quit()
#     return
