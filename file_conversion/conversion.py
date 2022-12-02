# conversion.py
# Converts text files to HDF5 format
# Avery Lee: averyl36@yahoo.com
# November 15, 2022

# The MgII data files contain 958112 rows with 450 numbers each

# Import statements
import h5py
import pandas as pd

# SPECIFY YOUR FILE HERE
datadir = '.'
fi = "test"

# Builds strings to represent both versions of your file
txt = datadir + "/" + fi + ".txt"
h5 = datadir + "/" + fi + ".h5"

# Convert file to dataframe
data = []
with open(txt) as f:
    for line in f.readlines():
        data.append(list(map(float, line.split())))
    f.close()
df = pd.DataFrame(data)
print(df.shape)
print(max(df))

# Add to HDF file
with h5py.File(h5, 'w') as hdf:
    hdf.create_dataset('dataset1', data=df)
