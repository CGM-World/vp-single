# -*- coding: utf-8 -*-
"""
one cloud model           
"""

# Imports needed for project
import os
import sys
import time
import h5py
import random
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import tensorflow as tf
import keras_tuner as kt
from tensorflow import keras
from tensorflow.keras import layers
from tensorflow.keras.models import Model
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score
matplotlib.use('Agg')       # Needed for saving plot

# Set seed values
seed_value = 42
os.environ['PYTHONHASHSEED']=str(seed_value)
random.seed(seed_value)
np.random.seed(seed_value)
tf.compat.v1.set_random_seed(seed_value)

# Configure a new global `tensorflow` session
session_conf = tf.compat.v1.ConfigProto(intra_op_parallelism_threads=1, inter_op_parallelism_threads=1)
sess = tf.compat.v1.Session(graph=tf.compat.v1.get_default_graph(), config=session_conf)
tf.compat.v1.keras.backend.set_session(sess)

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

# Set variables to control corresponding file and directory names
model_id = '1cloud_dropout'
sys.stdout = open('%s.out' % (model_id), 'w')
sys.stderr = open('%s.err' % (model_id), 'w')
date = time.strftime('%Y-%m-%d_%Hh_%Mm_%Ss', time.localtime(int(time.time())))  
datadir = '.'

# Load our data
with h5py.File('%s/data.h5' % datadir, 'r') as hdf:
    spec1 = np.array(hdf.get('MgII2796'))
    spec2 = np.array(hdf.get('MgII2803'))
    labels = np.array(hdf.get('labels'))

print(spec1.shape)
print(spec2.shape)
print(labels.shape)

# Prepare data for training
indices = list(range(len(labels))) 
parameters = ['velocity', 'column density', 'b parameter']
spec = np.concatenate((np.expand_dims(spec1, axis=2), np.expand_dims(spec2, axis=2)), axis=2)

x_train, x_test, y_train, y_test, indices_train, indices_test = train_test_split(spec, labels, indices, test_size=0.1, random_state=42)

indices_test = np.expand_dims(indices_test, axis=1)
means_list = []
stds_list = []

for i in range(3):
    tmp_mean = np.mean(y_train[:, i])
    tmp_std = np.std(y_train[:, i])
    print('mean for', parameters[i], ':', tmp_mean)
    print('std for', parameters[i], ':', tmp_std)
    y_train[:, i] -= tmp_mean
    y_train[:, i] /= tmp_std
    means_list.append(tmp_mean)
    stds_list.append(tmp_std)

# Define hyperparameter variables
num_input = 450
num_classes = 3
epochs = 10
batch_size = 32

# Define variables for creating checkpoints
filepath = "./checkpoint-{epoch:02d}-{val_loss:.2f}.h5"
checkpoint = tf.keras.callbacks.ModelCheckpoint(filepath, monitor='val_loss',  save_freq='epoch', mode='auto')
callbacks_list = [checkpoint]

# Function to plot the data from the model training
def plot_history(history):
  plt.ylabel('Loss')
  plt.xlabel('Epoch')
  plt.xticks(np.arange(0, epochs+1, (epochs/10)))
  plt.plot(history['loss'], label="training", marker='o')
  plt.plot(history['val_loss'], label="validation", marker='o')
  plt.legend()
  plt.savefig("loss.png")

# Build our model
def build_model():

    # Tuning parameters
    f = 256
    ks = 3
    f2 = 32
    ks2 = 1
    du = 1024
    du2 = 1024

    model = keras.models.Sequential()

    # Add convolutional layers. Tuning for filters and kernel
    model.add(layers.Conv1D(f,ks,use_bias=False,padding='same'))
    model.add(layers.Activation('relu'))
    model.add(layers.BatchNormalization())
    model.add(layers.Conv1D(f2,ks2,use_bias=False,padding='same'))
    model.add(layers.Activation('relu'))
    model.add(layers.BatchNormalization())
    
    # Max pool and flatten convolutional output
    model.add(layers.MaxPooling1D())
    model.add(layers.Flatten())

    # Add dense layer. Tuning for units
    model.add(layers.Dense(du))
    model.add(layers.Activation('relu'))
    model.add(layers.BatchNormalization())

    # Add dropout layer
    model.add(layers.Dropout(0.5))

    # Add one more dense layer. Tuning for units
    model.add(layers.Dense(du2))
    model.add(layers.Activation('relu'))
    model.add(layers.BatchNormalization())

    # Add the final output layer, compile, and return
    model.add(layers.Dense(3))
    model.compile(optimizer=keras.optimizers.RMSprop(lr=1e-7, rho=0.9), loss="mse", metrics=["accuracy"])
    return model

# Create the model
print("Model")
model = build_model()

# Fit the model with the data and plot its progress
print("Fit")
history = model.fit(x_train, y_train, batch_size=batch_size, epochs=epochs, verbose=1, callbacks=callbacks_list,validation_split=0.2)
history = pd.DataFrame(history.history)
plot_history(history)

# Make predictions
print("Predict")
final_pred = model.predict(x_test, batch_size=batch_size, verbose=1)

for i in range(3):
    final_pred[:, i] *= stds_list[i]
    final_pred[:, i] += means_list[i]

final_pred_windex = np.concatenate((indices_test, final_pred), axis=1)

f = open('%s_withheld_set_predictions.txt' % model_id, 'w')
for i in range(len(final_pred_windex[:, 0])):
    for j in range(len(final_pred_windex[0])):
        f.write('%.6f ' % final_pred_windex[i, j])
    f.write('\n')
f.close()

# Save the model
model.save('./%s.h5' % (model_id))

fig, ax = plt.subplots(1, num_classes, figsize=(20, 4))
i = 0
titles = ["Velocity1", "logN1", "b1"]
for ii in range(num_classes):
    R2 = r2_score(y_test[:, ii], final_pred[:, ii])
    ax[ii].set_title(titles[ii])
    ax[ii].set_ylabel('the predicted values')
    ax[ii].set_xlabel('the true values')
    ax[ii].annotate('R2='+('%0.3f' % R2)+'', xy=(0.05, 0.9), xycoords='axes fraction',
                    bbox=dict(boxstyle="round", fc="w"), size=14)
    ax[ii].plot(y_test[:, ii], final_pred[:, ii], 'bo', mfc='white', alpha=0.5)
    ax[ii].plot(y_test[:, ii], y_test[:, ii], 'k-')  # where the trend should be, plot after everything else
fig.tight_layout()
fig.savefig('%s.png' % (model_id))

# Create new directory and move this iteration's output files there
os.system("mkdir %s" % date)
os.system("mv %s.out %s.err %s.png %s.h5 %s_withheld_set_predictions.txt loss.png checkpoint* error.out result.out %s" % (model_id, model_id, model_id, model_id, model_id, date))

sys.stderr.close()
sys.stderr = sys.__stderr__
