# -*- coding: utf-8 -*-
"""
one cloud model           
"""

import numpy as np
import os
import h5py
import pandas as pd
# import plotting tools
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import random
import tensorflow as tf
from tensorflow.keras import layers  # , regularizers
from sklearn.model_selection import train_test_split
# import sklearn
from sklearn.metrics import r2_score
from tensorflow.keras.models import Model
from tensorflow import keras
import sys
import time

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

model_id = '1cloud_dropout'

seconds = int(time.time())  # seconds since epoch
date = time.strftime('%Y-%m-%d_%Hh_%Mm_%Ss', time.localtime(seconds))  # date and time for easier human readability

sys.stdout = open('%s.out' % (model_id), 'w')
sys.stderr = open('%s.err' % (model_id), 'w')
datadir = '.'

with h5py.File('%s/data.h5' % datadir, 'r') as hdf:
    spec1 = np.array(hdf.get('MgII2796'))
    spec2 = np.array(hdf.get('MgII2803'))
    labels = np.array(hdf.get('labels'))

# labels = np.loadtxt('%s/labels.txt' % datadir)

number_of_samples = len(labels)
indices = list(range(number_of_samples))
parameters = ['velocity', 'column density', 'b parameter']

spec = np.concatenate((np.expand_dims(spec1, axis=2), np.expand_dims(spec2, axis=2)), axis=2)
# (number_of_samples, 450, 2)

x_train, x_test, y_train, y_test, indices_train, indices_test = train_test_split(spec, labels, indices, test_size=0.1,
                                                                                 random_state=42)

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

# x_train.shape, x_test.shape, y_train.shape, y_test.shape

num_input = 450
num_classes = 3
epochs = 30
batch_size = 32

inputs = layers.Input(shape=(450, 2))

x = layers.Conv1D(filters=16, kernel_size=3, use_bias=False, padding='same')(inputs)
x = layers.BatchNormalization()(x)
x = layers.ReLU()(x)
x = layers.Conv1D(filters=32, kernel_size=3, use_bias=False, padding='same')(x)
x = layers.BatchNormalization()(x)
x = layers.ReLU()(x)
x = layers.Flatten()(x)
x = layers.Dropout(0.15)(x)
x = layers.Dense(units=1024, use_bias=False)(x)
x = layers.BatchNormalization()(x)
x = layers.ReLU()(x)
x = layers.Dense(units=1024, use_bias=False)(x)
x = layers.BatchNormalization()(x)
x = layers.ReLU()(x)
x = layers.Dense(units=1024, use_bias=False)(x)
x = layers.BatchNormalization()(x)
x = layers.ReLU()(x)
out1 = layers.Dense(units=3, use_bias=False)(x)

model = Model(inputs=[inputs], outputs=out1)

filepath = "./checkpoint-{epoch:02d}-{val_loss:.2f}.h5"
checkpoint = tf.keras.callbacks.ModelCheckpoint(filepath, monitor='val_loss',  save_freq='epoch', mode='auto')
callbacks_list = [checkpoint]

model.compile(loss='mse', optimizer=keras.optimizers.RMSprop(lr=1e-7, rho=0.9))

## Function to plot the data from the model training ##
def plot_history(history):
  plt.ylabel('Loss')
  plt.xlabel('Epoch')
  plt.xticks(range(0, len(history['loss'] + 1)))
  plt.plot(history['loss'], label="training", marker='o')
  plt.plot(history['val_loss'], label="validation", marker='o')
  plt.legend()
  plt.savefig("loss.png")

history = model.fit(x_train, y_train, batch_size=batch_size, epochs=epochs, verbose=1, callbacks=callbacks_list,
          validation_data=(x_test, y_test))

print(pd.DataFrame(history.history))
history = pd.DataFrame(history.history)
plot_history(history)

# keras.utils.vis_utils.plot_model(model, to_file='%s_design.png'%model_id, show_shapes=True, show_layer_names=True)
final_pred = model.predict(x_test, batch_size=batch_size, verbose=1)    # use model.predict to make predictions for
                                                                        # some subset of the data

model.save('./%s.h5' % (model_id))

for i in range(3):
    final_pred[:, i] *= stds_list[i]
    final_pred[:, i] += means_list[i]

final_pred_windex = np.concatenate((indices_test, final_pred), axis=1)

f = open('%s_withheld_set_predictions.txt' % model_id, 'w')
for i in range(len(final_pred_windex[:, 0])):
    for j in range(len(final_pred_windex[0])):
        # print(final_pred_windex[i,j])
        # print(type(final_pred_windex[i,j]))
        f.write('%.6f ' % final_pred_windex[i, j])
    f.write('\n')
f.close()

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
