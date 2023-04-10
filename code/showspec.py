import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt

datafile, i = sys.argv[1], int(sys.argv[2])         # System to be plotted

with h5py.File(datafile, "r") as hdf:
    MgII2796 = np.array(hdf.get("MgII2796"))
    MgII2803 = np.array(hdf.get("MgII2803"))
    labels = np.array(hdf.get("labels"))

f2796 = MgII2796[i]      # Pull proper row of flux
f2803 = MgII2803[i]
label = labels[i]
x = np.linspace(-250, 250, len(f2796), endpoint=True)  # Velocities

fig = plt.figure(figsize=(10.0, 10.0))
plt.subplots_adjust(hspace=0.3)
ax1 = fig.add_subplot(2, 1, 2)
ax2 = fig.add_subplot(2, 1, 1)

fig.suptitle("Spectrum #" + str(i) + ", " + str(label))
ax1.set_xlabel("Velocity (km/s)")
ax1.set_ylabel("Normalized Flux")
ax1.title.set_text("MgII2796")
ax1.set_ylim(-0.1, 1.2*np.max(f2796))
ax1.step(x, f2796, color="black", linewidth=0.5, where="mid")
ax1.plot(x, f2796*0., color="blue", linewidth=0.5)
ax1.plot(x, np.ones(len(f2796)), color="blue", linewidth=0.5)

ax2.set_xlabel("Velocity (km/s)")
ax2.set_ylabel("Normalized Flux")
ax2.title.set_text("MgII2803")
ax2.set_ylim(-0.1, 1.2*np.max(f2803))
ax2.step(x, f2803, color="black", linewidth=0.5, where="mid")
ax2.plot(x, f2803*0., color="blue", linewidth=0.5)
ax2.plot(x, np.ones(len(f2803)), color="blue", linewidth=0.5)

plt.show()
