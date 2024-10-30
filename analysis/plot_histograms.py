import h5py as h5
import matplotlib.pyplot as plt
import numpy as np
import scipy as scp

f = h5.File("../data/testos.h5", 'r')

fig, ax = plt.subplots()
fig2, ax2 = plt.subplots()


def autocorrelation(array, n_t):
    mean = np.mean(array)
    mean_sqr = np.mean(array*array)



for key in f.keys():
    sim = f[key]
    hist = sim["histogram"]
    energy = sim["energy"]
    ax.stairs(hist, label = key)
    ax2.plot(energy[:], label = key)



ax.legend()
plt.show()