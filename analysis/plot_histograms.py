from logging import exception

import h5py as h5
import matplotlib.pyplot as plt
import numpy as np
import scipy as scp



filename = "../data/big-swap.h5"

size = 100*100

time_series_start = 0
time_series_end = 1000

n_e_bins = 300
n_m_bins = 300


f = h5.File(filename, 'r')

n_processes = len(f.keys())

time_series_length = len(f['0/K_series'])
print(time_series_length)

energy_time_series = np.zeros([n_processes, time_series_length])
magnetization_time_series = np.zeros([n_processes, time_series_length])

Ks = np.array([])

k_fig, k_ax = plt.subplots()
e_fig, e_ax = plt.subplots()
m_fig, m_ax = plt.subplots()
e_hist_fig, e_hist_ax = plt.subplots()
m_hist_fig, m_hist_ax = plt.subplots()


for key in f.keys():
    sim = f[key]
    K_series = sim["K_series"]
    k_ax.plot(K_series[time_series_start:time_series_end], label = key)
    for K in K_series[:]:
        if not Ks.__contains__(K):
            Ks = np.append(Ks, K)

print(Ks)
if len(Ks) != n_processes:
    raise exception("Error, number of distinct Ks not equal to number of processes.")

Ks = np.sort(Ks)

for ik in range(n_processes):
    K = Ks[ik]
    for key in f.keys():
        sim = f[key]
        K_series = sim["K_series"][:]
        E_series = sim["energy_series"][:]
        M_series = sim["magnetization_series"][:]
        Kmask = np.where(K_series == K)

        energy_time_series[ik][Kmask] = E_series[Kmask]
        magnetization_time_series[ik][Kmask] = M_series[Kmask]

magnetization_time_series = magnetization_time_series * 1.0/size

max_energy = np.max(energy_time_series)
min_energy = np.min(energy_time_series)
max_mag = np.max(magnetization_time_series)
min_mag = np.min(magnetization_time_series)




for ik in range(n_processes):
    K = Ks[ik]
    #e_ax.plot(energy_time_series[ik][time_series_start:time_series_end], label = Ks[ik])
    #m_ax.plot(magnetization_time_series[ik][time_series_start:time_series_end], label = Ks[ik])
    #e_hist = np.histogram(energy_time_series[ik], bins=n_e_bins, range=(min_energy, max_energy))
    #m_hist = np.histogram(magnetization_time_series[ik], bins=n_m_bins, range=(min_mag, max_mag))
    #e_hist_ax.stairs(*e_hist, label = Ks[ik])
    #m_hist_ax.stairs(*m_hist, label = Ks[ik])



k_ax.set_title('K diffusion of different processes')
k_ax.set_xlabel('Simulation time')
k_ax.set_ylabel('K')
k_ax.legend()

e_ax.set_title('Energy at different K')
e_ax.set_xlabel('Simulation time')
e_ax.set_ylabel('Energy')
e_ax.legend()

m_ax.set_title('Magnetization at different K')
m_ax.set_xlabel('Simulation time')
m_ax.set_ylabel('Magnetization')
m_ax.legend()

e_hist_ax.set_title('Energy histograms for different K')
e_hist_ax.set_xlabel('Energy')
e_hist_ax.set_ylabel('Number of occurrences')
e_hist_ax.legend()

m_hist_ax.set_title('Magnetization histograms for different K')
m_hist_ax.set_xlabel('Magnetization')
m_hist_ax.set_ylabel('Number of occurrences')
m_hist_ax.legend()

#plt.show()






