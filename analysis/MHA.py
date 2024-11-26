from logging import exception

import h5py as h5
import matplotlib.pyplot as plt
import numpy as np
import scipy as scp

filename = "../data/manual_spacing-3.h5"

tol_sqr = 10e-12
n_iter = 10000


def get_mean_energies_sqr(betas, energy_series, raw_log_weights):
    energies_sqr = np.zeros_like(betas)
    for [k, beta] in enumerate(betas):
        energies_sqr[k] = calculate_obs_avg(beta, energy_series, energy_series*energy_series, raw_log_weights)
    return energies_sqr
def get_mean_energies(betas, energy_series, raw_log_weights):
    energies = np.zeros_like(betas)
    for [k, beta] in enumerate(betas):
        energies[k] = calculate_obs_avg(beta, energy_series, energy_series, raw_log_weights)
    return energies
def plot_heat_capacity(betas, energy_series, raw_log_weights, sample_betas = None):
    energies = get_mean_energies(betas, energy_series, raw_log_weights)
    energies_sqr = get_mean_energies_sqr(betas, energy_series, raw_log_weights)
    plt.plot(betas, betas*betas*(energies_sqr - energies*energies))
    if sample_betas is not None:
        betas = sample_betas
        energies = get_mean_energies(betas, energy_series, raw_log_weights)
        energies_sqr = get_mean_energies_sqr(betas, energy_series, raw_log_weights)
        plt.plot(betas, betas*betas*(energies_sqr - energies*energies), 'o')
    plt.show()

def calculate_obs_avg(beta, energy_series, observable_series, raw_log_weights, log_z = None):
    if log_z is None:
        log_z = calculate_log_z(beta, energy_series, raw_log_weights)
    obs_min = np.min(observable_series)
    log_obs = np.log(observable_series - obs_min + 1)
    logs_obs_mean = log_sum(log_obs - beta*energy_series - raw_log_weights) - log_z
    return np.exp(logs_obs_mean) + obs_min - 1

def find_partition_logs(betas, energy_series):
    new_log_zs = np.ones_like(betas)
    for n in range(n_iter):
        old_log_zs = new_log_zs.copy()
        raw_log_weights = calculate_raw_log_weights(betas, energy_series, old_log_zs)
        for k in range(len(betas)):
            new_log_zs[k] = calculate_log_z(betas[k], energy_series, raw_log_weights)
        print(new_log_zs)
        deltas = np.exp(new_log_zs - old_log_zs) - 1
        delta_sqr = np.sum(deltas**2)
        print(delta_sqr)
        if (delta_sqr <= tol_sqr):
            print(n_iter)
            break

    return new_log_zs
def calculate_log_z(beta, energy_series, raw_log_weights):
    return log_sum((-beta*energy_series - raw_log_weights).flatten())

def calculate_raw_log_weights(betas, energy_series, log_zs):
    raw_log_weights = np.zeros_like(energy_series)
    sizes = np.zeros(len(betas))
    for i in range(len(betas)):
        sizes[i] = len(energy_series[i])
    for [i, series] in enumerate(energy_series):
        for [j, energy] in enumerate(series):
            raw_log_weights[i, j] = log_sum(np.log(sizes) - log_zs - betas*energy)

    return raw_log_weights

def log_sum(log_array):
    maximum = np.max(log_array)
    small_array = np.delete(log_array, np.argmax(log_array))
    small_array -= maximum
    return maximum + np.log1p(np.sum(np.exp(small_array)))


f = h5.File(filename, 'r')

n_processes = len(f.keys())

time_series_length = len(f['0/K_series'])
print(time_series_length)

energy_time_series = np.zeros([n_processes, time_series_length])
magnetization_time_series = np.zeros([n_processes, time_series_length])

Ks = np.array([])


for key in f.keys():
    sim = f[key]
    K_series = sim["K_series"]
    for K in K_series[:]:
        if not Ks.__contains__(K):
            Ks = np.append(Ks, K)

if len(Ks) != n_processes:
    raise exception("Error, number of distinct Ks not equal to number of processes.")

Ks = np.sort(Ks)
print(Ks)

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

z_file = h5.File("manual_spacing_log_zs.h5", 'r')
log_zs = z_file.get('log_zs')
log_zs = log_zs[:]

raw_weights = calculate_raw_log_weights(Ks, energy_time_series, log_zs)

print(calculate_obs_avg(Ks[-1], energy_time_series, energy_time_series, raw_weights))

plot_heat_capacity(np.linspace(Ks[0], Ks[-1], 200), energy_time_series, raw_weights, Ks)


#log_zs = find_partition_logs(Ks, energy_time_series)

#log_zs = find_sim_log_partition_functions(Ks, energy_time_series)
#out = h5.File("manual_spacing_log_zs.h5", 'w')
#dset = out.create_dataset(name='log_zs', data=log_zs)


#print(Ks)
#plt.legend()
#plt.show()