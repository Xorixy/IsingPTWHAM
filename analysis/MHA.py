from logging import exception

import h5py as h5
import matplotlib.pyplot as plt
import numpy as np
import scipy as scp

filename = "../data/manual_spacing.h5"

size = 100*100

time_series_start = 0
time_series_end = 1000

n_e_bins = 300
n_m_bins = 300

tol = 10e-9
def multiple_histogram_analysis(betas, energy_series, observable_series):
    1

def multiple_histogram_analysis(betas, energy_series, observable_series, log_zs):
    1


def find_sim_log_partition_functions(betas, energy_series):
    new_log_zs = np.zeros(len(betas))
    print(new_log_zs)
    for _ in range (5):
        old_log_zs = new_log_zs.copy()
        for b in range(len(betas)):
            print(b)
            new_log_zs[b] = approximate_log_quantity(betas[b], 0, betas, energy_series, old_log_zs)
        print(new_log_zs)

def approximate_log_quantity(beta, log_z, betas, energy_series, log_zs, observable_series = None):
    if observable_series is None:
        observable_series = np.ones(energy_series.shape)
    log_denominator_series = np.zeros(energy_series.shape)
    for [i, series] in enumerate(energy_series):
        for [j, energy] in enumerate(series):
            log_denominator_series[i, j] = log_denominator_weight(beta, energy, betas, len(series), log_zs)
    term_exponents = (np.log(observable_series) - np.log(log_denominator_series)).flatten()
    largest_term = np.max(term_exponents)
    term_exponents -= largest_term
    term_exponents = np.delete(term_exponents, np.argmax(term_exponents))
    sum_log = np.log1p(np.sum(np.exp(term_exponents)))
    return -log_z + largest_term + sum_log
def log_denominator_weight(beta, energy, betas, sizes, log_zs):
    term_exponents = np.log(sizes) - log_zs - betas*energy
    largest_term = np.max(term_exponents)
    term_exponents -= largest_term
    term_exponents = np.delete(term_exponents, np.argmax(term_exponents))
    sum_log = np.log1p(np.sum(np.exp(term_exponents)))
    return beta*energy + largest_term + sum_log


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

find_sim_log_partition_functions(Ks, energy_time_series)
