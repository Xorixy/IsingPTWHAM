from logging import exception

import h5py as h5
import matplotlib.pyplot as plt
import numpy as np
import scipy as scp

filename = "../data/manual_spacing-3.h5"

#size = 100*100

time_series_start = 0
time_series_end = 1000

n_e_bins = 300
n_m_bins = 300



def mean_energy_function(betas, energy_series, log_zs, n_points):
    raw_weights = raw_log_weights(betas, energy_series, log_zs)
    sample_betas = np.linspace(betas[0], betas[-1], n_points)
    sample_means = np.zeros(len(sample_betas))
    for k in range(len(sample_betas)):
        log_z = multiple_histogram_analysis(energy_series, sample_betas[k], 0, raw_weights)
        sample_means[k] = multiple_histogram_analysis(energy_series, sample_betas[k], log_z, raw_weights, energy_series)
        sample_means[k] = log_z
    plt.plot(sample_betas, sample_means)
    plt.plot(betas, log_zs, 'o')
    plt.show()



def multiple_histogram_analysis(energy_series, beta, log_z, raw_log_weights, observable_series = None):
    if observable_series is None:
        observable_series = np.ones(energy_series.shape)
    log_observable_series = np.log(observable_series)
    inv_log_weights = beta * energy_series + raw_log_weights
    term_exponents = (log_observable_series - inv_log_weights).flatten()
    largest_term = np.max(term_exponents)
    term_exponents -= largest_term
    term_exponents = np.delete(term_exponents, np.argmax(term_exponents))
    sum_log = np.log1p(np.sum(np.exp(term_exponents)))
    return -log_z + largest_term + sum_log

def raw_log_weights(betas, energy_series, log_zs):
    raw_log_weights = np.zeros(energy_series.shape)
    for [i, series] in enumerate(energy_series):
        for [j, energy] in enumerate(series):
            raw_log_weights[i, j] = log_denominator_weight(0, energy, betas, len(series), log_zs)
    return raw_log_weights

tol_sqr = 10**-12
def find_sim_log_partition_functions(betas, energy_series):
    new_log_zs = np.zeros(len(betas))
    print(new_log_zs)
    delta_sqr = 1.0
    n_iter = 10000
    log_zs_iter_raw = np.zeros([n_iter, len(betas)])
    n_iter_actual = n_iter
    for _ in range (n_iter):
        old_log_zs = new_log_zs.copy()
        for b in range(len(betas)):
            #print(b)
            new_log_zs[b] = approximate_log_quantity(betas[b], 0, betas, energy_series, old_log_zs)
        deltas = np.exp(new_log_zs - old_log_zs) - 1
        delta_sqr = np.sum(deltas**2)
        print(_)
        print(new_log_zs)
        print(deltas)
        print(delta_sqr)
        log_zs_iter_raw[_] = new_log_zs.copy()
        if delta_sqr < tol_sqr:
            n_iter_actual = _ + 1
            break
    log_zs_iter = log_zs_iter_raw[0:n_iter_actual-1]
    print(1/np.arange(1, n_iter_actual + 1, 1))
    for i in range(len(Ks)):
        plt.plot(1/np.arange(1, n_iter_actual, 1), log_zs_iter[:,i], label=str(Ks[i]))
    return new_log_zs

def approximate_log_quantity(beta, log_z, betas, energy_series, log_zs, observable_series = None):
    if observable_series is None:
        observable_series = np.ones(energy_series.shape)
    log_denominator_series = np.zeros(energy_series.shape)
    for [i, series] in enumerate(energy_series):
        for [j, energy] in enumerate(series):
            log_denominator_series[i, j] = log_denominator_weight(beta, energy, betas, len(series), log_zs)
    term_exponents = (np.log(observable_series) - log_denominator_series).flatten()
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

z_file = h5.File("manual_spacing_log_zs.h5", 'r')
log_zs = z_file.get('log_zs')
log_zs = log_zs[:]



mean_energy_function(Ks, energy_time_series, log_zs, 40)


#log_zs = find_sim_log_partition_functions(Ks, energy_time_series)
#out = h5.File("manual_spacing_log_zs.h5", 'w')
#dset = out.create_dataset(name='log_zs', data=log_zs)


#print(Ks)
#plt.legend()
#plt.show()