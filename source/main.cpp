#include "../include/ising.h"
#include <complex>
#include <fmt/ranges.h>
#include "../include/stats.h"

void calculate_omega(std::vector<int> sizes, double K) {
    ising::Ising is(sizes, K, 0);
    for (int in = 0 ; in < 1000000 ; in++) {
        is.run_sim_step(true);
    }
    auto hist = is.get_histogram();
    //fmt::print("{}\n", fmt::join(hist, ", "));
    std::vector<double> fracs;
    int sum = 0;
    for (int i = 0 ; i < hist.size() ; i++) {
        sum += hist[i];
    }
    for (int i = 0 ; i < hist.size() ; i++) {
        fracs.push_back(hist[i]*1.0/sum);
    }
    //fmt::print("{}\n", fmt::join(fracs, ", "));
    std::vector<double> omega;
    double fracsum = 0.0;
    for (int i = 0 ; i < hist.size() ; i++) {
        omega.push_back(fracs[i]*std::exp(2*i*K));
        fracsum += omega[i];
    }
    for (int i = 0 ; i < hist.size() ; i++) {
        omega[i] /= fracsum;
    }
    fmt::print("{}\n", fmt::join(omega, ", "));
}

int main() {
    int size = 10;
    int t_save = size;
    std::vector<int> sizes = {size, size};
    ising::Ising is(sizes, 1.0, 0);
    for (int in = 0 ; in < 100000000 ; in++) {
        if (in % t_save == 0)
            is.run_sim_step(true);
        else
            is.run_sim_step(false);
    }
    auto time_series = is.get_energy_time_series();
    fmt::print("Relaxation time: {}\n", stats::relaxation_time(time_series));
}

