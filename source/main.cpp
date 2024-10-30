#include "../include/ising.h"
#include <complex>
#include <fmt/ranges.h>
#include "../include/stats.h"
#include "../include/pt.h"
#include "../include/io.h"
#include <nlohmann/json.hpp>
#include "../include/rnd.h"

void calculate_omega(std::vector<int> sizes, double K) {
    ising::Ising is(sizes, K, 0);
    for (int in = 0 ; in < 1000000 ; in++) {
        is.run_step(true);
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
    io::outfile = io::try_to_open_file(settings::io::outfile, false);
    pt::ParallelIsing pi(settings::constants::sizes, settings::constants::Ks, 0.01, settings::random::seed);
    pi.run_sim(100000, 100000, 1000);
    pi.save_data();
    io::outfile.writeDataset<int>(1, "test");
}

