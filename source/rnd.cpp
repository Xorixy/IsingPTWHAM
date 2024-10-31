//
// Created by david on 2024-02-28.
//
#include "../include/rnd.h"
#include <valarray>

void rnd::seed(std::optional<uint64_t> number) {
    if(number.has_value()) {
        pcg_extras::seed_seq_from<pcg64> seq(number.value());
        internal::prng.seed(seq);
    }
}

long long unsigned int rnd::next_event(const double p, pcg64 &prng) {
    if (p >= 1.0) return 1;
    if (p <= 0) return std::numeric_limits<long long unsigned int>::max();
    const double P = uniform(0.0, 1.0, prng);
    const double n = std::log(1.0 - P)/std::log(1.0 - p);
    if (n + 1 > std::numeric_limits<long long unsigned int>::max()) {
        return std::numeric_limits<long long unsigned int>::max();
    }
    return static_cast<long long unsigned int>(n + 1);
}

long long unsigned int rnd::next_event(const double p) {
    return next_event(p, internal::prng);
}

