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

int rnd::next_event(const double p, pcg64 &prng) {
    assert(p > 0.0 && p < 1.0);
    const double P = uniform(0.0, 1.0, prng);
    const double n = std::log(1.0 - P)/std::log(1.0 - p);
    return 1 + static_cast<int>(n);
}

int rnd::next_event(const double p) {
    return next_event(p, internal::prng);
}

