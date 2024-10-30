#pragma once
#include "ising.h"
#include <vector>
#include <omp.h>
#include <chrono>

namespace pt {

    class ParallelIsing {
        private:
            std::vector<ising::Ising> m_isings;
            double m_swap;
            int m_successful_swaps {0};
        public:
            ParallelIsing(const std::vector<int>& sizes, const std::vector<double> & Ks, double p_swap, uint64_t seed);
            void run_step(bool save_data);
            void run_step(long long int n_steps, bool save_data);
            void run_sim(long long int n_steps, long long int n_therm, long int n_save);
            void run_sim(long long int n_steps, long long int n_therm, long int n_save, double p_swap);
            void swap_random();
            void save_data();
    };
}
