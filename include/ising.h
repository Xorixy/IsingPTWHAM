#pragma once
#include <vector>
#include <fmt/core.h>
#include <cassert>
#include <pcg_random.hpp>
#include "rnd.h"
#include "io.h"

namespace ising {
    class State {
        private:
            std::vector<unsigned char> m_spins;
            std::vector<unsigned char*> m_neighbours;
            std::vector<int> m_neighbour_ids;
            std::vector<int> m_sizes;
            int m_size { 0 };
            int m_dim;
            int m_energy { 0 };

            void create_neighbours();
        public:
            explicit State(const std::vector<int>& sizes);
            int flip_random_spin(double K, pcg64& prng) noexcept;
            [[nodiscard]] int get_energy() const noexcept;
            [[nodiscard]] int get_size() const noexcept;
            [[nodiscard]] int get_dim() const noexcept;
            void print() const noexcept;
            void recalculate_energy() noexcept;
            void print_neighbours() const noexcept;
            void print_neighbour_spins() const noexcept;
            void set_spin(int point, unsigned char spin);
    };
    class Ising {
        private:
            State m_state;
            double m_K;
            pcg64 m_prng;
            std::vector<long long unsigned int> m_histogram;
            std::vector<int> m_energy_time_series;
        public:
            Ising(const std::vector<int>& sizes, double K, const pcg64 & prng);
            Ising(const std::vector<int>& sizes, double K, uint64_t seed);
            [[nodiscard]] std::vector<int> const& get_energy_time_series() const noexcept;
            [[nodiscard]] std::vector<long long unsigned int> const& get_histogram() const noexcept;
            [[nodiscard]] int get_size() const noexcept;
            [[nodiscard]] int get_energy() const noexcept;
            [[nodiscard]] double get_K() const noexcept;

            void reserve_time_series(int size) noexcept;
            void run_step(bool save_data) noexcept;
            void run_step(long long int n_steps, bool save_data) noexcept;

            static void swap_states(Ising& ising_one, Ising& ising_two) noexcept;
            void print_state() const noexcept;
            void print_neighbours() const noexcept;
            void print_neighbour_spins() const noexcept;
            void set_spin(int point, unsigned char spin);
            void save_data(const std::string& prefix);
            void save_data();
    };
}