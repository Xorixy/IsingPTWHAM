#pragma once
#include "ising.h"
#include <mpi/mpi.h>

namespace mpi_pt {
    inline int world_rank = -1;
    inline int world_size = -1;

    enum class WalkDirection {
        up,
        down,
        none,
    };

    void open_file();

    void try_swap(ising::Ising & local_ising, std::vector<int>& world_order);

    void run_mpi_simulation();
    void init_Ks();
    void update_swap_counter(unsigned long long int &swap_counter);
    void save_direction_numbers(WalkDirection &dir, std::vector<int>& world_order, std::vector<int> &n_up, std::vector<int> &n_down);
    void run_thermalization(ising::Ising & local_ising, std::vector<int>& world_order, long long unsigned int &swap_counter);
    void run_equilibration(ising::Ising & local_ising, std::vector<int>& world_order, long long unsigned int &swap_counter);
    void run_production(ising::Ising & local_ising, std::vector<int>& world_order, long long unsigned int &swap_counter);
    void save_data(ising::Ising & local_ising);
}