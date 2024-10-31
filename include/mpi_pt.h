#pragma once
#include "ising.h"
#include <mpi/mpi.h>

namespace mpi_pt {
    inline int world_rank = -1;
    inline int world_size = -1;

    void open_file();

    void try_swap(ising::Ising & local_ising, std::vector<int>& world_order);

    void run_mpi_simulation();

    void update_swap_counter(unsigned long long int &swap_counter);

    void run_thermalization(ising::Ising & local_ising, std::vector<int>& world_order, long long unsigned int &swap_counter);
    void run_simulation(ising::Ising & local_ising, std::vector<int>& world_order, long long unsigned int &swap_counter);
    void save_data(ising::Ising & local_ising);
}