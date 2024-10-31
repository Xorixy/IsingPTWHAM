#include "../include/mpi_pt.h"

void mpi_pt::try_swap(ising::Ising & local_ising, std::vector<int>& world_order) {
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    if (world_order.size() < 2) return;

    std::vector<int> energies(world_order.size(), 0);

    if (world_rank != 0) {
        const int energy = local_ising.get_energy();
        MPI_Send(&energy, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
    } else {
        energies[0] = local_ising.get_energy();
        for (int i = 1; i < world_order.size(); i++) {
            MPI_Recv(&energies[i], 1, MPI_INT, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    if (world_rank == 0) {
        int n_swap = rnd::uniform(0, static_cast<int>(world_order.size() - 2));
        int delta_e = energies[world_order[n_swap]] - energies[world_order[n_swap + 1]];
        double delta_K = settings::constants::Ks[n_swap] - settings::constants::Ks[n_swap + 1];
        if (rnd::uniform(0.0, 1.0) < exp(delta_K * delta_K)) {
            std::swap(world_order[n_swap], world_order[n_swap + 1]);
        }
        for (int i = 0; i < world_order.size(); i++) {
            if (world_order[i] == 0) {
                local_ising.set_K(settings::constants::Ks[i]);
            } else {
                MPI_Send(&settings::constants::Ks[i], 1, MPI_DOUBLE, world_order[i], 0, MPI_COMM_WORLD);
            }
        }
    } else {
        double new_K;
        MPI_Recv(&new_K, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        local_ising.set_K(new_K);
    }
}

void mpi_pt::update_swap_counter(unsigned long long int &swap_counter) {
    if (world_rank == 0) {
        swap_counter = rnd::next_event(settings::constants::p_swap);
    }
    MPI_Bcast(&swap_counter, 1, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
}

void mpi_pt::run_thermalization(ising::Ising &local_ising, std::vector<int> &world_order, long long unsigned int &swap_counter) {
    long long unsigned int n_therm = settings::constants::n_therm;
    while (n_therm >= swap_counter) {
        n_therm -= swap_counter;
        local_ising.run_step(swap_counter, false);
        try_swap(local_ising, world_order);
        update_swap_counter(swap_counter);
    }
    swap_counter -= n_therm;
    local_ising.run_step(n_therm, false);
}

void mpi_pt::run_simulation(ising::Ising &local_ising, std::vector<int> &world_order, long long unsigned int &swap_counter) {
    long long unsigned int n_steps = settings::constants::n_steps;
    long long unsigned int n_save = settings::constants::n_save;
    if (n_save == 1) {
        while (n_steps >= swap_counter) {
            n_steps -= swap_counter;
            local_ising.run_step(swap_counter, true);
            try_swap(local_ising, world_order);
            update_swap_counter(swap_counter);
        }
        swap_counter -= n_steps;
        local_ising.run_step(n_steps, true);
    } else {
        while (n_steps >= n_save) {
            while (n_save >= swap_counter) {
                n_save -= swap_counter;
                local_ising.run_step(swap_counter, false);
                try_swap(local_ising, world_order);
                update_swap_counter(swap_counter);
            }
            local_ising.run_step(n_save - 1, false);
            local_ising.run_step(1, true);
            swap_counter -= n_save;
            n_save = settings::constants::n_save;
        }
    }
}



void mpi_pt::run_mpi_simulation() {
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    std::vector<int> world_order(world_size, 0);
    if(world_rank == 0) {
        for(int i = 0; i < world_size; i++) { world_order[i] = i; }
    }
    const std::vector<double> Ks = settings::constants::Ks;
    double K = Ks[world_rank];
    ising::Ising local_ising(settings::constants::sizes, K, settings::random::seed + world_rank);

    long long unsigned int swap_counter = 0;
    run_thermalization(local_ising, world_order, swap_counter);
    run_simulation(local_ising, world_order, swap_counter);
    save_data(local_ising);
}

void mpi_pt::save_data(ising::Ising &local_ising) {
    for (int rank = 0; rank < world_size; ++rank) {
        if(world_rank == rank) {
            fmt::print("Process {} saving data.\n", world_rank);
            local_ising.save_data(fmt::format("{}", rank));
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}