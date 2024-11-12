#include "../include/mpi_pt.h"

void mpi_pt::open_file() {
    if (world_rank == 0) {
        io::outfile = io::try_to_open_file(settings::io::outfile, true);
        settings::io::outfile = io::outfile.getFilePath();
    }
    int outfile_size = settings::io::outfile.size();
    MPI_Bcast(&outfile_size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (world_rank != 0)
        settings::io::outfile.resize(outfile_size);
    MPI_Bcast(const_cast<char*>(settings::io::outfile.data()), outfile_size, MPI_CHAR, 0, MPI_COMM_WORLD);
    if (world_rank != 0) {
        io::outfile = io::try_to_open_file(settings::io::outfile, false);
    }
}


void mpi_pt::try_swap(ising::Ising & local_ising, std::vector<int>& world_order) {
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
        //fmt::print("Trying to swap\n");
        int n_swap = rnd::uniform(0, static_cast<int>(world_order.size() - 2));
        int delta_e = energies[world_order[n_swap]] - energies[world_order[n_swap + 1]];
        double delta_K = settings::constants::Ks[n_swap] - settings::constants::Ks[n_swap + 1];
        if (rnd::uniform(0.0, 1.0) < exp(delta_K * delta_e)) {
            //fmt::print("Swap successful\n");
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

    while (n_therm >= swap_counter && n_therm != 0) {
        n_therm -= swap_counter;
        local_ising.run_step(swap_counter, false);
        try_swap(local_ising, world_order);
        update_swap_counter(swap_counter);
    }
    swap_counter -= n_therm;
    local_ising.run_step(n_therm, false);
}

void mpi_pt::run_production(ising::Ising &local_ising, std::vector<int> &world_order, long long unsigned int &swap_counter) {
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
        while (n_steps >= n_save && n_steps != 0) {
            n_steps -= n_save;
            while (n_save >= swap_counter) {
                n_save -= swap_counter;
                local_ising.run_step(swap_counter, false);
                try_swap(local_ising, world_order);
                update_swap_counter(swap_counter);
            }
            local_ising.run_step((n_save > 0)*(n_save - 1), false);
            local_ising.run_step(1, true);
            swap_counter -= n_save;
            n_save = settings::constants::n_save;
        }
    }
}

void mpi_pt::save_direction_numbers(WalkDirection &dir, std::vector<int>& world_order, std::vector<int> &n_up, std::vector<int> &n_down) {
    int dir_int = 0;
    if (dir == WalkDirection::up) {
        dir_int = 1;
    } else if (dir == WalkDirection::down) {
        dir_int = -1;
    }
    std::vector<int> dirs;
    if (world_rank == 0) {
        dirs = std::vector(world_size, 0);
    }
    MPI_Gather(&dir_int, 1, MPI_INT, &dirs[0], 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (world_rank == 0) {
        for (int i = 0; i < dirs.size(); i++) {
            if (dirs[world_order[i]] == 1) {
                n_up[i]++;
            } else if (dirs[world_order[i]] == -1) {
                n_down[i]++;
            }
        }
    }
}


void mpi_pt::run_equilibration(ising::Ising &local_ising, std::vector<int> &world_order, long long unsigned int &swap_counter) {
    long long unsigned int n_therm = settings::constants::n_therm;
    WalkDirection dir = WalkDirection::none;
    if (local_ising.get_K() == settings::constants::K_min) {
        dir = WalkDirection::up;
    } else if (local_ising.get_K() == settings::constants::K_max) {
        dir = WalkDirection::down;
    }

    std::vector<int> n_up(world_size, 0);
    std::vector<int> n_down(world_size, 0);
    int n_equilibration_runs = 1;
    for (int _ = 0 ; _ < n_equilibration_runs; _++) {
        n_therm = settings::constants::n_therm;
        while (n_therm >= swap_counter && n_therm != 0) {
            n_therm -= swap_counter;
            local_ising.run_step(swap_counter, false);
            save_direction_numbers(dir, world_order, n_up, n_down);
            try_swap(local_ising, world_order);
            if (local_ising.get_K() == settings::constants::K_min) {
                dir = WalkDirection::up;
            } else if (local_ising.get_K() == settings::constants::K_max) {
                dir = WalkDirection::down;
            }
            update_swap_counter(swap_counter);
        }
        swap_counter -= n_therm;
        local_ising.run_step(n_therm, false);
    }
    if (world_rank == 0) {
        std::vector<double> f(world_size, 0.0);
        for (int i = 0 ; i < world_size ; i++) {
            f[i] = n_up[i]*1.0/(n_up[i] + n_down[i]);
        }
        fmt::print("f : {}\n", fmt::join(f, ","));
    }
}

void mpi_pt::run_mpi_simulation() {
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    open_file();

    std::vector<int> world_order(world_size, 0);
    if(world_rank == 0) {
        for(int i = 0; i < world_size; i++) { world_order[i] = i; }
    }
    init_Ks();

    double K = settings::constants::Ks[world_rank];
    ising::Ising local_ising(settings::constants::sizes, K, settings::random::seed + world_rank);

    long long unsigned int swap_counter = 0;
    update_swap_counter(swap_counter);
    fmt::print("Process {} starting thermalization.\n", world_rank);
    run_equilibration(local_ising, world_order, swap_counter);
    fmt::print("Process {} thermalization done.\nStarting simulation.\n", world_rank);
    run_production(local_ising, world_order, swap_counter);
    fmt::print("Process {} simulation done.\n", world_rank);
    MPI_Barrier(MPI_COMM_WORLD);
    save_data(local_ising);
}


void mpi_pt::save_data(ising::Ising &local_ising) {
    for (int rank = 0; rank < world_size; ++rank) {
        if(world_rank == rank) {
            fmt::print("Process {} saving data.\n", world_rank);
            int n_tries = 0;
            while(true) {
                try {
                    local_ising.save_data(fmt::format("{}", rank));
                    break;
                } catch (std::exception &e) {
                    fmt::print("Failed to save data with error: {}\n", e.what());
                    fmt::print("Waiting 3 seconds and trying again...\n");
                    std::this_thread::sleep_for(std::chrono::seconds(3));
                    n_tries++;
                    if (n_tries >= 10) throw e;
                }
            }
            fmt::print("Process {} data saved.\n", world_rank);
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

void mpi_pt::init_Ks() {
    if (world_size == settings::constants::Ks.size()) {
        return;
    }
    std::vector<double> Ks(world_size, 0);
    if (world_size == 1) {
        Ks[0] = settings::constants::K_min;
        settings::constants::Ks = Ks;
    }
    else {
        for (int i = 0 ; i < world_size ; i++) {
            Ks[i] = settings::constants::K_min + (settings::constants::K_max - settings::constants::K_min)/(world_size - 1)*i;
        }
        Ks[world_size - 1] = settings::constants::K_max;
        settings::constants::Ks = Ks;
        return;
        double K_min = settings::constants::K_min;
        double K_max = settings::constants::K_max;
        double R = pow(K_max/K_min, 1.0/(world_size - 1));
        for (int k = 0 ; k < world_size; k++) {
            Ks[world_size - k - 1] = K_max + K_min - K_min*pow(R, k);
        }
        Ks[0] = settings::constants::K_min;
        Ks[world_size - 1] = settings::constants::K_max;
        settings::constants::Ks = Ks;
    }
}
