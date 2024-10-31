#include "../include/ising.h"
#include "../include/stats.h"
#include "../include/pt.h"
#include "../include/io.h"
#include "../include/mpi_pt.h"
#include <mpi/mpi.h>


void run_single_core() {
    fmt::print("Starting simulation with {} threads\n", std::getenv("OMP_NUM_THREADS"));
    pt::ParallelIsing pi(settings::constants::sizes, settings::constants::Ks, settings::constants::p_swap, settings::random::seed);
    pi.run_sim(settings::constants::n_steps, settings::constants::n_therm, settings::constants::n_save);
    pi.save_data();
}




int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);
    io::load_settings(settings::io::settings_path);
    pcg_extras::seed_seq_from<pcg64> seq(settings::random::seed);
    rnd::internal::prng.seed(seq);
    mpi_pt::run_mpi_simulation();
    MPI_Finalize();
    return 0;
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    std::vector<int> world_order(world_size, 0);
    if (world_rank == 0) {
        for (int i = 0 ; i < world_size ; i++) {
            world_order[i] = i;
        }
    }

    double K = settings::constants::Ks.at(world_rank);
    fmt::print("Hello from processor {} of {}. I have K value {}\n", world_rank, world_size, K);
    ising::Ising local_ising(settings::constants::sizes, K, settings::random::seed + world_rank);
    local_ising.reserve_time_series(settings::constants::n_steps/static_cast<unsigned long long>(settings::constants::n_save));
    if (settings::constants::n_save >= 1) {
        fmt::print("Process {} starting thermalization\n", world_rank);

        int swap_counter = 0;
        if (world_rank == 0) {
            swap_counter = rnd::next_event(settings::constants::p_swap);
        }
        MPI_Bcast(&swap_counter, 1, MPI_INT, 0, MPI_COMM_WORLD);

        long long unsigned int n_therm = settings::constants::n_therm;
        while(n_therm > swap_counter) {
            n_therm -= swap_counter;
            local_ising.run_step(swap_counter, false);
        }


        fmt::print("Process {} thermalization done.\nStarting simulation\n", world_rank);

        if (settings::constants::n_save == 1) {
            local_ising.run_step(settings::constants::n_steps, true);

        } else {

            while(settings::constants::n_steps > settings::constants::n_save) {
                local_ising.run_step(settings::constants::n_save - 1, false);
                local_ising.run_step(1, true);
                settings::constants::n_steps -= settings::constants::n_save;
            }
        }
        fmt::print("Process {} simulation done.\n", world_rank);
    }
    mpi_pt::try_swap(local_ising, world_order);
    local_ising.run_step(1, true);
    for (int rank = 0; rank < world_size; ++rank) {
        if(world_rank == rank) {
            fmt::print("Process {} saving data.\n", world_rank);
            local_ising.save_data(fmt::format("{}", rank));
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }


    MPI_Finalize();

    return 0;
}

