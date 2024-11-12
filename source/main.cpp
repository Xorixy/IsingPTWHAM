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
}

