#include "../include/pt.h"

pt::ParallelIsing::ParallelIsing(const std::vector<int> &sizes, const std::vector<double> &Ks, double p_swap, uint64_t seed)
    : m_swap { p_swap }
{
    for (int ik = 0 ; ik < Ks.size(); ++ik) {
        m_isings.emplace_back(sizes, Ks[ik], seed + ik*Ks.size());
    }
}

void pt::ParallelIsing::run_step(bool save_data) {
    run_step(1, save_data);
}


void pt::ParallelIsing::run_step(long long int n_steps, bool save_data) {
    if (omp_get_thread_num() > 1) {
        #pragma omp parallel for
        for (auto &is : m_isings) {
            is.run_step(n_steps, save_data);
        }
    } else {
        for (auto &is : m_isings) {
            is.run_step(n_steps, save_data);
        }
    }
}

void pt::ParallelIsing::swap_random() {
    int n_swap = rnd::uniform(0, static_cast<int>(m_isings.size()) - 2);
    int delta_e = m_isings[n_swap + 1].get_energy() - m_isings[n_swap].get_energy();
    double delta_K = m_isings[n_swap + 1].get_K() - m_isings[n_swap].get_K();
    if (rnd::uniform(0.0, 1.0) < exp(delta_K*delta_e)) {
        //fmt::print("Swapping");
        ising::Ising::swap_states(m_isings[n_swap], m_isings[n_swap+1]);
        m_successful_swaps++;
    }
}

void pt::ParallelIsing::run_sim(long long int n_steps, long long int n_therm, long int n_save) {
    if (m_swap > 0 && m_swap <= 1.0) {
        run_sim(n_steps, n_therm, n_save, m_swap);
    } else {
        std::chrono::time_point<std::chrono::steady_clock> start_time;
        start_time = std::chrono::steady_clock::now();
        fmt::print("Starting thermalization...");

        run_step(n_therm, false);

        std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start_time;
        fmt::print("Done in {} s\n", elapsed.count());
        start_time = std::chrono::steady_clock::now();
        fmt::print("Starting simulation...");
        auto save_counter = n_save;
        while (n_steps - save_counter >= 0) {
            n_steps -= save_counter;
            run_step(save_counter - 1, false);
            run_step(1, true);
        }
        elapsed = std::chrono::steady_clock::now() - start_time;
        fmt::print("Done in {} s\n", elapsed.count());
    }
}


void pt::ParallelIsing::run_sim(long long int n_steps, long long int n_therm, long int n_save, double p_swap) {
    std::chrono::time_point<std::chrono::steady_clock> start_time;
    start_time = std::chrono::steady_clock::now();
    fmt::print("Starting thermalization...");
    long int swap_counter = rnd::next_event(p_swap);
    while (n_therm - swap_counter >= 0) {
        n_therm -= swap_counter;
        run_step(swap_counter, false);
        swap_random();
        swap_counter = rnd::next_event(p_swap);
    }
    if (n_therm > 0) {
        run_step(n_therm, false);
        swap_counter -= n_therm;
    }
    std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start_time;
    fmt::print("Done in {} s\n", elapsed.count());
    start_time = std::chrono::steady_clock::now();
    fmt::print("Starting simulation...");
    auto save_counter = n_save;
    while (n_steps - swap_counter >= 0) {
        n_steps -= swap_counter;
        while(swap_counter - save_counter >= 0) {
            run_step(save_counter - 1, false);
            run_step(1, true);
            swap_counter -= save_counter;
            save_counter = n_save;
        }
        save_counter -= swap_counter;
        run_step(swap_counter, false);
        swap_random();
        swap_counter = rnd::next_event(p_swap);
    }
    elapsed = std::chrono::steady_clock::now() - start_time;
    fmt::print("Done in {} s\n", elapsed.count());
    fmt::print("Number of successful swaps {}\n", m_successful_swaps);

}

void pt::ParallelIsing::save_data() {
    for (auto &ising : m_isings) {
        ising.save_data();
    }
}




