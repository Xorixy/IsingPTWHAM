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
    for (auto &is : m_isings) {
        is.run_step(n_steps, save_data);
    }
}

void pt::ParallelIsing::swap_random() {
    int n_swap = rnd::uniform(0, static_cast<int>(m_isings.size()) - 2);
    ising::Ising::swap_states(m_isings[n_swap], m_isings[n_swap+1]);
}


void pt::ParallelIsing::run_sim(long long int n_steps, long long int n_therm, long int n_save) {
    long int swap_counter = rnd::next_event(m_swap);
    while (n_therm - swap_counter >= 0) {
        n_therm -= swap_counter;
        run_step(swap_counter, false);
        swap_random();
        swap_counter = rnd::next_event(m_swap);
    }
    if (n_therm > 0) {
        run_step(n_therm, false);
        swap_counter -= n_therm;
    }
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
        swap_counter = rnd::next_event(m_swap);
    }
}

void pt::ParallelIsing::save_data() {
    for (auto &ising : m_isings) {
        ising.save_data();
    }
}




