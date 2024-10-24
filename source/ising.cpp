#include "../include/ising.h"

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

ising::State::State(const std::vector<int>& sizes)
    : m_sizes (sizes)
    , m_dim {static_cast<int>(sizes.size())}
{
    m_size = 1;
    for (const auto s : m_sizes) {
        assert(s > 0);
        m_size *= s;
    }
    m_spins = std::vector<unsigned char>(m_size, 0);
    m_neighbours = std::vector<unsigned char*>(2*m_dim*m_size);
    #ifndef NDEBUG
    for (const auto i : m_spins) {
        assert(!i);
    }
    #endif
    create_neighbours();
}

void ising::State::create_neighbours() {
    std::vector iter_vec(m_dim, 0);
    std::vector neighbour_vec(m_dim, 0);
    std::vector multiplier_vec(m_dim, 1);
    for (int iN = 0 ; iN < m_dim - 1 ; iN++) {
        for (int id = iN+1 ; id < m_dim ; id++) {
            multiplier_vec.at(id) *= m_sizes.at(iN);
        }
    }
    for (int _ = 0 ; _ < m_size ; _++) {
        int iter_int = 0;
        for (int ix = 0 ; ix < m_dim ; ix++) {
            iter_int += multiplier_vec.at(ix) * iter_vec.at(ix);
        }
        for (int id = 0 ; id < m_dim ; id++) {
            for (const int delta : {-1, 1}) {
                neighbour_vec = iter_vec;
                neighbour_vec.at(id) = (neighbour_vec.at(id) + delta + m_sizes.at(id)) % m_sizes.at(id);
                int neighbour_int = 0;

                for (int ix = 0 ; ix < m_dim ; ix++) {
                    neighbour_int += multiplier_vec.at(ix) * neighbour_vec.at(ix);
                }
                int shift = m_dim * ( delta == -1);
                m_neighbours.at(2*m_dim*iter_int + shift + id) = &m_spins.at(neighbour_int);
            }
        }
        for (int id = 0 ; id < m_dim ; id++) {
            iter_vec.at(id)++;
            if (iter_vec.at(id) == m_sizes.at(id))
                iter_vec.at(id) = 0;
            else
                break;
        }
    }
}

int ising::State::flip_random_spin(double K, pcg64& prng) noexcept {
    int point_flip = rnd::uniform(0, m_dim-1, prng);
    int e_diff = 0;
    unsigned char spin = m_spins[point_flip];
    for (int in = 0 ; in < 2*m_dim ; in++) {
        e_diff += 2*(spin == *m_neighbours[2*m_dim*point_flip + in]) - 1;
    }
    if (sgn(K*e_diff) >= 0 || rnd::uniform(0.0, 1.0) < exp(K*e_diff)) {
        m_spins[point_flip] = 1 - m_spins[point_flip];
        return e_diff;
    }
    return 0;
}

int ising::State::get_energy() const noexcept {
    return m_energy;
}

int ising::State::get_size() const noexcept {
    return m_size;
}


ising::Ising::Ising(const std::vector<int> &sizes, double K, pcg64 &prng)
    : m_state(sizes)
    , m_K { K }
    , m_prng { prng }
    , m_histogram(m_state.get_size())
    , m_energy_time_series(0)
{}

/*
ising::Ising::Ising(const std::vector<int> &sizes, double K, uint64_t seed)
    : m_state(sizes)
    , m_K{K}
    , m_prng()
    , m_histogram(1 + sizes.size()*m_state.get_size())
    , m_energy_time_series(0)
{
    pcg_extras::seed_seq_from<pcg64> seq(seed);
    m_prng.seed(seq);
}
*/


std::vector<int> ising::Ising::get_energy_time_series() const noexcept {
    return m_energy_time_series;
}

std::vector<long long> ising::Ising::get_histogram() const noexcept {
    return m_histogram;
}

int ising::Ising::get_size() const noexcept {
    return m_state.get_size();
}


void ising::Ising::swap_states(Ising &ising_one, Ising &ising_two) noexcept {
    assert(ising_one.get_size() == ising_two.get_size());
    std::swap(ising_one.m_state, ising_two.m_state);
}

void ising::Ising::reserve_time_series(const int size) noexcept {
    m_energy_time_series.reserve(size);
}

void ising::Ising::run_sim_step(bool save_data) {
    m_state.flip_random_spin(m_K, m_prng);
    if (save_data) {
        const int state_energy = m_state.get_energy();
        m_histogram.at((state_energy + m_histogram.size())/2)++;
        m_energy_time_series.push_back(state_energy);
    }
}
