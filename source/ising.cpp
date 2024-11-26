#include "../include/ising.h"

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

ising::State::State(const std::vector<int>& sizes)
    : m_sizes (sizes)
    , m_dim {static_cast<int>(sizes.size())}
    , m_magnetization( 0 )
{
    m_size = 1;
    for (const auto s : m_sizes) {
        assert(s > 0);
        m_size *= s;
    }
    m_magnetization = -m_size;
    m_spins = std::vector<unsigned char>(m_size, 0);
    m_neighbours = std::vector<unsigned char*>(2*m_dim*m_size);
    m_neighbour_ids = std::vector<int>(2*m_dim*m_size);
    #ifndef NDEBUG
    for (const auto i : m_spins) {
        assert(!i);
    }
    #endif
    create_neighbours();
    recalculate_energy();
}

void ising::State::set_spin(int point, unsigned char spin) {
    m_spins[point] = spin;
}


void ising::State::recalculate_energy() noexcept {
    m_energy = 0;
    for (int ip = 0 ; ip < m_spins.size(); ip++) {
        unsigned char spin = m_spins[ip];
        for (int in = 0 ; in < 2*m_dim; in++) {
            m_energy -= 2*(spin == *m_neighbours[2*m_dim*ip + in]) - 1;
        }
    }
    assert(m_energy % 2 == 0);
    m_energy /= 2;
}

void ising::State::print_neighbours() const noexcept {
    for (int ip = 0 ; ip < m_spins.size(); ip++) {
        for (int in = 0 ; in < 2*m_dim; in++) {
            fmt::print("{}, ", m_neighbour_ids[2*m_dim*ip + in]);
        }
        fmt::print("\n");
    }
}

void ising::State::print_neighbour_spins() const noexcept {
    for (int ip = 0 ; ip < m_spins.size(); ip++) {
        for (int in = 0 ; in < 2*m_dim; in++) {
            unsigned char spin = *m_neighbours[2*m_dim*ip + in];
            if (spin == 1)
                fmt::print("  1, ");
            else if (spin == 0)
                fmt::print(" -1, ");
            else
                fmt::print("  ?, ");
        }
        fmt::print("\n");
    }
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
                m_neighbour_ids.at(2*m_dim*iter_int + shift + id) = neighbour_int;
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

int ising::State::flip_random_spin(const double K, pcg64& prng) noexcept {
    int point_flip = rnd::uniform(0, m_size-1, prng);
    int e_diff = 0;
    unsigned char spin = m_spins[point_flip];
    for (int in = 0 ; in < 2*m_dim ; in++) {
        e_diff += 4*(spin == *m_neighbours[2*m_dim*point_flip + in]) - 2;
    }
    if (sgn(-K*e_diff) >= 0 || rnd::uniform(0.0, 1.0) < exp(-K*e_diff)) {
        m_magnetization += 2 - 4*m_spins[point_flip];
        m_spins[point_flip] = 1 - m_spins[point_flip];
        m_energy += e_diff;
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

int ising::State::get_dim() const noexcept {
    return m_dim;
}

int ising::State::get_magnetization() const noexcept {
    return m_magnetization;
}


void ising::State::print() const noexcept {
    fmt::print("energy : {}, magnetization : {}\n", m_energy, m_magnetization);
    if (m_dim == 1) {
        for (int is = 0 ; is < m_size ; is++) {
            if (m_spins.at(is) == 1)
                fmt::print("  1  ");
            else if (m_spins.at(is) == 0)
                fmt::print(" -1 ");
            else
                fmt::print("  ?  ");
        }
        fmt::print("\n");
    }
    if (m_dim == 2 && m_sizes[0] == 2 && m_sizes[1] == 2) {
        for (int i = 0 ; i < 2 ; i++) {
            for (int j = 0 ; j < 2 ; j++) {
                if (m_spins.at(2*i + j) == 1)
                    fmt::print("  1  ");
                else if (m_spins.at(2*i + j) == 0)
                    fmt::print(" -1 ");
                else
                    fmt::print("  ?  ");
            }
            fmt::print("\n");
        }
    }
}

ising::Ising::Ising(const std::vector<int> &sizes, const double K, const pcg64 &prng)
    : m_state(sizes)
    , m_K { K }
    , m_prng { prng }
    , m_energy_time_series(0)
    , m_magnetization_time_series(0)
{}


ising::Ising::Ising(const std::vector<int> &sizes, double K, uint64_t seed)
    : m_state(sizes)
    , m_K{K}
    , m_energy_time_series(0)
    , m_magnetization_time_series(0)
{
    pcg_extras::seed_seq_from<pcg64> seq(seed);
    m_prng.seed(seq);
}



std::vector<int> const& ising::Ising::get_energy_time_series() const noexcept {
    return m_energy_time_series;
}

int ising::Ising::get_size() const noexcept {
    return m_state.get_size();
}

int ising::Ising::get_energy() const noexcept {
    return m_state.get_energy();
}

void ising::Ising::swap_states(Ising &ising_one, Ising &ising_two) noexcept {
    assert(ising_one.get_size() == ising_two.get_size());
    std::swap(ising_one.m_state, ising_two.m_state);
}

void ising::Ising::reserve_time_series(const int size) noexcept {
    m_energy_time_series.reserve(size);
    m_magnetization_time_series.reserve(size);
    m_K_time_series.reserve(size);
}

void ising::Ising::run_step(const bool save_data) noexcept {
    m_state.flip_random_spin(m_K, m_prng);
    if (save_data) {
        const int state_energy = m_state.get_energy();
        m_energy_time_series.push_back(state_energy);

        const int state_magnetization = m_state.get_magnetization();
        m_magnetization_time_series.push_back(state_magnetization);

        m_K_time_series.push_back(m_K);
    }
}

void ising::Ising::run_step(const long long unsigned int n_steps, const bool save_data) noexcept {
    for (long long int _ = 0 ; _ < n_steps ; _++) {
        run_step(save_data);
    }
}

void ising::Ising::print_state() const noexcept {
    m_state.print();
}
void ising::Ising::print_neighbours() const noexcept {
    m_state.print_neighbours();
}

void ising::Ising::print_neighbour_spins() const noexcept {
    m_state.print_neighbour_spins();
}

void ising::Ising::set_spin(int point, unsigned char spin) {
    m_state.set_spin(point, spin);
}

void ising::Ising::save_data(const std::string &prefix) {
    auto energy_series_name = prefix + "/energy_series";
    auto magnetization_series_name = prefix + "/magnetization_series";
    auto K_series_name = prefix + "/K_series";

    io::save_vector(m_energy_time_series, energy_series_name);
    io::save_vector(m_magnetization_time_series, magnetization_series_name);
    io::save_vector(m_K_time_series, K_series_name);
}

void ising::Ising::save_data() {
    save_data(std::to_string(m_K));
}

double ising::Ising::get_K() const noexcept {
    return m_K;
}

void ising::Ising::set_K(const double K) {
    m_K = K;
}

