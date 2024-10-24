#pragma once
#include <optional>
#include <pcg_random.hpp>
#include <random>

namespace rnd {
    namespace internal {
        inline pcg64 prng; // The pseudorandom number generator
    }

    void seed(std::optional<uint64_t> number = std::nullopt);

    template<typename T>
    requires std::is_arithmetic_v<T>
    [[nodiscard]] T uniform(T min, T max) {
        if constexpr(std::is_integral_v<T>) {
            std::uniform_int_distribution<T> dist(min, max);
            return dist(internal::prng);
        } else if constexpr(std::is_floating_point_v<T>) {
            std::uniform_real_distribution<T> dist(min, max);
            return dist(internal::prng);
        }
    }

    template<typename T>
    requires std::is_arithmetic_v<T>
    [[nodiscard]] T uniform(T min, T max, pcg64& prng) {
        if constexpr(std::is_integral_v<T>) {
            std::uniform_int_distribution<T> dist(min, max);
            return dist(internal::prng);
        } else if constexpr(std::is_floating_point_v<T>) {
            std::uniform_real_distribution<T> dist(min, max);
            return dist(internal::prng);
        }
    }
}