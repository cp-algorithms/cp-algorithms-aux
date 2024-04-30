#ifndef CP_ALGO_RANDOM_RNG_HPP
#define CP_ALGO_RANDOM_RNG_HPP
#include <chrono>
#include <random>
namespace cp_algo::random {
    uint64_t rng() {
        static std::mt19937_64 rng(std::chrono::steady_clock::now().time_since_epoch().count());
        return rng();
    }
}
#endif // CP_ALGO_RANDOM_RNG_HPP
