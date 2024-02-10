#ifndef CP_ALGO_RANDOM_RNG_HPP
#define CP_ALGO_RANDOM_RNG_HPP
#include <chrono>
#include <random>
namespace cp_algo::random {
    std::mt19937_64 rng(std::chrono::steady_clock::now().time_since_epoch().count()); 
}
#endif // CP_ALGO_RANDOM_RNG_HPP