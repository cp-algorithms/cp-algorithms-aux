// Licensed under GPL-3.0-only with the online-judge exception.
// See cp-algo/LICENSE and cp-algo/ONLINE-JUDGE-EXCEPTION.
#ifndef CP_ALGO_RANDOM_RNG_HPP
#define CP_ALGO_RANDOM_RNG_HPP
#include <chrono>
#include <random>
namespace cp_algo::random {
    std::mt19937_64 gen(
        std::chrono::steady_clock::now().time_since_epoch().count()
    );
    uint64_t rng() {
        return gen();
    }
}
#endif // CP_ALGO_RANDOM_RNG_HPP
