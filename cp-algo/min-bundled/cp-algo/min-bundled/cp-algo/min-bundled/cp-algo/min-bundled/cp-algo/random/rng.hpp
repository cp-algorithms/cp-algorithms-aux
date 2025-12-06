#line 1 "cp-algo/min-bundled/cp-algo/min-bundled/cp-algo/min-bundled/cp-algo/random/rng.hpp"
#line 1 "cp-algo/min-bundled/cp-algo/min-bundled/cp-algo/random/rng.hpp"
#line 1 "cp-algo/min-bundled/cp-algo/random/rng.hpp"
#line 1 "cp-algo/random/rng.hpp"
#include <chrono>
#include <random>
namespace cp_algo::random{std::mt19937_64 gen(
std::chrono::steady_clock::now().time_since_epoch().count()
);
uint64_t rng(){return gen();}}