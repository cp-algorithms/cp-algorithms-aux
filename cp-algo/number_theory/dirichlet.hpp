#ifndef CP_ALGO_NUMBER_THEORY_DIRICHLET_HPP
#define CP_ALGO_NUMBER_THEORY_DIRICHLET_HPP
#include <numeric>
#include <cstdint>
#include <vector>
#include <cmath>
namespace cp_algo::math {    
    auto floor_stats(int64_t n) {
        auto rt_n = int(sqrtl(n));
        return std::pair{rt_n, 2 * rt_n - (n / rt_n == rt_n)};
    }

    struct interval {
        int lo, hi;
    };

    void exec_on_blocks(int64_t n, auto &&callback) {
        auto [rt_n, num_floors] = floor_stats(n);

        auto to_ord = [&](int64_t k) {
            return k <= rt_n ? int(k) : num_floors - int(n / k) + 1;
        };

        callback({1, 1}, {1, 1}, {1, num_floors});

        for (int k = 2; k <= num_floors; ++k) {
            if(k > rt_n) {
                int z = num_floors - k + 1;
                for (int x = 2; ; x++) {
                    int y_lo_ord = std::max(x, z) + 1;
                    int y_hi_ord = to_ord(n / (x * z));
                    if (y_hi_ord < y_lo_ord) break;
                    callback({x, x}, {y_lo_ord, y_hi_ord}, {k, k});
                    callback({y_lo_ord, y_hi_ord}, {x, x}, {k, k});
                }
            }

            callback({1, 1}, {k, k}, {k, num_floors});
            callback({k, k}, {1, 1}, {k, num_floors});

            if(k <= rt_n) {
                int x = k;
                for (int y = 2; y < k; ++y) {
                    int z_lo_ord = to_ord(1LL * x * y);
                    int z_hi_ord = to_ord(n / x);
                    if (z_hi_ord < z_lo_ord) break;
                    callback({x, x}, {y, y}, {z_lo_ord, z_hi_ord});
                    callback({y, y}, {x, x}, {z_lo_ord, z_hi_ord});
                }
                int z_lo_ord = to_ord(1LL * x * x);
                callback({x, x}, {x, x}, {z_lo_ord, num_floors});
            }
        }
    }

    auto Dirichlet_mul(auto &&F, auto &&G, int64_t n) {
        auto m = size(F);
        std::decay_t<decltype(F)> H(m+1);
        exec_on_blocks(n, [&](interval x, interval y, interval z) {
            auto sum_x = F[x.hi] - F[x.lo - 1];
            auto sum_y = G[y.hi] - G[y.lo - 1];
            auto t = sum_x * sum_y;
            H[z.lo] += t;
            H[z.hi + 1] -= t;
        });
        std::partial_sum(begin(H), end(H), begin(H));
        H.pop_back();
        return H;
    }

    auto Dirichlet_div(auto &&H, auto &&G, int64_t n) {
        H.push_back(0);
        std::adjacent_difference(begin(H), end(H), begin(H));
        auto m = size(G);
        std::decay_t<decltype(G)> F(m);
        std::vector<bool> assigned(m);
        auto Gi = G[1].inv();
        exec_on_blocks(n, [&](interval x, interval y, interval z) {
            if (!assigned[x.hi]) {
                F[x.hi] = F[x.lo - 1] + H[z.lo] * Gi;
                assigned[x.hi] = true;
            }
            auto sum_x = F[x.hi] - F[x.lo - 1];
            auto sum_y = G[y.hi] - G[y.lo - 1];
            auto t = sum_y * sum_x;
            H[z.lo] -= t;
            H[z.hi + 1] += t;
        });
        return F;
    }

}
#endif // CP_ALGO_NUMBER_THEORY_DIRICHLET_HPP
