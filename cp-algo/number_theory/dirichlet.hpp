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
        auto operator <=>(const interval&) const = default;
    };

    // callback(k, prefix) such that:
    // (F * G)[k] = prefix + (F[k] - F[k-1]) * G[1] + (G[k] - G[k-1]) * F[1]
    // Uses H as a buffer for (F * G)[k], then overrides with callback results
    void exec_on_blocks(int64_t n, auto &H, auto const& F, auto const& G, auto &&callback) {
        auto [rt_n, num_floors] = floor_stats(n);

        auto to_ord = [&](int64_t k) {
            return k <= rt_n ? int(k) : num_floors - int(n / k) + 1;
        };

        auto call = [&](interval x, interval y, interval z) {
            auto sum_x = F[x.hi] - F[x.lo - 1];
            auto sum_y = G[y.hi] - G[y.lo - 1];
            auto t = sum_x * sum_y;
            H[z.lo] += t;
            if (z.hi < num_floors) {
                H[z.hi + 1] -= t;
            }
        };

        auto prefix = F[1] * G[1];
        for (int k = 2; k <= num_floors; ++k) {
            if(k > rt_n) {
                int z = num_floors - k + 1;
                for (int x = 2; ; x++) {
                    int y_lo_ord = std::max(x, z) + 1;
                    int y_hi_ord = to_ord(n / (x * z));
                    if (y_hi_ord < y_lo_ord) break;
                    call({x, x}, {y_lo_ord, y_hi_ord}, {k, k});
                    call({y_lo_ord, y_hi_ord}, {x, x}, {k, k});
                }
            }

            H[k] = callback(k, prefix += H[k]);
            prefix += (F[k] - F[k-1]) * G[1] + (G[k] - G[k-1]) * F[1];

            if(k <= rt_n) {
                int x = k;
                for (int y = 2; y < k; ++y) {
                    int z_lo_ord = to_ord(1LL * x * y);
                    int z_hi_ord = to_ord(n / x);
                    if (z_hi_ord < z_lo_ord) break;
                    call({x, x}, {y, y}, {z_lo_ord, z_hi_ord});
                    call({y, y}, {x, x}, {z_lo_ord, z_hi_ord});
                }
                int z_lo_ord = to_ord(1LL * x * x);
                call({x, x}, {x, x}, {z_lo_ord, num_floors});
            }
        }
    }

    auto Dirichlet_mul(auto &&F, auto &&G, int64_t n) {
        auto m = size(F);
        std::decay_t<decltype(F)> H(m);
        H[1] = F[1] * G[1];
        exec_on_blocks(n, H, F, G, [&](auto k, auto prefix) {
            return prefix + (F[k] - F[k - 1]) * G[1] + (G[k] - G[k - 1]) * F[1];
        });
        return H;
    }

    auto Dirichlet_div(auto &&H, auto &&G, int64_t n) {
        auto m = size(G);
        std::decay_t<decltype(G)> F(m);
        auto Gi = G[1].inv();
        F[1] = Gi * H[1];
        exec_on_blocks(n, F, F, G, [&](auto k, auto prefix) {
            return F[k-1] + (Gi * (H[k] - prefix - (G[k] - G[k-1]) * F[1]));
        });
        return F;
    }

}
#endif // CP_ALGO_NUMBER_THEORY_DIRICHLET_HPP
