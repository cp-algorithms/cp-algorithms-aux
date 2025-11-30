#ifndef CP_ALGO_NUMBER_THEORY_DIRICHLET_HPP
#define CP_ALGO_NUMBER_THEORY_DIRICHLET_HPP
#include <algorithm>
#include <cstdint>
#include <ranges>
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

    // callback(k) when:
    //     (F * G)[k] = H[k] + (F[k] - F[k-1]) * G[1] + (G[k] - G[k-1]) * F[1]
    // Return the value to be saved in H[k]
    enum exec_mode { standard, reverse };
    template<exec_mode mode = standard>
    void exec_on_blocks(int64_t n, auto &H, auto const& F, auto const& G, auto &&callback) {
        auto [rt_n, num_floors] = floor_stats(n);

        auto to_ord = [&](int64_t k) {
            return k <= rt_n ? int(k) : num_floors - int(n / k) + 1;
        };

        auto call = [&](interval x, interval y, interval z) {
            auto Fx = F[x.hi] - F[x.lo - 1];
            auto Fy = F[y.hi] - F[y.lo - 1];
            decltype(Fx) Gx, Gy;
            if constexpr (mode == standard) {
                Gy = G[y.hi] - G[y.lo - 1];
                Gx = G[x.hi] - G[x.lo - 1];
            } else {
                Gy = G[y.lo - 1] - G[y.hi];
                Gx = G[x.lo - 1] - G[x.hi];
            }
            auto t = Fx * Gy;
            if(x != y) [[likely]] {
                t += Fy * Gx;
            }
            H[z.lo] += t;
            if (z.hi < num_floors) [[likely]] {
                H[z.hi + 1] -= t;
            }
        };
        for (int k = 2; k <= num_floors; ++k) {
            if(k > rt_n) {
                int z = num_floors - k + 1;
                for (int x = 2; ; x++) {
                    int y_lo_ord = std::max(x, z) + 1;
                    int y_hi_ord = to_ord(n / (x * z));
                    if (y_hi_ord < y_lo_ord) break;
                    call({x, x}, {y_lo_ord, y_hi_ord}, {k, k});
                }
            }

            H[k] = callback(k);

            if(k <= rt_n) {
                int x = k;
                for (int y = 2; y < k; ++y) {
                    int z_lo_ord = to_ord(1LL * x * y);
                    int z_hi_ord = to_ord(n / x);
                    if (z_hi_ord < z_lo_ord) break;
                    call({x, x}, {y, y}, {z_lo_ord, z_hi_ord});
                }
                int z_lo_ord = to_ord(1LL * x * x);
                call({x, x}, {x, x}, {z_lo_ord, num_floors});
            }
        }
    }

    auto Dirichlet_mul(auto const& F, auto const& G, int64_t n) {
        auto m = size(F);
        std::decay_t<decltype(F)> H(m);
        H[1] = F[1] * G[1];
        exec_on_blocks(n, H, F, G, [&](auto k) {
            return H[k] + (F[k] - F[k-1]) * G[1] + (G[k] - G[k-1]) * F[1];
        });
        partial_sum(begin(H), end(H), begin(H));
        return H;
    }

    void Dirichlet_div_inplace(auto &H, auto const& G, int64_t n) {
        auto Gi = G[1].inv();
        H[0] -= H[0];
        adjacent_difference(begin(H), end(H), begin(H));
        H[1] *= Gi;
        exec_on_blocks<reverse>(n, H, H, G, [&](auto k) {
            return (Gi * (H[k] - (G[k] - G[k-1]) * H[1])) + H[k-1];
        });
    }

    auto Dirichlet_div(auto const& H, auto const& G, int64_t n) {
        auto m = size(G);
        auto F = H | std::views::take(m) | std::ranges::to<std::decay_t<decltype(G)>>();
        Dirichlet_div_inplace(F, G, n);
        return F;
    }
}
#endif // CP_ALGO_NUMBER_THEORY_DIRICHLET_HPP
