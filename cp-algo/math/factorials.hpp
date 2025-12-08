#ifndef CP_ALGO_MATH_FACTORIALS_HPP
#define CP_ALGO_MATH_FACTORIALS_HPP
#pragma GCC push_options
#pragma GCC target("avx2")
#include "../util/checkpoint.hpp"
#include "../util/bump_alloc.hpp"
#include "../util/simd.hpp"
#include "../math/combinatorics.hpp"
#include "../number_theory/modint.hpp"
#include <ranges>

namespace cp_algo::math {
    template<bool use_bump_alloc = false, int maxn = -1>
    auto facts(auto const& args) {
        static_assert(!use_bump_alloc || maxn > 0, "maxn must be set if use_bump_alloc is true");
        constexpr int max_mod = 1'000'000'000;
        constexpr int accum = 4;
        constexpr int simd_size = 8;
        constexpr int block = 1 << 18;
        constexpr int subblock = block / simd_size;
        using base = std::decay_t<decltype(args[0])>;
        static_assert(modint_type<base>, "Base type must be a modint type");
        using T = std::array<int, 2>;
        using alloc = std::conditional_t<use_bump_alloc,
            bump_alloc<T, 30 * maxn>,
            big_alloc<T>>;
        std::basic_string<T, std::char_traits<T>, alloc> odd_args_per_block[max_mod / subblock];
        std::basic_string<T, std::char_traits<T>, alloc> reg_args_per_block[max_mod / subblock];
        constexpr int limit_reg = max_mod / 64;
        int limit_odd = 0;

        big_vector<base> res(size(args), 1);
        const int mod = base::mod();
        const int imod = -math::inv2(mod);
        for(auto [i, xy]: std::views::zip(args, res) | std::views::enumerate) {
            auto [x, y] = xy;
            int t = x.getr();
            if(t >= mod / 2) {
                t = mod - t - 1;
                y = t % 2 ? 1 : mod-1;
            }
            auto pw = 32ull * (t + 1);
            while(t > limit_reg) {
                limit_odd = std::max(limit_odd, (t - 1) / 2);
                odd_args_per_block[(t - 1) / 2 / subblock].push_back({int(i), (t - 1) / 2});
                t /= 2;
                pw += t;
            }
            reg_args_per_block[t / subblock].push_back({int(i), t});
            y *= pow_fixed<base, 2>(int(pw % (mod - 1)));
        }
        checkpoint("init");
        base bi2x32 = pow_fixed<base, 2>(32).inv();
        auto process = [&](int limit, auto &args_per_block, auto step, auto &&proj) {
            base fact = 1;
            for(int b = 0; b <= limit; b += accum * block) {
                u32x8 cur[accum];
                static std::array<u32x8, subblock> prods[accum];
                for(int z = 0; z < accum; z++) {
                    for(int j = 0; j < simd_size; j++) {
                        cur[z][j] = uint32_t(b + z * block + j * subblock);
                        cur[z][j] = proj(cur[z][j]);
                        prods[z][0][j] = cur[z][j] + !cur[z][j];
                        prods[z][0][j] = uint32_t(uint64_t(prods[z][0][j]) * bi2x32.getr() % mod);
                    }
                }
                for(int i = 1; i < block / simd_size; i++) {
                    for(int z = 0; z < accum; z++) {
                        cur[z] += step;
                        prods[z][i] = montgomery_mul(prods[z][i - 1], cur[z], mod, imod);
                    }
                }
                checkpoint("inner loop");
                for(int z = 0; z < accum; z++) {
                    for(int j = 0; j < simd_size; j++) {
                        int bl = b + z * block + j * subblock;
                        for(auto [i, x]: args_per_block[bl / subblock]) {
                            res[i] *= fact * prods[z][x - bl][j];
                        }
                        fact *= base(prods[z].back()[j]);
                    }
                }
                checkpoint("mul ans");
            }
        };
        process(limit_reg, reg_args_per_block, 1, std::identity{});
        process(limit_odd, odd_args_per_block, 2, [](uint32_t x) {return 2 * x + 1;});
        auto invs = bulk_invs<base>(res);
        for(auto [i, x]: res | std::views::enumerate) {
            if (args[i] >= mod / 2) {
                x = invs[i];
            }
        }
        checkpoint("inv ans");
        return res;
    }
}
#pragma GCC pop_options
#endif // CP_ALGO_MATH_FACTORIALS_HPP
