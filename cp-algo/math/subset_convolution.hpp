#ifndef CP_ALGO_MATH_SUBSET_CONVOLUTION_HPP
#define CP_ALGO_MATH_SUBSET_CONVOLUTION_HPP
#include "../util/simd.hpp"
#include "../util/big_alloc.hpp"
#include "../util/bit.hpp"
#include "../util/checkpoint.hpp"
#include <array>
#include <ranges>
#include <algorithm>
#include <bit>
#include <cstring>
CP_ALGO_SIMD_PRAGMA_PUSH
namespace cp_algo::math {
    const size_t max_logn = 20;
    
    enum transform_dir { forw, inv };
    
    template<auto N, transform_dir direction>
    inline void xor_transform(auto &&a) {
        [[gnu::assume(N <= 1 << 30)]];
        if constexpr (N <= 32) {
            for(size_t i = 1; i < N; i *= 2) {
                for(size_t j = 0; j < N; j += 2 * i) {
                    for(size_t k = j; k < j + i; k++) {
                        for(size_t z = 0; z < max_logn; z++) {
                            auto x = a[k][z] + a[k + i][z];
                            auto y = a[k][z] - a[k + i][z];
                            a[k][z] = x;
                            a[k + i][z] = y;
                        }
                    }
                }
            }
        } else {
            constexpr auto half = N / 2;
            xor_transform<half, direction>(&a[0]);
            xor_transform<half, direction>(&a[half]);
            for (size_t i = 0; i < half; i++) {
                #pragma GCC unroll max_logn
                for(size_t z = 0; z < max_logn; z++) {
                    auto x = a[i][z] + a[i + half][z];
                    auto y = a[i][z] - a[i + half][z];
                    a[i][z] = x;
                    a[i + half][z] = y;
                }
            }
        }
    }
    
    template<transform_dir direction>
    inline void xor_transform(auto &&a, auto n) {
        with_bit_floor(n, [&]<auto NN>() {
            assert(NN == n);
            xor_transform<NN, direction>(a);
        });
    }
    
    template<transform_dir direction = forw>
    inline void xor_transform(auto &&a) {
        xor_transform<direction>(a, std::size(a));
    }

    // Generic rank vectors processor with variadic inputs
    // Assumes output[0] = 0, caller is responsible for handling rank 0
    // Returns the output array
    auto on_rank_vectors(auto &&cb, auto const& ...inputs) {
        static_assert(sizeof...(inputs) >= 1, "on_rank_vectors requires at least one input");
        
        // Create tuple of input references once
        auto input_tuple = std::forward_as_tuple(inputs...);
        auto const& first_input = std::get<0>(input_tuple);
        using base = std::decay_t<decltype(first_input[0])>;
        big_vector<base> out(std::size(first_input));
        
        auto N = std::size(first_input);
        constexpr size_t K = 4;
        N = std::max(N, 2 * K);
        const size_t n = std::bit_width(N) - 1;
        const size_t T = std::min<size_t>(n - 3, 2);
        const size_t bottoms = 1 << (n - T - 1);
        const auto M = std::size(first_input);
        
        // Create array buffers for each input
        auto create_buffers = [bottoms]<typename... Args>(const Args&...) {
            return std::make_tuple(
                big_vector<std::array<typename std::decay_t<Args>::value_type, max_logn>>(bottoms)...
            );
        };
        auto buffers = std::apply(create_buffers, input_tuple);
        
        checkpoint("alloc buffers");
        big_vector<uint32_t> counts(2 * bottoms);
        for(size_t i = 1; i < 2 * bottoms; i++) {
            counts[i] = (uint32_t)std::popcount(i);
        }
        checkpoint("prepare");
        
        for(size_t top = 0; top < N / 2; top += bottoms) {
            // Clear all buffers
            std::apply([bottoms](auto&... bufs) {
                (..., memset(bufs.data(), 0, sizeof(bufs[0]) * bottoms));
            }, buffers);
            checkpoint("memset");
            
            // Initialize buffers from inputs
            std::apply([&](auto const&... inps) {
                std::apply([&](auto&... bufs) {
                    auto init_one = [&](auto const& inp, auto& buf) {
                        for(size_t i = 0; i < M; i += 2 * bottoms) {
                            bool parity = __builtin_parity(uint32_t((i >> 1) & top));
                            size_t limit = std::min(M, i + 2 * bottoms) - i;
                            uint32_t count = (uint32_t)std::popcount(i) - 1;
                            for(size_t bottom = (i == 0); bottom < limit; bottom++) {
                                if (parity) {
                                    buf[bottom >> 1][count + counts[bottom]] -= inp[i + bottom];
                                } else {
                                    buf[bottom >> 1][count + counts[bottom]] += inp[i + bottom];
                                }
                            }
                        }
                    };
                    (init_one(inps, bufs), ...);
                }, buffers);
            }, input_tuple);
            
            checkpoint("init");
            std::apply([](auto&... bufs) {
                (..., xor_transform(bufs));
            }, buffers);
            checkpoint("transform");
            
            assert(bottoms % K == 0);
            for(size_t i = 0; i < bottoms; i += K) {
                std::apply([&](auto&... bufs) {
                    auto extract_one = [&](auto& buf) {
                        std::array<u64x4, max_logn> aa;
                        for(size_t j = 0; j < max_logn; j++) {
                            for(size_t z = 0; z < K; z++) {
                                aa[j][z] = buf[i + z][j].getr();
                            }
                        }
                        return aa;
                    };
                    
                    auto aa_tuple = std::make_tuple(extract_one(bufs)...);
                    std::apply(cb, aa_tuple);
                    
                    // Write results back: only first array needs to be written
                    auto& first_buf = std::get<0>(std::forward_as_tuple(bufs...));
                    const auto& first_aa = std::get<0>(aa_tuple);
                    for(size_t j = 0; j < max_logn; j++) {
                        for(size_t z = 0; z < K; z++) {
                            first_buf[i + z][j].setr((uint32_t)first_aa[j][z]);
                        }
                    }
                }, buffers);
            }
            
            checkpoint("dot");
            auto& first_buf = std::get<0>(buffers);
            xor_transform<inv>(first_buf);
            checkpoint("transform");
            
            // Gather results from first buffer

            for(size_t i = 0; i < M; i += 2 * bottoms) {
                bool parity = __builtin_parity(uint32_t((i >> 1) & top));
                size_t limit = std::min(M, i + 2 * bottoms) - i;
                uint32_t count = (uint32_t)std::popcount(i) - 1;
                for(size_t bottom = (i == 0); bottom < limit; bottom++) {
                    if (parity) {
                        out[i + bottom] -= first_buf[bottom >> 1][count + counts[bottom]];
                    } else {
                        out[i + bottom] += first_buf[bottom >> 1][count + counts[bottom]];
                    }
                }
            }
            checkpoint("gather");
        }
        const base ni = base(N / 2).inv();
        for(auto& x : out) {x *= ni;}
        return out;
    }

    template<typename base>
    big_vector<base> subset_convolution(std::span<base> inpa, std::span<base> inpb) {
        big_vector<base> outpa;
        with_bit_floor(std::size(inpa), [&]<auto N>() {
            constexpr size_t lgn = std::bit_width(N) - 1;
            [[gnu::assume(lgn <= max_logn)]];
            outpa = on_rank_vectors([](auto &a, auto const& b) {
                std::decay_t<decltype(a)> res = {};
                const auto mod = base::mod();
                const auto imod = math::inv2(-mod);
                const auto r4 = u64x4() + uint64_t(-1) % mod + 1;
                auto add = [&](size_t i) {
                    if constexpr (lgn) for(size_t j = 0; i + j + 1 < lgn; j++) {
                        res[i + j + 1] += (u64x4)_mm256_mul_epu32(__m256i(a[i]), __m256i(b[j]));
                    }
                };
                if constexpr (lgn) for(size_t i = 0; i < lgn; i++) { add(i); }
                if constexpr (lgn) if constexpr (lgn) for(size_t k = 0; k < lgn; k++) {
                    res[k] = montgomery_reduce(res[k], mod, imod);
                    res[k] = montgomery_mul(res[k], r4, mod, imod);
                    a[k] = res[k] >= mod ? res[k] - mod : res[k];
                }
            }, inpa, inpb);
            
            outpa[0] = inpa[0] * inpb[0];
            for(size_t i = 1; i < std::size(inpa); i++) {
                outpa[i] += inpa[i] * inpb[0] + inpa[0] * inpb[i];
            }
            checkpoint("fix 0");
        });
        return outpa;
    }

    template<typename base>
    big_vector<base> subset_exp(std::span<base> inpa) {
        if (size(inpa) == 1) {
            return big_vector<base>{1};
        }
        size_t N = std::size(inpa);
        auto out0 = subset_exp(std::span(inpa).first(N / 2));
        auto out1 = subset_convolution<base>(out0, std::span(inpa).last(N / 2));
        out0.insert(end(out0), begin(out1), end(out1));
        cp_algo::checkpoint("extend out");
        return out0;
    }

    template<typename base>
    big_vector<big_vector<base>> subset_compose(big_vector<std::span<base>> fd, std::span<base> inpa) {
        if (size(inpa) == 1) {
            big_vector<big_vector<base>> res(size(fd), {base(0)});
            big_vector<base> pw(size(fd[0]), 1);
            for (size_t i = 1; i < size(fd[0]); i++) {
                pw[i] = pw[i - 1] * inpa[0];
            }
            for (size_t i = 0; i < size(fd); i++) {
                for (size_t j = 0; j < size(fd[i]); j++) {
                    res[i][0] += pw[j] * fd[i][j];
                }
            }
            cp_algo::checkpoint("base case");
            return res;
        }
        size_t N = std::size(inpa);
        big_vector<base> fdk(size(fd[0]));
        for (size_t i = 0; i + 1 < size(fdk); i++) {
            fdk[i] = fd.back()[i + 1] * base(i + 1);
        }
        fd.push_back(fdk);
        cp_algo::checkpoint("fdk");
        auto deeper = subset_compose(fd, std::span(inpa).first(N / 2));
        for(size_t i = 0; i + 1 < size(fd); i++) {
            auto next = subset_convolution<base>(deeper[i + 1], std::span(inpa).last(N / 2));
            deeper[i].insert(end(deeper[i]), begin(next), end(next));
        }
        deeper.pop_back();
        cp_algo::checkpoint("combine");
        return deeper;
    }

    template<typename base>
    big_vector<base> subset_compose(std::span<base> f, std::span<base> inpa) {
        return subset_compose(big_vector{f}, inpa)[0];
    }
}
#pragma GCC pop_options
#endif // CP_ALGO_MATH_SUBSET_CONVOLUTION_HPP
