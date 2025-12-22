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
    const size_t logn = 20;
    
    enum transform_dir { forw, inv };
    
    template<auto N, transform_dir direction>
    inline void or_transform(auto &&a) {
        [[gnu::assume(N <= 1ull << 30)]];
        if constexpr (N <= 32) {
            for(size_t i = 1; i < N; i *= 2) {
                for(size_t j = 0; j < N; j += 2 * i) {
                    for(size_t k = j; k < j + i; k++) {
                        for(size_t z = 0; z < logn; z++) {
                            if constexpr (direction == forw) {
                                a[k + i][z] += a[k][z];
                            } else {
                                a[k + i][z] -= a[k][z];
                            }
                        }
                    }
                }
            }
        } else {
            constexpr auto half = N / 2;
            or_transform<half, direction>(&a[0]);
            or_transform<half, direction>(&a[half]);
            for (size_t i = 0; i < half; i++) {
                #pragma GCC unroll logn
                for(size_t z = 0; z < logn; z++) {
                    if constexpr (direction == forw) {
                        a[i + half][z] += a[i][z];
                    } else {
                        a[i + half][z] -= a[i][z];
                    }
                }
            }
        }
    }
    
    template<transform_dir direction>
    inline void or_transform(auto &&a, auto n) {
        with_bit_floor(n, [&]<auto NN>() {
            assert(NN == n);
            or_transform<NN, direction>(a);
        });
    }
    
    template<transform_dir direction = forw>
    inline void or_transform(auto &&a) {
        or_transform<direction>(a, std::size(a));
    }

    // Generic rank vectors processor with variadic inputs
    // Assumes output[0] = 0, caller is responsible for handling rank 0
    // Returns the output array
    auto on_rank_vectors(auto &&cb, auto const& ...inputs) {
        static_assert(sizeof...(inputs) >= 1, "on_rank_vectors requires at least one input");
        
        // Create tuple of input references once
        auto input_tuple = std::forward_as_tuple(inputs...);
        auto const& first_input = std::get<0>(input_tuple);
        
        auto out = first_input;
        using base = std::decay_t<decltype(first_input[0])>;
        std::ranges::fill(out, base(0));
        
        auto N = std::size(first_input);
        constexpr size_t K = 4;
        N = std::max(N, K);
        const size_t n = std::bit_width(N) - 1;
        const size_t T = std::min<size_t>(n - 2, 4);
        const size_t bottoms = 1 << (n - T);
        const auto M = std::size(first_input);
        
        // Create array buffers for each input
        auto create_buffers = [bottoms]<typename... Args>(const Args&...) {
            return std::make_tuple(
                big_vector<std::array<typename std::decay_t<Args>::value_type, logn>>(bottoms)...
            );
        };
        auto buffers = std::apply(create_buffers, input_tuple);
        
        big_vector<uint32_t> counts(N);
        for(size_t i = 1; i < N; i++) {
            counts[i] = (uint32_t)std::popcount(i);
        }
        checkpoint("prepare");
        
        for(size_t top = 0; top < N; top += bottoms) {
            // Clear all buffers
            std::apply([bottoms](auto&... bufs) {
                (..., memset(bufs.data(), 0, sizeof(bufs[0]) * bottoms));
            }, buffers);
            
            // Initialize buffers from inputs
            std::apply([&](auto const&... inps) {
                std::apply([&](auto&... bufs) {
                    auto init_one = [&](auto const& inp, auto& buf) {
                        for(size_t mask = top; ; mask = (mask - bottoms) & top) {
                            size_t limit = std::min(M, mask + bottoms) - mask;
                            uint32_t count = counts[mask / bottoms] - 1;
                            for(size_t bottom = (mask == 0); bottom < limit; bottom++) {
                                size_t i = bottom | mask;
                                buf[bottom][count + counts[bottom]] += inp[i];
                            }
                            if (!mask) break;
                        }
                    };
                    (init_one(inps, bufs), ...);
                }, buffers);
            }, input_tuple);
            
            checkpoint("init");
            std::apply([](auto&... bufs) {
                (..., or_transform(bufs));
            }, buffers);
            checkpoint("transform");
            
            assert(bottoms % K == 0);
            for(size_t i = 0; i < bottoms; i += K) {
                std::apply([&](auto&... bufs) {
                    auto extract_one = [&](auto& buf) {
                        std::array<u64x4, logn> aa;
                        for(size_t j = 0; j < logn; j++) {
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
                    for(size_t j = 0; j < logn; j++) {
                        for(size_t z = 0; z < K; z++) {
                            first_buf[i + z][j].setr((uint32_t)first_aa[j][z]);
                        }
                    }
                }, buffers);
            }
            
            checkpoint("dot");
            auto& first_buf = std::get<0>(buffers);
            or_transform<inv>(first_buf);
            checkpoint("transform");
            
            // Gather results from first buffer
            for(size_t mask = top; mask < N; mask = (mask + bottoms) | top) {
                bool parity = __builtin_parity(uint32_t(mask ^ top));
                size_t limit = std::min(M, mask + bottoms) - mask;
                uint32_t count = counts[mask / bottoms] - 1;
                for(size_t bottom = (mask == 0); bottom < limit; bottom++) {
                    size_t i = bottom | mask;
                    if (parity) {
                        out[i] -= first_buf[bottom][count + counts[bottom]];
                    } else {
                        out[i] += first_buf[bottom][count + counts[bottom]];
                    }
                }
            }
            checkpoint("gather");
        }
        return out;
    }

    template<typename base>
    auto subset_convolution(auto const& inpa, auto const& inpb) {
        auto outpa = on_rank_vectors([](auto &a, auto const& b) {
            std::decay_t<decltype(a)> res = {};
            const auto mod = base::mod();
            const auto imod = math::inv2(-mod);
            const auto r4 = u64x4() + uint64_t(-1) % mod + 1;
            for(size_t i = 0; i < logn; i++) {
                for(size_t j = 0; i + j + 1 < logn; j++) {
                    res[i + j + 1] += (u64x4)_mm256_mul_epu32(__m256i(a[i]), __m256i(b[j]));
                }
                if (i == logn / 2) {
                    for(size_t k = logn / 2; k < logn; k++) {
                        res[k] = res[k] >= base::modmod8() ? res[k] - base::modmod8() : res[k];
                    }
                }
            }
            for(size_t k = 0; k < logn; k++) {
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
        return outpa;
    }

    template<typename base>
    auto subset_exp(auto const& inpa) {
        auto outpa = on_rank_vectors([](auto &p) {
            std::decay_t<decltype(p)> q = {};
            const auto mod = base::mod();
            const auto imod = math::inv2(-mod);
            const auto r2 = uint32_t(-1) % mod + 1;
            const auto r4 = uint64_t(-1) % mod + 1;
            static const auto invs = [&]() {
                std::array<uint64_t, logn> invs;
                for(size_t i = 0; i < logn; i++) {
                    invs[i] = (base(i + 1).inv() * base(r4)).getr();
                }
                return invs;
            }();
            for(size_t k = 0; k < logn; k++) {
                q[k] = p[k] = montgomery_mul(p[k], u64x4() + (k + 1) * r2 % mod, mod, imod);
                for(size_t i = 0; i < k; i++) {
                    q[k] += (u64x4)_mm256_mul_epu32((__m256i)p[i], (__m256i)q[k - i - 1]);
                    if (i == logn / 2) {
                        q[k] = q[k] >= base::modmod8() ? q[k] - base::modmod8() : q[k];
                    }
                }
                q[k] = montgomery_reduce(q[k], mod, imod);
                q[k] = montgomery_mul(q[k], u64x4() + invs[k], mod, imod);
                q[k] = q[k] >= mod ? q[k] - mod : q[k];
            }
            p = q;
        }, inpa);
        outpa[0] = base(1);
        return outpa;
    }
}
#pragma GCC pop_options
#endif // CP_ALGO_MATH_SUBSET_CONVOLUTION_HPP
