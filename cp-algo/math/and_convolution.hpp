#ifndef CP_ALGO_MATH_AND_CONVOLUTION_HPP
#define CP_ALGO_MATH_AND_CONVOLUTION_HPP
#include "../number_theory/modint.hpp"
#include "../util/bit.hpp"
#include "../util/checkpoint.hpp"

namespace cp_algo::math {
    enum transform_dir { forward, inverse };
    // Recursive Zeta / Möbius (AND) transform for size N (power of two)
    template<auto N, transform_dir direction>
    void and_transform(auto &&a) {
        if constexpr (N == 1) {
            return;
        } else {
            constexpr auto half = N / 2;
            and_transform<half, direction>(&a[0]);
            and_transform<half, direction>(&a[half]);
            for (uint32_t i = 0; i < half; i++) {
                if constexpr (direction == forward) {
                    a[i] += a[i + half];
                } else {
                    a[i] -= a[i + half];
                }
            }
        }
    }
    // Wrapper that deduces N at compile time via with_bit_floor
    template<transform_dir direction>
    inline void and_transform(auto &&a, auto n) {
        with_bit_floor(n, [&]<auto NN>() {
            assert(NN == n);
            and_transform<NN, direction>(a);
        });
    }
    template<transform_dir direction = forward>
    inline void and_transform(auto &&a) {
        and_transform<direction>(a, std::size(a));
    }
    // In-place AND convolution on sequences of equal length (power of two)
    void and_convolution_inplace(auto &a, auto &b) {
        auto N = static_cast<uint32_t>(std::size(a));
        and_transform(a);
        and_transform(b);
        checkpoint("transform");
        for (uint32_t i = 0; i < N; i++) {
            a[i] *= b[i];
        }
        checkpoint("dot");
        and_transform<inverse>(a);
        checkpoint("itransform");
    }

    // Returns AND convolution of a and b; pads to next power of two
    auto and_convolution(auto a, auto b) {
        auto n = std::bit_ceil(std::max(std::size(a), std::size(b)));
        a.resize(n);
        b.resize(n);
        and_convolution_inplace(a, b);
        return a;
    }
}
#endif // CP_ALGO_MATH_AND_CONVOLUTION_HPP
