#ifndef CP_ALGO_MATH_XOR_CONVOLUTION_HPP
#define CP_ALGO_MATH_XOR_CONVOLUTION_HPP
#include "../number_theory/modint.hpp"
#include "../util/bit.hpp"
#include "../util/checkpoint.hpp"
#include <cassert>
#include <algorithm>
#include <vector>

namespace cp_algo::math {
    // Recursive FWHT (XOR) transform for size N (power of two)
    template<auto N>
    void xor_transform(auto &&a) {
        if constexpr (N == 1) {
            return;
        } else {
            constexpr auto half = N / 2;
            xor_transform<half>(&a[0]);
            xor_transform<half>(&a[half]);
            for (uint32_t i = 0; i < half; i++) {
                auto x = a[i] + a[i + half];
                auto y = a[i] - a[i + half];
                a[i] = x;
                a[i + half] = y;
            }
        }
    }

    // FWHT wrapper that deduces N at compile time via with_bit_floor
    inline void xor_transform(auto &&a, auto n) {
        with_bit_floor(n, [&]<auto NN>() {
            assert(NN == n);
            xor_transform<NN>(a);
        });
    }

    inline void xor_transform(auto &&a) {
        xor_transform(a, std::size(a));
    }

    // In-place XOR convolution on sequences of equal length (power of two)
    void xor_convolution_inplace(auto &a, auto &b) {
        auto N = static_cast<uint32_t>(std::size(a));
        xor_transform(a);
        xor_transform(b);
        checkpoint("transform");
        for (uint32_t i = 0; i < N; i++) {
            a[i] *= b[i];
        }
        checkpoint("dot");
        xor_transform(a);
        checkpoint("transform");
        using base = std::decay_t<decltype(a[0])>;
        base ni = base(N).inv();
        for (auto &it : a) {
            it *= ni;
        }
        checkpoint("mul_inv");
    }

    // Returns XOR convolution of a and b; pads to next power of two
    auto xor_convolution(auto a, auto b) {
        auto n = std::bit_ceil(std::max(std::size(a), std::size(b)));
        a.resize(n);
        b.resize(n);
        xor_convolution_inplace(a, b);
        return a;
    }
}
#endif // CP_ALGO_MATH_XOR_CONVOLUTION_HPP
