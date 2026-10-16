#ifndef CP_ALGO_MATH_XOR_CONVOLUTION_HPP
#define CP_ALGO_MATH_XOR_CONVOLUTION_HPP
#include "../util/bit.hpp"
#include "../util/checkpoint.hpp"

namespace cp_algo::math {
    // Recursive FWHT (XOR) transform for size N (power of two)
    template<auto N>
    void xor_transform(auto &&a) {
        if constexpr (N == 1) {
            return;
        } else if constexpr (N == 2) {
            auto x = a[0] + a[1], y = a[0] - a[1];
            a[0] = x;
            a[1] = y;
        } else {
            constexpr auto q = N / 4;
            for (size_t j = 0; j < 4; j++) {
                xor_transform<q>(&a[j * q]);
            }
            // Combine two stages without storing the intermediate butterflies.
            for (size_t i = 0; i < q; i++) {
                auto x0 = a[i] + a[i + q], x1 = a[i] - a[i + q];
                auto x2 = a[i + 2 * q] + a[i + 3 * q], x3 = a[i + 2 * q] - a[i + 3 * q];
                a[i] = x0 + x2;
                a[i + q] = x1 + x3;
                a[i + 2 * q] = x0 - x2;
                a[i + 3 * q] = x1 - x3;
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
