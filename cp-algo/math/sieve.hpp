#ifndef CP_ALGO_MATH_SIEVE_HPP
#define CP_ALGO_MATH_SIEVE_HPP
#include "../structures/bit_array.hpp"
#include "../structures/bit_array_util.hpp"
#include "../util/bit.hpp"
#include "../util/checkpoint.hpp"
#include <cstdint>
#include <cstddef>
#include <vector>
#include <span>
#include <algorithm>
#include <cassert>

CP_ALGO_BIT_PRAGMA_PUSH
namespace cp_algo::math {
    using cp_algo::structures::dynamic_bit_array;
    using cp_algo::structures::bit_array;

    constexpr uint32_t period = 210;
    constexpr uint32_t coprime = 48;
    constexpr auto coprime210 = [](auto x) {
        return x % 2 && x % 3 && x % 5 && x % 7;
    };

    // Residues coprime to 210
    constexpr auto res210 = []() {
        std::array<uint8_t, coprime> res;
        int idx = 0;
        for(uint8_t i = 1; i < period; i += 2) {
            if (coprime210(i)) {
                res[idx++] = i;
            }
        }
        return res;
    }();

    // Maps residue mod 210 to pre-upper_bound index in res210
    constexpr auto state210 = []() {
        std::array<uint8_t, period> state;
        uint8_t idx = 0;
        for(uint8_t i = 0; i < period; i++) {
            state[i] = idx;
            idx += coprime210(i);
        }
        return state;
    }();

    // Add to reach next coprime residue
    constexpr auto add210 = []() {
        std::array<uint8_t, period> add;
        for(uint8_t i = 0; i < period; i++) {
            add[i] = 1;
            while (!coprime210(i + add[i])) {
                add[i]++;
            }
        }
        return add;
    }();

    constexpr auto gap210 = []() {
        std::array<uint8_t, coprime> gap;
        for(uint8_t i = 0; i < coprime; i++) {
            gap[i] = add210[res210[i]];
        }
        return gap;
    }();

    // Convert value to ordinal (index in compressed bit array)
    constexpr uint32_t to_ord(uint32_t x) {
        return (x / period) * coprime + state210[x % period];
    }

    // Convert ordinal to value
    constexpr uint32_t to_val(uint32_t x) {
        return (x / coprime) * period + res210[x % coprime];
    }

    constexpr size_t sqrt_threshold = 1 << 16;
    constexpr auto sqrt_prime_bits = []() {
        const int size = sqrt_threshold / 4;
        bit_array<size> prime;
        prime.set_all();
        prime.reset(to_ord(1));
        for(uint32_t i = 11; to_ord(i * i) < size; i += add210[i % period]) {
            if (prime[to_ord(i)]) {
                for(uint32_t k = i; to_ord(i * k) < size; k += add210[k % period]) {
                    prime.reset(to_ord(i * k));
                }
            }
        }
        return prime;
    }();

    constexpr size_t num_primes = []() {
        size_t cnt = 0;
        for(uint32_t i = 11; i < sqrt_threshold; i += add210[i % period]) {
            cnt += sqrt_prime_bits[to_ord(i)];
        }
        return cnt;
    }();
    constexpr auto sqrt_primes = []() {
        std::array<uint32_t, num_primes> primes;
        size_t j = 0;
        for(uint32_t i = 11; i < sqrt_threshold; i += add210[i % period]) {
            if (sqrt_prime_bits[to_ord(i)]) {
                primes[j++] = i;
            }
        }
        return primes;
    }();

    struct wheel_t {
        dynamic_bit_array mask;
        uint32_t product;
    };

    constexpr auto make_wheel(big_vector<uint32_t> primes, uint32_t product) {
        assert(product % period == 0 && product / period * coprime % dynamic_bit_array::width == 0);
        wheel_t wheel;
        wheel.product = product;
        wheel.mask.resize(product / period * coprime);
        wheel.mask.set_all();
        for(auto p: primes) {
            for (uint32_t k = 1; p * k < product; k += add210[k % period]) {
                wheel.mask.reset(to_ord(p * k));
            }
        }
        return wheel;
    }
    
    constexpr void sieve_dense(auto &prime, uint32_t l, uint32_t r, wheel_t const& wheel) {
        if (l >= r) return;
        const uint32_t width = (uint32_t)dynamic_bit_array::width; // 64
        uint32_t wl = l / width;
        uint32_t wr = (r + width - 1) / width;
        uint32_t N  = (uint32_t)wheel.mask.words;
        auto loop = [&](uint32_t i, uint32_t block) {
            auto p_ptr = std::assume_aligned<32>(&prime.word(i));
            auto m_ptr = std::assume_aligned<32>(&wheel.mask.word(0));
            #pragma GCC unroll coprime
            for (uint32_t j = 0; j < block; j++) {
                p_ptr[j] &= m_ptr[j];
            }
        };
        while (wl + N <= wr) {
            loop(wl, N);
            wl += N;
        }
        loop(wl, wr - wl);
    }

    template <class BitArray>
    constexpr void sieve210(BitArray& prime, uint32_t l, uint32_t r, size_t i, int state) {
        static const auto ord_step = []() {
            big_vector<std::array<uint32_t, 2 * coprime>> ord_steps(num_primes);
            for (uint32_t i = 0; i < size(sqrt_primes); i++) {
                auto p = sqrt_primes[i];
                auto &ords = ord_steps[i];
                auto last = to_ord(p);
                for(uint32_t j = 0; j < coprime; j++) {
                    auto next = to_ord(p * (res210[j] + gap210[j]));
                    ords[j] = ords[j + coprime] = next - last;
                    last = next;
                }
            }
            return ord_steps;
        }();
        auto advance = [&]() {
            prime.reset(std::exchange(l, l + ord_step[i][state++]));
        };
        uint32_t p = sqrt_primes[i];
        while (l + p * coprime <= r) {
            #pragma GCC unroll coprime
            for (size_t j = 0; j < coprime; j++) {
                advance();
            }
            state -= coprime;
        }
        while (l < r) {
            advance();
        }
    }

    // Primes smaller or equal than N
    constexpr dynamic_bit_array sieve210(uint32_t N) {
        N++;
        dynamic_bit_array prime(to_ord(N));
        prime.set_all();
        static const auto [wheels, medium_primes_begin] = []() {
            constexpr size_t max_wheel_size = 1 << 20;
            uint32_t product = period * dynamic_bit_array::width / 4;
            big_vector<uint32_t> current;
            big_vector<wheel_t> wheels;
            for(size_t i = 0; i < size(sqrt_primes); i++) {
                uint32_t p = sqrt_primes[i];
                if (product * p > max_wheel_size) {
                    wheels.push_back(make_wheel(current, product));
                    current = {p};
                    product = period * dynamic_bit_array::width / 4 * p;
                    if (product > max_wheel_size) {
                        checkpoint("make wheels");
                        return std::pair{wheels, i};
                    }
                } else {
                    current.push_back(p);
                    product *= p;
                }
            }
            assert(false);
        }();
        static constexpr uint32_t dense_block = 1 << 25;
        for(uint32_t start = 0; start < N; start += dense_block) {
            uint32_t r = std::min(start + dense_block, N);  
            for(auto const& wheel: wheels) {
                auto l = start / wheel.product * wheel.product;
                sieve_dense(prime, to_ord(l), to_ord(r), wheel);
            }
        }
        checkpoint("dense sieve");
        static constexpr uint32_t sparse_block = 1 << 24;
        for(uint32_t start = 0; start < N; start += sparse_block) {
            uint32_t r = std::min(start + sparse_block, N);
            for(size_t i = medium_primes_begin; i < size(sqrt_primes); i++) {
                auto p = sqrt_primes[i];
                if(p * p >= r) break;
                auto k = std::max(start / p, p);
                if (!coprime210(k)) {k += add210[k % 210];}
                sieve210(prime, to_ord(p * k), to_ord(r), i, state210[k % 210]);
            }
        }
        checkpoint("sparse sieve");
        for(size_t i = 0; i < std::min(prime.words, sqrt_prime_bits.words); i++) {
            prime.word(i) = sqrt_prime_bits.word(i);
        }
        return prime;
    }
}
#pragma GCC pop_options
#endif // CP_ALGO_MATH_SIEVE_HPP
